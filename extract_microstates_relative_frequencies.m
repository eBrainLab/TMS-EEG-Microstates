function data_struct = extract_microstates_relative_frequencies(data_struct, microstates, baseline_indices, time_window_ranges, features2extract)
% Function to extract Relative Occurrence Frequencies (ROF) and
% Relative Transition Frequencies (RTF) from microstate data.

%% Relative Occurrence Frequencies (ROF)
if any(strcmp(features2extract, 'ROF'))
    fprintf('Starting Relative Occurrence Frequencies (ROF) computation...\n');
    % Iterate over each participant in the data structure
    for ii = 1:length(data_struct)
        % Extract raw microstate sequences for the current participant
        microstate_struct = data_struct(ii).microstate_data;

        % Get the field names
        fields = fieldnames(microstate_struct);

        for f = 1:length(fields)
            fieldName = fields{f};

            microstate_data = microstate_struct.(fieldName);

            % Concatenate raw data sequences into a single string array
            string_data = horzcat(microstate_data{:});

            % Determine the number of trials (columns in the concatenated data)
            ntrials = size(string_data, 2);

            % Process each microstate separately
            for ms_idx = 1:length(microstates)
                microstate = microstates(ms_idx);
                % Identify the occurrences of the current microstate across all trials
                index = string_data == microstate;

                % Calculate the occurrence proportion at each time point
                occurrence = sum(index, 2) / ntrials;

                % Store the occurrence data in the structure
                data_struct(ii).(fieldName).occurrences.(microstate) = occurrence;
            end

            % Combine all microstate occurrences into a matrix
            % (rows: time points, columns: microstates)
            all_occ = cell2mat(arrayfun(@(s) data_struct(ii).(fieldName).occurrences.(s), microstates, 'UniformOutput', false));

            % Initialize matrix for storing CLR-transformed occurrences
            clr_occ = zeros(size(all_occ));

            % Apply centered log-ratio (CLR) transformation for each time point
            delta = 0.5; % Small constant to avoid log(0)
            for jj = 1:size(all_occ, 1)
                clr_occ(jj, :) = centered_log_ratio(all_occ(jj, :), delta);
            end

            % Compute the median CLR values for the baseline period (for baseline correction)
            baseline_shift_clr = median(clr_occ(baseline_indices, :), 1);

            % Perform baseline correction on the CLR-transformed data
            baseline_corr_clr = clr_occ - repmat(baseline_shift_clr, size(clr_occ, 1), 1);

            % Store the CLR-transformed Baseline-corrected data
            for ms_idx = 1:length(microstates)
                data_struct(ii).(fieldName).occurrences_clr_bc.(microstates(ms_idx)) = baseline_corr_clr(:, ms_idx);
            end
        end
    end
    fprintf('Completed Relative Occurrence Frequencies (ROF) computation.\n');
end

%% Relative Transition Frequencies (RTF)
if any(strcmp(features2extract, 'RTF'))
    fprintf('Starting Relative Transition Frequencies (RTF) computation...\n');
    % Define state labels and mapping to numeric codes
    num_states = length(microstates);
    state_label_to_code = containers.Map(microstates, 1:num_states);

    % Define time vector and transition times
    time = linspace(-1000, 1000, 2001); % Time points for the data
    dt = time(2) - time(1); % Time step
    transition_times = time(1:end-1) + dt/2; % Times corresponding to transitions, length 2000

    % Compute time indices for the specified time windows
    window_names = fieldnames(time_window_ranges);
    time_window_indices = struct();
    for w = 1:length(window_names)
        window_name = window_names{w};
        time_range = time_window_ranges.(window_name); % [start_time, end_time]
        start_time = time_range(1);
        end_time = time_range(2);
        time_indices = (transition_times >= start_time) & (transition_times <= end_time);
        time_window_indices.(window_name) = time_indices;
    end

    % For each participant
    for i = 1:length(data_struct)
        % Extract raw microstate sequences for the current participant
        microstate_struct = data_struct(i).microstate_data;

        % Get the field names
        fields = fieldnames(microstate_struct);

        for f = 1:length(fields)
            fieldName = fields{f};

            % Initialize transitions_time_series
            transitions_time_series = zeros(num_states, num_states, length(transition_times)); % 5x5x2000
            num_trials = length(microstate_struct.(fieldName));

            % Process each trial
            for trial_idx = 1:num_trials
                state_labels_trial = data_struct(i).microstate_data.(fieldName){trial_idx}; % 2001x1 cell array
                % Convert to numeric codes
                state_numbers_trial = cellfun(@(x) state_label_to_code(x), state_labels_trial);
                % Get s1 and s2
                s1 = state_numbers_trial(1:end-1);
                s2 = state_numbers_trial(2:end);
                % Find transitions where s1 ~= s2
                transition_indices = find(s1 ~= s2);
                s1_transitions = s1(transition_indices);
                s2_transitions = s2(transition_indices);
                t_transitions = transition_indices;
                % Update transitions_time_series
                for k = 1:length(transition_indices)
                    s1_t = s1_transitions(k);
                    s2_t = s2_transitions(k);
                    t = t_transitions(k);
                    transitions_time_series(s1_t, s2_t, t) = transitions_time_series(s1_t, s2_t, t) + 1;
                end
            end
            % Average over trials
            transitions_time_series = transitions_time_series / num_trials;
            % Store transitions_time_series in data_struct
            data_struct(i).(fieldName).transitions_time_series = transitions_time_series;
            % Compute averages within specified time windows
            data_struct(i).(fieldName).transition_averages = struct();
            for w = 1:length(window_names)
                window_name = window_names{w};
                time_indices = time_window_indices.(window_name); % logical array
                % Initialize transition averages matrix
                transition_averages = zeros(num_states, num_states);
                for s1 = 1:num_states
                    for s2 = 1:num_states
                        if s1 ~= s2
                            % Get the data for this transition and time window
                            data = squeeze(transitions_time_series(s1, s2, time_indices));
                            % Compute average over time window
                            average_value = mean(data, 'omitnan');
                            transition_averages(s1, s2) = average_value;
                        end
                    end
                end
                % Store in data_struct
                data_struct(i).(fieldName).transition_averages.(window_name) = transition_averages;
            end
            % Baseline correction: subtract 'baseline' averages from other windows
            if isfield(data_struct(i).(fieldName).transition_averages, 'baseline')
                baseline_averages = data_struct(i).(fieldName).transition_averages.baseline;
                window_names_no_baseline = setdiff(window_names, 'baseline');
                for w = 1:length(window_names_no_baseline)
                    window_name = window_names_no_baseline{w};
                    data_struct(i).(fieldName).transition_averages_bc.(window_name) = ...
                        data_struct(i).(fieldName).transition_averages.(window_name) - baseline_averages;
                end
            else
                warning('Baseline window not found for subject %d. Skipping baseline correction.', i);
            end
        end
    end
    fprintf('Completed Relative Transition Frequencies (RTF) computation.\n');
end
% Display a message indicating the computation is complete
disp("Feature extraction complete.");
end
