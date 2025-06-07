function data_struct = extract_microstates_frequencies(data_struct, microstates, ...
    baseline_indices, time_window_ranges, features2extract)
% -------------------------------------------------------------------------
% EXTRACT_MICROSTATES_FREQUENCIES - Extract occurrence and transition frequencies
%
% This function extracts two key measures from microstate time series data:
%   1. Relative Occurrence Frequencies (ROF) - baseline-corrected CLR values
%   2. Relative Transition Frequencies (RTF) - transition probabilities
%
% Syntax:
%   data_struct = extract_microstates_frequencies(data_struct, microstates, ...
%       baseline_indices, time_window_ranges, features2extract)
%
% Inputs:
%   data_struct        - Structure array containing microstate data for each subject
%   microstates        - Cell array of microstate labels (e.g., ["A","B","C","D","E"])
%   baseline_indices   - Indices corresponding to the baseline period
%   time_window_ranges - Structure with time window definitions
%                        (e.g., struct('baseline', [-1000,-10], 'post_tms', [20,500]))
%   features2extract   - Cell array specifying features ('ROF', 'RTF', or both)
%
% Outputs:
%   data_struct - Updated structure with extracted features added to each subject
%
% Example:
%   data_struct = extract_microstates_frequencies(data_struct, ["A","B","C","D","E"], ...
%       1:991, time_windows, {'ROF', 'RTF'});
%
% See also: CENTERED_LOG_RATIO, TEST_ROF_ONE_CONDITION, TEST_RTF
% -------------------------------------------------------------------------

%% Extract Relative Occurrence Frequencies (ROF)
if any(strcmp(features2extract, 'ROF'))
    fprintf('Starting Relative Occurrence Frequencies computation (ROF)...\n');
    
    % Process each participant
    for ii = 1:length(data_struct)
        % Extract microstate data for current participant
        microstate_struct = data_struct(ii).microstate_data;
        fields = fieldnames(microstate_struct);
        
        % Process each field (stimulation site)
        for f = 1:length(fields)
            fieldName = fields{f};
            microstate_data = microstate_struct.(fieldName);
            
            % Concatenate all trials into single array
            string_data = horzcat(microstate_data{:});
            ntrials = size(string_data, 2);
            
            % Calculate occurrence proportion for each microstate
            for ms_idx = 1:length(microstates)
                microstate = microstates(ms_idx);
                
                % Find occurrences of current microstate
                index = string_data == microstate;
                
                % Calculate proportion at each time point
                occurrence = sum(index, 2) / ntrials;
                
                % Store occurrence data
                data_struct(ii).(fieldName).occurrences.(microstate) = occurrence;
            end
            
            % Combine all microstate occurrences into matrix
            all_occ = cell2mat(arrayfun(@(s) ...
                data_struct(ii).(fieldName).occurrences.(s), ...
                microstates, 'UniformOutput', false));
            
            % Apply CLR transformation
            clr_occ = zeros(size(all_occ));
            delta = 0.5;  % Multiplicative replacement parameter
            
            for jj = 1:size(all_occ, 1)
                clr_occ(jj, :) = centered_log_ratio(all_occ(jj, :), delta);
            end
            
            % Baseline correction
            baseline_shift_clr = median(clr_occ(baseline_indices, :), 1);
            baseline_corr_clr = clr_occ - repmat(baseline_shift_clr, size(clr_occ, 1), 1);
            
            % Store baseline-corrected CLR data
            for ms_idx = 1:length(microstates)
                data_struct(ii).(fieldName).occurrences_clr_bc.(microstates(ms_idx)) = ...
                    baseline_corr_clr(:, ms_idx);
            end
        end
    end
    fprintf('Completed Relative Occurrence Frequencies computation (ROF).\n');
end

%% Extract Relative Transition Frequencies (RTF)
if any(strcmp(features2extract, 'RTF'))
    fprintf('Starting Relative Transition Frequencies computation (RTF)...\n');
    
    % Initialize state mapping
    num_states = length(microstates);
    state_label_to_code = containers.Map(microstates, 1:num_states);
    
    % Define time vectors
    time = linspace(-1000, 1000, 2001);  % Time points
    dt = time(2) - time(1);              % Time step
    transition_times = time(1:end-1) + dt/2;  % Transition midpoints
    
    % Compute time window indices
    window_names = fieldnames(time_window_ranges);
    time_window_indices = struct();
    
    for w = 1:length(window_names)
        window_name = window_names{w};
        time_range = time_window_ranges.(window_name);
        start_time = time_range(1);
        end_time = time_range(2);
        time_indices = (transition_times >= start_time) & (transition_times <= end_time);
        time_window_indices.(window_name) = time_indices;
    end
    
    % Process each participant
    for i = 1:length(data_struct)
        microstate_struct = data_struct(i).microstate_data;
        fields = fieldnames(microstate_struct);
        
        % Process each field (stimulation site)
        for f = 1:length(fields)
            fieldName = fields{f};
            
            % Initialize transition time series
            transitions_time_series = zeros(num_states, num_states, length(transition_times));
            num_trials = length(microstate_struct.(fieldName));
            
            % Process each trial
            for trial_idx = 1:num_trials
                % Get state sequence for this trial
                state_labels_trial = data_struct(i).microstate_data.(fieldName){trial_idx};
                
                % Convert to numeric codes
                state_numbers_trial = cellfun(@(x) state_label_to_code(x), state_labels_trial);
                
                % Find transitions (where consecutive states differ)
                s1 = state_numbers_trial(1:end-1);
                s2 = state_numbers_trial(2:end);
                transition_indices = find(s1 ~= s2);
                
                % Update transition counts
                for k = 1:length(transition_indices)
                    from_state = s1(transition_indices(k));
                    to_state = s2(transition_indices(k));
                    time_idx = transition_indices(k);
                    transitions_time_series(from_state, to_state, time_idx) = ...
                        transitions_time_series(from_state, to_state, time_idx) + 1;
                end
            end
            
            % Average over trials
            transitions_time_series = transitions_time_series / num_trials;
            data_struct(i).(fieldName).transitions_time_series = transitions_time_series;
            
            % Calculate average transitions for each time window
            data_struct(i).(fieldName).transition_averages = struct();
            
            for w = 1:length(window_names)
                window_name = window_names{w};
                time_indices = time_window_indices.(window_name);
                transition_averages = zeros(num_states, num_states);
                
                % Average transitions within window (excluding self-transitions)
                for s1 = 1:num_states
                    for s2 = 1:num_states
                        if s1 ~= s2
                            data = squeeze(transitions_time_series(s1, s2, time_indices));
                            transition_averages(s1, s2) = mean(data, 'omitnan');
                        end
                    end
                end
                
                data_struct(i).(fieldName).transition_averages.(window_name) = transition_averages;
            end
            
            % Baseline correction
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
    fprintf('Completed Relative Transition Frequencies computation (RTF).\n');
end

%% Display completion message
features_extracted = [];
if any(strcmp(features2extract, 'ROF')), features_extracted = [features_extracted, "ROF"]; end
if any(strcmp(features2extract, 'RTF')), features_extracted = [features_extracted, "RTF"]; end

fprintf("Feature extraction complete. Extracted: %s\n", strjoin(features_extracted, ', '));

end