function plot_microstate_occurrence(data_struct, current_stim_site, p_values_matrix, microstates, common_time_array, figure_time_range, figure_title, save_figure)
% PLOT_MICROSTATE_OCCURRENCE Visualizes microstate occurrence with confidence intervals and significant time markers.
%
% This function generates a plot displaying the relative occurrence frequency
% of microstates over time, including 95% confidence intervals for the averages
% and markers indicating time points with statistically significant effects.
%
% Inputs:
%   - data_struct: A struct array containing the occurrence data for all subjects.
%     Each subject's data should have fields for each stimulation site and microstate.
%   - current_stim_site: A string (or array of two strings) representing the
%     stimulation site(s) for which the data is being plotted. If two sites are
%     provided, differences between the two are computed.
%   - p_values_matrix: A matrix (num_states x num_timepoints) of p-values for
%     statistical tests associated with the microstates at each time point.
%   - microstates: A cell array of strings, each representing the name of a microstate.
%   - common_time_array: A vector of time points (in milliseconds) common to all
%     subjects and conditions.
%   - figure_time_range: A two-element vector specifying the range of time
%     (e.g., [-200, 500]) to display on the plot, in milliseconds.
%   - figure_title: A string specifying the title of the figure.
%   - save_figure: A logical value (true/false). If true, the figure will be
%     saved as a PNG file with a resolution of 500 dpi.
%
% Outputs:
%   The function generates a figure displaying:
%     1. Smoothed relative occurrence frequency of each microstate.
%     2. 95% confidence intervals as shaded regions.
%     3. Horizontal bars indicating time points with significant effects
%        (above the plot for positive effects, below for negative effects).
%
% Notes:
%   - The function automatically excludes a range of time points near
%     [-10, 20 ms] from the analysis and visualization.
%   - The figure is customized with Helvetica font, rotated x-axis labels,
%     and specific grid exclusions.
%   - The function now handles missing data for specific stimulation sites.

% Determine the indices for the desired time range (-200 to 500 ms)
time_range_start = figure_time_range(1);
time_range_end = figure_time_range(2);
time_indices = find(common_time_array >= time_range_start & common_time_array <= time_range_end);

% Ensure current_stim_site is properly handled
if ischar(current_stim_site) || isstring(current_stim_site)
    % If it's a string, keep it as is
elseif ~iscell(current_stim_site)
    % If it's not a cell or string, make it a cell
    current_stim_site = {current_stim_site};
end

num_subjects = length(data_struct);
num_timepoints = length(common_time_array);
num_states = length(microstates);

% Preallocate matrices for averages and p-values
average_matrix = zeros(num_states, num_timepoints);

% Initialize arrays to store the data of each field for each subject
data_all_subjects = struct();
for i = 1:num_states
    current_state = microstates{i};
    data_all_subjects.(current_state) = NaN(2001, num_subjects); % Use NaN to identify missing data
end

% First, identify which subjects have the required data
valid_subjects = [];
for sub = 1:num_subjects
    if iscell(current_stim_site)
        % We have multiple stimulation sites (cell array)
        if numel(current_stim_site) > 1
            % Check if both stimulation sites exist for this subject
            site1 = current_stim_site{1};
            site2 = current_stim_site{2};
            state1 = microstates{1};
            
            if isfield(data_struct(sub), site1) && isfield(data_struct(sub), site2)
                % Check if the first microstate exists in both sites
                site1_data = data_struct(sub).(site1);
                site2_data = data_struct(sub).(site2);
                
                if isfield(site1_data, 'occurrences_clr_bc') && isfield(site2_data, 'occurrences_clr_bc')
                    if isfield(site1_data.occurrences_clr_bc, state1) && ...
                       isfield(site2_data.occurrences_clr_bc, state1)
                        valid_subjects = [valid_subjects, sub];
                    end
                end
            end
        end
    else
        % Single stimulation site (string)
        if isfield(data_struct(sub), current_stim_site) 
            % Check if the occurrences_clr_bc field exists
            site_data = data_struct(sub).(current_stim_site);
            if isfield(site_data, 'occurrences_clr_bc')
                % Check if the first microstate exists
                if isfield(site_data.occurrences_clr_bc, microstates{1})
                    valid_subjects = [valid_subjects, sub];
                end
            end
        end
    end
end

num_valid_subjects = length(valid_subjects);
if num_valid_subjects == 0
    error('No valid subjects found with required stimulation sites and data.');
end

fprintf('Found %d valid subjects out of %d total subjects.\n', num_valid_subjects, num_subjects);

% Preallocate data_microstate for each iteration to avoid dimension errors
data_microstate = NaN(num_valid_subjects, num_timepoints);

% Loop through each valid subject and store the data
for idx = 1:num_valid_subjects
    sub = valid_subjects(idx);
    
    for i = 1:num_states
        current_state = microstates{i};
        
        if iscell(current_stim_site) && numel(current_stim_site) > 1
            % Using two stimulation sites (difference calculation)
            site1 = current_stim_site{1};
            site2 = current_stim_site{2};
            
            site1_data = data_struct(sub).(site1).occurrences_clr_bc.(current_state);
            site2_data = data_struct(sub).(site2).occurrences_clr_bc.(current_state);
            
            data_all_subjects.(current_state)(:, idx) = site1_data - site2_data;
        else
            % Using a single stimulation site
            if iscell(current_stim_site)
                site = current_stim_site{1}; % If it's a cell but with only one element
            else
                site = current_stim_site;    % If it's a string
            end
            
            data_all_subjects.(current_state)(:, idx) = ...
                data_struct(sub).(site).occurrences_clr_bc.(current_state);
        end
    end
end

for i = 1:num_states
    % Reset data_microstate for each microstate
    data_microstate = NaN(num_valid_subjects, num_timepoints);
    current_state = microstates{i};
    
    for idx = 1:num_valid_subjects
        sub = valid_subjects(idx);
        
        if iscell(current_stim_site) && numel(current_stim_site) > 1
            % Using two stimulation sites (difference calculation)
            site1 = current_stim_site{1};
            site2 = current_stim_site{2};
            
            site1_data = data_struct(sub).(site1).occurrences_clr_bc.(current_state);
            site2_data = data_struct(sub).(site2).occurrences_clr_bc.(current_state);
            
            data_microstate(idx, :) = (site1_data - site2_data)';
        else
            % Using a single stimulation site
            if iscell(current_stim_site)
                site = current_stim_site{1}; % If it's a cell but with only one element
            else
                site = current_stim_site;    % If it's a string
            end
            
            data_microstate(idx, :) = data_struct(sub).(site).occurrences_clr_bc.(current_state)';
        end
    end
    
    % Calculate average across valid subjects
    average_state = mean(data_microstate, 1, 'omitnan');
    average_matrix(i, :) = average_state;
end

% Calculate the average and standard deviation for each field, using only valid data
averages = struct();
std_devs = struct();
for i = 1:num_states
    current_state = microstates{i};
    averages.(current_state) = mean(data_all_subjects.(current_state), 2, 'omitnan'); % Average over subjects
    std_devs.(current_state) = std(data_all_subjects.(current_state), 0, 2, 'omitnan'); % Standard deviation over subjects
end

% Calculate the standard error of the mean (SEM)
standard_errors = struct();
for i = 1:num_states
    current_state = microstates{i};
    % Count valid (non-NaN) data points for each timepoint
    valid_count = sum(~isnan(data_all_subjects.(current_state)), 2);
    % Avoid division by zero
    valid_count(valid_count == 0) = NaN;
    standard_errors.(current_state) = std_devs.(current_state) ./ sqrt(valid_count);
end

% Calculate the critical t-value for 95% confidence interval
alpha = 0.05;
degrees_of_freedom = num_valid_subjects - 1;
t_critical = tinv(1 - alpha/2, degrees_of_freedom); % Two-tailed test

% Calculate the 95% confidence intervals
confidence_intervals = struct();
for i = 1:num_states
    current_state = microstates{i};
    confidence_intervals.(current_state) = t_critical * standard_errors.(current_state);
end

% Smooth the averaged data and the confidence intervals with less smoothing
smoothed_averages = struct();
smoothed_conf_intervals = struct();
for i = 1:num_states
    current_state = microstates{i};
    smoothed_averages.(current_state) = smooth(averages.(current_state), 0.02, 'loess'); % Using less loess smoothing
    smoothed_conf_intervals.(current_state) = smooth(confidence_intervals.(current_state), 0.02, 'loess'); % Smooth CI as well
end

% Define the exclusion range
exclude_start = -10;
exclude_end = 20;

% Find the indices within 'common_time_array' that fall into this range
exclude_indices = find(common_time_array >= exclude_start & common_time_array <= exclude_end);

% Set the average and confidence interval data to NaN for these indices
for i = 1:num_states
    current_state = microstates{i};
    smoothed_averages.(current_state)(exclude_indices) = NaN;
    smoothed_conf_intervals.(current_state)(exclude_indices) = NaN;
end

% Plot the smoothed averaged data for the specified time range with thicker lines
figure;
set(gcf, 'Position', [100, 100, 2000, 1000]); % [left, bottom, width, height]
hold on;

% Get distinct colors for each state
colors = lines(num_states);

% First, plot the shaded areas representing the 95% confidence intervals
for i = 1:num_states
    current_state = microstates{i};

    % Get the smoothed average and confidence interval for the time indices
    avg = smoothed_averages.(current_state)(time_indices)';
    conf_interval = smoothed_conf_intervals.(current_state)(time_indices)';

    % Define the upper and lower bounds of the shaded area
    upper = avg + conf_interval;
    lower = avg - conf_interval;

    % Identify non-NaN regions
    pre_indices = find(common_time_array(time_indices) < exclude_start);
    post_indices = find(common_time_array(time_indices) > exclude_end);

    % Plot shaded area for pre-exclusion range
    if ~isempty(pre_indices)
        non_nan_indices = ~isnan(avg(pre_indices)) & ~isnan(upper(pre_indices)) & ~isnan(lower(pre_indices));
        if any(non_nan_indices)
            x_pre = [common_time_array(time_indices(pre_indices(non_nan_indices))), fliplr(common_time_array(time_indices(pre_indices(non_nan_indices))))];
            y_pre = [upper(pre_indices(non_nan_indices)), fliplr(lower(pre_indices(non_nan_indices)))];
            fill(x_pre, y_pre, colors(i, :), 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        end
    end

    % Plot shaded area for post-exclusion range
    if ~isempty(post_indices)
        non_nan_indices = ~isnan(avg(post_indices)) & ~isnan(upper(post_indices)) & ~isnan(lower(post_indices));
        if any(non_nan_indices)
            x_post = [common_time_array(time_indices(post_indices(non_nan_indices))), fliplr(common_time_array(time_indices(post_indices(non_nan_indices))))];
            y_post = [upper(post_indices(non_nan_indices)), fliplr(lower(post_indices(non_nan_indices)))];
            fill(x_post, y_post, colors(i, :), 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        end
    end
end

% Then, plot the curves
for i = 1:num_states
    current_state = microstates{i};
    plot(common_time_array(time_indices), smoothed_averages.(current_state)(time_indices), ...
        'Color', colors(i, :), 'DisplayName', current_state, 'LineWidth', 4);
end

% Calculate the global minimum and maximum values across all states
global_min = inf;
global_max = -inf;

for i = 1:num_states
    current_state = microstates{i};
    data_values = smoothed_averages.(current_state)(time_indices);
    global_min = min(global_min, min(data_values, [], 'omitnan'));
    global_max = max(global_max, max(data_values, [], 'omitnan'));
end

% Identify states with any p-values < 0.05
states_with_effects = find(any(p_values_matrix < 0.05, 2));

% Number of states with significant effects
num_sig_states = length(states_with_effects);

% Initialize lists for positive and negative significant states
states_with_positive_effects = [];
states_with_negative_effects = [];

for idx = 1:num_sig_states
    s = states_with_effects(idx);
    % Calculate the mean effect during significant time points
    sig_indices = find(p_values_matrix(s, :) < 0.05);
    mean_effect = mean(average_matrix(s, sig_indices), 'omitnan');

    if mean_effect > 0
        states_with_positive_effects(end+1) = s;
    else
        states_with_negative_effects(end+1) = s;
    end
end

% Number of positive and negative significant states
num_pos_sig_states = length(states_with_positive_effects);
num_neg_sig_states = length(states_with_negative_effects);

% Define the vertical spacing for lines above and below the plot
bar_height = 0.2;  % Adjust as needed for spacing

% Assign y-positions for positive and negative significant lines
pos_y_positions = global_max + (1:num_pos_sig_states) * bar_height;
neg_y_positions = global_min - (1:num_neg_sig_states) * bar_height;

% Calculate required ylim
ymin = global_min - (num_neg_sig_states + 1) * bar_height;
ymax = global_max + (num_pos_sig_states + 1) * bar_height;

% Adjust the ylim to include data and lines
% ylim([ymin, ymax]);
ylim([-3, 3]);
xlim([time_range_start, time_range_end]);

% Set x-axis ticks to show labels every 50 ms
xticks(time_range_start:50:time_range_end);

% Plot horizontal lines for positive significant effects
for idx = 1:num_pos_sig_states
    s = states_with_positive_effects(idx);
    current_state = microstates{s};

    % Find significant p-values for this state
    significant = p_values_matrix(s, :) < 0.05;

    % Find indices where p < 0.05
    sig_indices = find(significant);

    if isempty(sig_indices)
        continue;  % Skip if no significant points found
    end

    % Identify the start and end of each continuous significant segment
    d = diff(sig_indices);
    split_points = find(d > 1);

    % Initialize start and end indices for segments
    segment_starts = [sig_indices(1), sig_indices(split_points + 1)];
    segment_ends = [sig_indices(split_points), sig_indices(end)];

    % Assign y-position for significant lines above the plot
    y_pos = pos_y_positions(idx);

    % Plot horizontal lines for each significant segment
    for seg = 1:length(segment_starts)
        start_time = common_time_array(segment_starts(seg));
        end_time = common_time_array(segment_ends(seg));

        % Plot the horizontal line
        line([start_time, end_time], [y_pos, y_pos], 'Color', colors(s, :), 'LineWidth', 8);
    end
end

% Plot horizontal lines for negative significant effects
for idx = 1:num_neg_sig_states
    s = states_with_negative_effects(idx);
    current_state = microstates{s};

    % Find significant p-values for this state
    significant = p_values_matrix(s, :) < 0.05;

    % Find indices where p < 0.05
    sig_indices = find(significant);

    if isempty(sig_indices)
        continue;  % Skip if no significant points found
    end

    % Identify the start and end of each continuous significant segment
    d = diff(sig_indices);
    split_points = find(d > 1);

    % Initialize start and end indices for segments
    segment_starts = [sig_indices(1), sig_indices(split_points + 1)];
    segment_ends = [sig_indices(split_points), sig_indices(end)];

    % Assign y-position for significant lines below the plot
    y_pos = neg_y_positions(idx);

    % Plot horizontal lines for each significant segment
    for seg = 1:length(segment_starts)
        start_time = common_time_array(segment_starts(seg));
        end_time = common_time_array(segment_ends(seg));

        % Plot the horizontal line
        line([start_time, end_time], [y_pos, y_pos], 'Color', colors(s, :), 'LineWidth', 8);
    end
end

% Adjust font properties
ax = gca;
ax.FontSize = 30;
ax.FontName = 'Helvetica';
ax.XTickLabelRotation = 0;

% Exclude specific grid lines
grid on;
xGridLines = ax.XTick;
yGridLines = ax.YTick;
xGridLines = setdiff(xGridLines, 600); % Exclude x=600
yGridLines = setdiff(yGridLines, 3);   % Exclude y=3
grid off;
for k = 1:length(xGridLines)
    xline(xGridLines(k), '-', 'Color', [0.8, 0.8, 0.8]);
end
for k = 1:length(yGridLines)
    yline(yGridLines(k), '-', 'Color', [0.8, 0.8, 0.8]);
end

% Label the axes and add a title
xlabel('Latency (ms)', 'FontSize', 35, 'FontName', 'Helvetica');
ylabel('Relative Occurrence Frequency', 'FontSize', 35, 'FontName', 'Helvetica');
title(figure_title, 'FontSize', 40, 'FontName', 'Helvetica');

% Add reference lines
yline(0, '--k');
xline(0, '--r', 'TMS', 'LineWidth', 5, ...
    'LabelHorizontalAlignment', 'center', ...
    'LabelVerticalAlignment', 'top', ...
    'LabelOrientation','horizontal', ...
    'FontSize', 35, 'FontName', 'Helvetica');

% Create the original microstate legend
lgd_microstate = legend(microstates, 'FontSize', 28, 'FontName', 'Helvetica');
title(lgd_microstate, 'Microstate');
lgd_microstate.Location = 'northeast';

hold off;

% Save the figure as a TIF file with specified resolution
if save_figure
    safe_figure_title = matlab.lang.makeValidName(figure_title);
    % exportgraphics(gcf, [safe_figure_title, '.png'], 'Resolution', 500);
    exportgraphics(gcf, [safe_figure_title, '.pdf'], 'ContentType', 'vector');
    close all;
end
end