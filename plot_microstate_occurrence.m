function plot_microstate_occurrence(data_struct, current_stim_site, p_values_matrix, ...
    microstates, common_time_array, figure_time_range, figure_title, save_figure)
% -------------------------------------------------------------------------
% PLOT_MICROSTATE_OCCURRENCE - Visualize microstate occurrence with statistics
%
% This function creates a publication-ready plot showing microstate relative
% occurrence frequencies over time, including confidence intervals and
% significance markers.
%
% Syntax:
%   plot_microstate_occurrence(data_struct, current_stim_site, p_values_matrix, ...
%       microstates, common_time_array, figure_time_range, figure_title, save_figure)
%
% Inputs:
%   data_struct        - Structure array with occurrence data for all subjects
%   current_stim_site  - String or cell array of stimulation site(s)
%                        If two sites provided, plots difference
%   p_values_matrix    - Matrix [states x timepoints] of p-values
%   microstates        - Cell array of microstate labels
%   common_time_array  - Vector of time points (ms)
%   figure_time_range  - Two-element vector [start, end] for x-axis (ms)
%   figure_title       - String for figure title
%   save_figure        - Boolean to save figure as PDF
%
% Features:
%   - Smoothed occurrence curves with 95% confidence intervals
%   - Horizontal bars indicating significant time periods
%   - Automatic exclusion of artifact period [-10, 20 ms]
%   - Handles missing data gracefully
%
% Example:
%   plot_microstate_occurrence(data_struct, 'lpfc', p_vals, ...
%       ["A","B","C","D","E"], -1000:1000, [-200,600], 'LPFC Stimulation', true);
%
% See also: EXTRACT_MICROSTATES_FREQUENCIES, TEST_ROF_ONE_CONDITION
% -------------------------------------------------------------------------

%% Initialize parameters
% Define time range indices
time_range_start = figure_time_range(1);
time_range_end = figure_time_range(2);
time_indices = find(common_time_array >= time_range_start & ...
                   common_time_array <= time_range_end);

% Ensure stimulation site is properly formatted
if ischar(current_stim_site) || isstring(current_stim_site)
    % Keep as is
elseif ~iscell(current_stim_site)
    current_stim_site = {current_stim_site};
end

% Get dimensions
num_subjects = length(data_struct);
num_timepoints = length(common_time_array);
num_states = length(microstates);

% Preallocate
average_matrix = zeros(num_states, num_timepoints);

% Initialize data storage
data_all_subjects = struct();
for i = 1:num_states
    current_state = microstates{i};
    data_all_subjects.(current_state) = NaN(2001, num_subjects);
end

%% Identify valid subjects
valid_subjects = [];
for sub = 1:num_subjects
    if iscell(current_stim_site) && numel(current_stim_site) > 1
        % Check for both stimulation sites
        site1 = current_stim_site{1};
        site2 = current_stim_site{2};
        state1 = microstates{1};
        
        if isfield(data_struct(sub), site1) && isfield(data_struct(sub), site2)
            site1_data = data_struct(sub).(site1);
            site2_data = data_struct(sub).(site2);
            
            if isfield(site1_data, 'occurrences_clr_bc') && ...
               isfield(site2_data, 'occurrences_clr_bc') && ...
               isfield(site1_data.occurrences_clr_bc, state1) && ...
               isfield(site2_data.occurrences_clr_bc, state1)
                valid_subjects = [valid_subjects, sub];
            end
        end
    else
        % Single stimulation site
        if isfield(data_struct(sub), current_stim_site)
            site_data = data_struct(sub).(current_stim_site);
            if isfield(site_data, 'occurrences_clr_bc') && ...
               isfield(site_data.occurrences_clr_bc, microstates{1})
                valid_subjects = [valid_subjects, sub];
            end
        end
    end
end

num_valid_subjects = length(valid_subjects);
if num_valid_subjects == 0
    error('No valid subjects found with required stimulation sites and data.');
end

fprintf('Found %d valid subjects out of %d total subjects.\n', ...
    num_valid_subjects, num_subjects);

%% Extract and process data
% Store data for all valid subjects
for idx = 1:num_valid_subjects
    sub = valid_subjects(idx);
    
    for i = 1:num_states
        current_state = microstates{i};
        
        if iscell(current_stim_site) && numel(current_stim_site) > 1
            % Calculate difference between sites
            site1 = current_stim_site{1};
            site2 = current_stim_site{2};
            
            site1_data = data_struct(sub).(site1).occurrences_clr_bc.(current_state);
            site2_data = data_struct(sub).(site2).occurrences_clr_bc.(current_state);
            
            data_all_subjects.(current_state)(:, idx) = site1_data - site2_data;
        else
            % Single site data
            if iscell(current_stim_site)
                site = current_stim_site{1};
            else
                site = current_stim_site;
            end
            
            data_all_subjects.(current_state)(:, idx) = ...
                data_struct(sub).(site).occurrences_clr_bc.(current_state);
        end
    end
end

% Calculate averages for plotting
for i = 1:num_states
    current_state = microstates{i};
    data_microstate = NaN(num_valid_subjects, num_timepoints);
    
    for idx = 1:num_valid_subjects
        sub = valid_subjects(idx);
        
        if iscell(current_stim_site) && numel(current_stim_site) > 1
            % Difference calculation
            site1 = current_stim_site{1};
            site2 = current_stim_site{2};
            
            site1_data = data_struct(sub).(site1).occurrences_clr_bc.(current_state);
            site2_data = data_struct(sub).(site2).occurrences_clr_bc.(current_state);
            
            data_microstate(idx, :) = (site1_data - site2_data)';
        else
            % Single site
            if iscell(current_stim_site)
                site = current_stim_site{1};
            else
                site = current_stim_site;
            end
            
            data_microstate(idx, :) = data_struct(sub).(site).occurrences_clr_bc.(current_state)';
        end
    end
    
    % Calculate average across subjects
    average_state = mean(data_microstate, 1, 'omitnan');
    average_matrix(i, :) = average_state;
end

%% Calculate statistics
% Calculate averages and standard deviations
averages = struct();
std_devs = struct();
for i = 1:num_states
    current_state = microstates{i};
    averages.(current_state) = mean(data_all_subjects.(current_state), 2, 'omitnan');
    std_devs.(current_state) = std(data_all_subjects.(current_state), 0, 2, 'omitnan');
end

% Calculate standard errors
standard_errors = struct();
for i = 1:num_states
    current_state = microstates{i};
    valid_count = sum(~isnan(data_all_subjects.(current_state)), 2);
    valid_count(valid_count == 0) = NaN;
    standard_errors.(current_state) = std_devs.(current_state) ./ sqrt(valid_count);
end

% Calculate 95% confidence intervals
alpha = 0.05;
degrees_of_freedom = num_valid_subjects - 1;
t_critical = tinv(1 - alpha/2, degrees_of_freedom);

confidence_intervals = struct();
for i = 1:num_states
    current_state = microstates{i};
    confidence_intervals.(current_state) = t_critical * standard_errors.(current_state);
end

%% Smooth data
smoothed_averages = struct();
smoothed_conf_intervals = struct();
for i = 1:num_states
    current_state = microstates{i};
    smoothed_averages.(current_state) = smooth(averages.(current_state), 0.02, 'loess');
    smoothed_conf_intervals.(current_state) = smooth(confidence_intervals.(current_state), 0.02, 'loess');
end

% Exclude artifact period
exclude_start = -10;
exclude_end = 20;
exclude_indices = find(common_time_array >= exclude_start & common_time_array <= exclude_end);

for i = 1:num_states
    current_state = microstates{i};
    smoothed_averages.(current_state)(exclude_indices) = NaN;
    smoothed_conf_intervals.(current_state)(exclude_indices) = NaN;
end

%% Create figure
figure;
set(gcf, 'Position', [100, 100, 2000, 1000]);
hold on;

% Get distinct colors
colors = lines(num_states);

% Plot confidence intervals (shaded areas)
for i = 1:num_states
    current_state = microstates{i};
    
    % Get data for time range
    avg = smoothed_averages.(current_state)(time_indices)';
    conf_interval = smoothed_conf_intervals.(current_state)(time_indices)';
    
    % Define bounds
    upper = avg + conf_interval;
    lower = avg - conf_interval;
    
    % Split into pre- and post-exclusion regions
    pre_indices = find(common_time_array(time_indices) < exclude_start);
    post_indices = find(common_time_array(time_indices) > exclude_end);
    
    % Plot pre-exclusion shaded area
    if ~isempty(pre_indices)
        non_nan_indices = ~isnan(avg(pre_indices)) & ~isnan(upper(pre_indices)) & ~isnan(lower(pre_indices));
        if any(non_nan_indices)
            x_pre = [common_time_array(time_indices(pre_indices(non_nan_indices))), ...
                     fliplr(common_time_array(time_indices(pre_indices(non_nan_indices))))];
            y_pre = [upper(pre_indices(non_nan_indices)), ...
                     fliplr(lower(pre_indices(non_nan_indices)))];
            fill(x_pre, y_pre, colors(i, :), 'FaceAlpha', 0.3, 'EdgeColor', 'none', ...
                 'HandleVisibility', 'off');
        end
    end
    
    % Plot post-exclusion shaded area
    if ~isempty(post_indices)
        non_nan_indices = ~isnan(avg(post_indices)) & ~isnan(upper(post_indices)) & ~isnan(lower(post_indices));
        if any(non_nan_indices)
            x_post = [common_time_array(time_indices(post_indices(non_nan_indices))), ...
                      fliplr(common_time_array(time_indices(post_indices(non_nan_indices))))];
            y_post = [upper(post_indices(non_nan_indices)), ...
                      fliplr(lower(post_indices(non_nan_indices)))];
            fill(x_post, y_post, colors(i, :), 'FaceAlpha', 0.3, 'EdgeColor', 'none', ...
                 'HandleVisibility', 'off');
        end
    end
end

% Plot main curves
for i = 1:num_states
    current_state = microstates{i};
    plot(common_time_array(time_indices), smoothed_averages.(current_state)(time_indices), ...
        'Color', colors(i, :), 'DisplayName', current_state, 'LineWidth', 4);
end

%% Add significance markers
% Find data range
global_min = inf;
global_max = -inf;

for i = 1:num_states
    current_state = microstates{i};
    data_values = smoothed_averages.(current_state)(time_indices);
    global_min = min(global_min, min(data_values, [], 'omitnan'));
    global_max = max(global_max, max(data_values, [], 'omitnan'));
end

% Identify significant states
states_with_effects = find(any(p_values_matrix < 0.05, 2));
num_sig_states = length(states_with_effects);

% Separate positive and negative effects
states_with_positive_effects = [];
states_with_negative_effects = [];

for idx = 1:num_sig_states
    s = states_with_effects(idx);
    sig_indices = find(p_values_matrix(s, :) < 0.05);
    mean_effect = mean(average_matrix(s, sig_indices), 'omitnan');
    
    if mean_effect > 0
        states_with_positive_effects(end+1) = s;
    else
        states_with_negative_effects(end+1) = s;
    end
end

% Set positions for significance bars
bar_height = 0.2;
num_pos_sig_states = length(states_with_positive_effects);
num_neg_sig_states = length(states_with_negative_effects);

pos_y_positions = global_max + (1:num_pos_sig_states) * bar_height;
neg_y_positions = global_min - (1:num_neg_sig_states) * bar_height;

% Plot positive significance bars
for idx = 1:num_pos_sig_states
    s = states_with_positive_effects(idx);
    significant = p_values_matrix(s, :) < 0.05;
    sig_indices = find(significant);
    
    if isempty(sig_indices), continue; end
    
    % Find continuous segments
    d = diff(sig_indices);
    split_points = find(d > 1);
    segment_starts = [sig_indices(1), sig_indices(split_points + 1)];
    segment_ends = [sig_indices(split_points), sig_indices(end)];
    
    y_pos = pos_y_positions(idx);
    
    % Plot each segment
    for seg = 1:length(segment_starts)
        start_time = common_time_array(segment_starts(seg));
        end_time = common_time_array(segment_ends(seg));
        line([start_time, end_time], [y_pos, y_pos], 'Color', colors(s, :), 'LineWidth', 8);
    end
end

% Plot negative significance bars
for idx = 1:num_neg_sig_states
    s = states_with_negative_effects(idx);
    significant = p_values_matrix(s, :) < 0.05;
    sig_indices = find(significant);
    
    if isempty(sig_indices), continue; end
    
    % Find continuous segments
    d = diff(sig_indices);
    split_points = find(d > 1);
    segment_starts = [sig_indices(1), sig_indices(split_points + 1)];
    segment_ends = [sig_indices(split_points), sig_indices(end)];
    
    y_pos = neg_y_positions(idx);
    
    % Plot each segment
    for seg = 1:length(segment_starts)
        start_time = common_time_array(segment_starts(seg));
        end_time = common_time_array(segment_ends(seg));
        line([start_time, end_time], [y_pos, y_pos], 'Color', colors(s, :), 'LineWidth', 8);
    end
end

%% Format figure
% Set axis properties
ylim([-3, 3]);
xlim([time_range_start, time_range_end]);
xticks(time_range_start:50:time_range_end);

% Font settings
ax = gca;
ax.FontSize = 20;
ax.FontName = 'Helvetica';
ax.XTickLabelRotation = 0;

% Custom grid
grid on;
xGridLines = ax.XTick;
yGridLines = ax.YTick;
xGridLines = setdiff(xGridLines, 600);
yGridLines = setdiff(yGridLines, 3);
grid off;

for k = 1:length(xGridLines)
    xline(xGridLines(k), '-', 'Color', [0.8, 0.8, 0.8]);
end
for k = 1:length(yGridLines)
    yline(yGridLines(k), '-', 'Color', [0.8, 0.8, 0.8]);
end

% Labels and title
xlabel('Latency (ms)', 'FontSize', 24, 'FontName', 'Helvetica');
ylabel('Relative Occurrence Frequency', 'FontSize', 24, 'FontName', 'Helvetica');
title(figure_title, 'FontSize', 30, 'FontName', 'Helvetica');

% Reference lines
yline(0, '--k');
xline(0, '--r', 'TMS', 'LineWidth', 5, ...
    'LabelHorizontalAlignment', 'center', ...
    'LabelVerticalAlignment', 'top', ...
    'LabelOrientation','horizontal', ...
    'FontSize', 28, 'FontName', 'Helvetica');

% Legend
lgd_microstate = legend(microstates, 'FontSize', 28, 'FontName', 'Helvetica');
title(lgd_microstate, 'Microstate');
lgd_microstate.Location = 'northeast';

hold off;

%% Save figure
if save_figure
    safe_figure_title = matlab.lang.makeValidName(figure_title);
    exportgraphics(gcf, [safe_figure_title, '.pdf'], 'ContentType', 'vector');
    close all;
end

end