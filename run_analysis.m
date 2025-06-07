%% Clear Workspace, Command Window, and Close All Figures
clear; clc; close all;

%% Constants and Directories
% Ensure that you are in the correct directory before running this script.
% Define the parent directory (update as needed)
parent_dir = "C:\Users\Amin\Documents\GitHub\TMS-EEG-Microstates";
cd(parent_dir);

% Time arrays and labels
common_time_array = -1000:1:1000;  % Time range from -1000 ms to 1000 ms
microstates = ["A", "B", "C", "D", "E"];  % Microstate labels
stim_sites = ["lpfc", "rpfc"];  % Stimulation sites ["lmc", "lpfc", "rmc", "rpfc"]
features2extract = {'ROF', 'RTF'}; % Features to extract (e.g., 'ROF', 'RTF')

% Baseline and stimulus periods
baseline_indices = 1:991;  % Indices for baseline period
stim_time = 1001;  % Time index for TMS onset
post_stim = 1000;  % Duration of post-stimulus period (in ms)
time_offset = 20;  % Offset to skip immediate post-stimulus artifacts
test_indices = stim_time + time_offset : stim_time + post_stim;  % Indices for testing

% Analysis parameters
nshuffles = 5000;  % Number of shuffles for permutation testing
cluster_percentile = 97.5;  % Percentile for cluster-based significance testing
apply_correction = 'bonferroni';  % Correction type ('fdr' or 'bonferroni')

% Transition intervals in milliseconds
time_window_ranges = struct(...
    'baseline', [-1000, -10], ...
    'post_early', [20, 500], ...
    'post_late', [500, 1000]);

% Flags for visualization and saving results
save_figure = true;
figure_time_range = [-200, 600];

%% Experimental Conditions
selected_dataset = "Dataset1";
selected_visit = "V1";
filter_subjects = true;  % Set to true to filter matching subjects between V1 and V2

main_analysis = false;
compare_sites = false;  % Compare different stimulation sites
compare_visits = false; % Compare V1 and V2 for a single site
single_site_for_visit_comparison = "lpfc"; % Site to use for visit comparison
compare_sham = true;

%% Data Loading and Preprocessing
% Load the selected visit data
data_file_selected = fullfile(parent_dir, 'data', sprintf('%s_%s.mat', selected_dataset, selected_visit));
% data_file_selected = fullfile(parent_dir, 'data', sprintf('%s.mat', selected_dataset));
load(data_file_selected);
data_struct_selected = data_struct;  % Store the selected visit data

% If the selected visit is V2, also load V1 data for subject matching
if strcmp(selected_visit, "V2")
    % Load V1 data
    data_file_v1 = fullfile(parent_dir, 'data', sprintf('%s_V1.mat', selected_dataset));
    load(data_file_v1);
    data_struct_v1 = data_struct;  % Store V1 data
    clear data_struct;  % Clear the loaded variable to avoid confusion
    
    % Manual approach to find matching subjects between V1 and V2
    % Get all subject IDs
    subjects_v1 = {data_struct_v1.subject};
    subjects_v2 = {data_struct_selected.subject};
    
    % Initialize empty arrays for matching indices
    idx_v1 = [];
    idx_v2 = [];
    
    % For each subject in V1, check if it exists in V2
    for i = 1:length(subjects_v1)
        for j = 1:length(subjects_v2)
            % Compare subjects (convert to string for consistent comparison)
            if isequal(subjects_v1{i}, subjects_v2{j})
                idx_v1 = [idx_v1, i];
                idx_v2 = [idx_v2, j];
                break;
            end
        end
    end
    
    fprintf('Found %d common subjects between V1 and V2.\n', length(idx_v1));
    
    % Pre-filter to only include common subjects in both datasets
    data_struct_v1 = data_struct_v1(idx_v1);
    data_struct_selected = data_struct_selected(idx_v2);
    
    % We will proceed with the selected visit data
    data_struct = data_struct_selected;
    
    % Store the V1 data for later use when filtering by stim site
    data_struct_v1_matched = data_struct_v1;
else
    % If we're working with V1 data, just keep the loaded data
    data_struct = data_struct_selected;
    clear data_struct_selected;  % Clear the temporary variable
end

% Extract relative frequencies of microstates
data_struct = extract_microstates_frequencies(data_struct, microstates, baseline_indices, time_window_ranges, features2extract);

%% Results Directory
results_dir = fullfile(parent_dir, 'results');
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

%% Test TMS-Induced Microstates
% Preallocate array for efficiency
num_subjects = length(data_struct);
num_timepoints = length(test_indices);
data_microstate = zeros(num_subjects, num_timepoints);
if compare_sites
    data_microstate_secondary = zeros(num_subjects, num_timepoints);
end

if main_analysis
% Loop through stimulation sites
for s = 1:length(stim_sites)
    current_stim_site = stim_sites(s);
    current_site_char = char(current_stim_site);

    % If V2 selected, filter subjects for this specific stim site
    if strcmp(selected_visit, "V2") && filter_subjects
        % Find subjects with data for this stim site in both V1 and V2
        valid_subjects_for_site = false(num_subjects, 1);

        for i = 1:num_subjects
            % Check if this subject has data for this site in both V1 and V2
            if isfield(data_struct(i), current_site_char) && ~isempty(data_struct(i).(current_site_char)) && ...
               isfield(data_struct_v1_matched(i), 'microstate_data') && ...
               isfield(data_struct_v1_matched(i).microstate_data, current_site_char) && ...
               ~isempty(data_struct_v1_matched(i).microstate_data.(current_site_char))
                valid_subjects_for_site(i) = true;
            end
        end

        % Create a site-specific filtered data structure
        data_struct_site_filtered = data_struct(valid_subjects_for_site);
        fprintf('For %s: Found %d subjects with data in both V1 and V2.\n', upper(current_site_char), sum(valid_subjects_for_site));
    else
        % If not V2 or not filtering, use all subjects
        data_struct_site_filtered = data_struct;
    end

    % Get the number of subjects for this analysis
    num_subjects_site = length(data_struct_site_filtered);

    % Analyze ROF (Relative Occurrence Frequencies) if selected
    if any(strcmp(features2extract, 'ROF'))
        for ss = 1:length(microstates)
            current_state = microstates(ss);

            % Prepare data for TFCE analysis
            data_microstate = zeros(num_subjects_site, num_timepoints);
            for i = 1:num_subjects_site
                if isfield(data_struct_site_filtered(i), current_site_char) && ~isempty(data_struct_site_filtered(i).(current_site_char))
                    data_microstate(i, :) = data_struct_site_filtered(i).(current_stim_site).occurrences_clr_bc.(current_state)(test_indices)';
                end
            end

            % Perform TFCE analysis
            [TFCE_Obs, TFCE_Perm, P_Values, Info, ROF_Results] = test_rof_one_condition(data_microstate, nshuffles);

            % Save results
            rof_file = fullfile(results_dir, sprintf('ROF_Results_%s_%s_%s_%s.mat', ...
                selected_dataset, selected_visit, upper(current_stim_site), current_state));
            save(rof_file, 'ROF_Results');
        end

        % Plot microstate occurrences
        if save_figure
            figure_title = sprintf('%s %s %s', selected_dataset, upper(current_stim_site), selected_visit);

            % Initialize a numerical p-values matrix filled with NaNs
            % Rows correspond to microstates, columns correspond to time points
            p_values_matrix = NaN(length(microstates), length(common_time_array));

            for ss = 1:length(microstates)
                % Extract current state
                current_state = microstates(ss);

                % Load p-values for the current state
                struct_name = fullfile(results_dir, sprintf('ROF_Results_%s_%s_%s_%s.mat', ...
                selected_dataset, selected_visit, upper(current_stim_site), current_state));
                loaded_data = load(struct_name);
                p_values = loaded_data.ROF_Results.P_Values;

                % Initialize a full p-value vector with NaNs
                full_p_values = NaN(1, length(common_time_array));

                % Assign p-values to the test period
                full_p_values(test_indices) = p_values;

                % Assign the full p-values vector to the matrix
                p_values_matrix(ss, :) = full_p_values;
            end

            plot_microstate_occurrence(data_struct_site_filtered, current_stim_site, ...
                p_values_matrix, microstates, common_time_array, ...
                figure_time_range, figure_title, save_figure);
        end
    end

    % Analyze RTF (Relative Transition Frequencies) if selected
    if any(strcmp(features2extract, 'RTF'))
        % Perform pairwise t-tests on transitions and compute Cohen's d
        RTF_Results = test_rtf(data_struct_site_filtered, current_stim_site, apply_correction);

        % Save RTF results
        rtf_file = fullfile(results_dir, sprintf('RTF_Results_%s_%s_%s.mat', ...
            selected_dataset, selected_visit, upper(current_stim_site)));
        save(rtf_file, 'RTF_Results');

        % Plot transition T-scores as heatmaps
        plot_title = sprintf('Microstates Transitions T-scores %s %s %s', ...
            selected_dataset, selected_visit, upper(current_stim_site));
        plot_transitions_t_scores_heatmap(microstates, RTF_Results, ...
            time_window_ranges.post_early, time_window_ranges.post_late, plot_title, save_figure);
    end
end
end

%%
if compare_sites
    % For comparison between sites, filter subjects that have data for both sites
    site1_char = char(stim_sites(1));
    site2_char = char(stim_sites(2));
    
    % If V2 selected, first filter subjects for both stim sites in both V1 and V2
    if strcmp(selected_visit, "V2") && filter_subjects
        % Find subjects with data for both sites in both V1 and V2
        valid_subjects_for_comparison = false(num_subjects, 1);
        
        for i = 1:num_subjects
            % Check if this subject has data for both sites in both V1 and V2
            if isfield(data_struct(i), site1_char) && ~isempty(data_struct(i).(site1_char)) && ...
               isfield(data_struct(i), site2_char) && ~isempty(data_struct(i).(site2_char)) && ...
               isfield(data_struct_v1_matched(i), 'microstate_data') && ...
               isfield(data_struct_v1_matched(i).microstate_data, site1_char) && ...
               ~isempty(data_struct_v1_matched(i).microstate_data.(site1_char)) && ...
               isfield(data_struct_v1_matched(i).microstate_data, site2_char) && ...
               ~isempty(data_struct_v1_matched(i).microstate_data.(site2_char))
                valid_subjects_for_comparison(i) = true;
            end
        end
        
        % Create filtered data structure for comparison
        data_struct_comparison = data_struct(valid_subjects_for_comparison);
        fprintf('For %s vs %s comparison: Found %d subjects with data in both V1 and V2.\n', ...
                upper(site1_char), upper(site2_char), sum(valid_subjects_for_comparison));
    else
        % Find subjects that have data for both sites (regardless of visit)
        valid_subjects_for_comparison = false(num_subjects, 1);
        
        for i = 1:num_subjects
            if isfield(data_struct(i), site1_char) && ~isempty(data_struct(i).(site1_char)) && ...
               isfield(data_struct(i), site2_char) && ~isempty(data_struct(i).(site2_char))
                valid_subjects_for_comparison(i) = true;
            end
        end
        
        data_struct_comparison = data_struct(valid_subjects_for_comparison);
        fprintf('For %s vs %s comparison: Found %d subjects with data for both sites.\n', ...
                upper(site1_char), upper(site2_char), sum(valid_subjects_for_comparison));
    end
    
    % Get the number of subjects for this comparison
    num_subjects_comparison = length(data_struct_comparison);
    
    % Analyze ROF (Relative Occurrence Frequencies) if selected
    if any(strcmp(features2extract, 'ROF'))
        for ss = 1:length(microstates)
            current_state = microstates(ss);

            % Prepare data for TFCE analysis
            data_microstate = zeros(num_subjects_comparison, num_timepoints);
            data_microstate_secondary = zeros(num_subjects_comparison, num_timepoints);
            
            for i = 1:num_subjects_comparison
                data_microstate(i, :) = data_struct_comparison(i).(stim_sites(1)).occurrences_clr_bc.(current_state)(test_indices)';
                data_microstate_secondary(i, :) = data_struct_comparison(i).(stim_sites(2)).occurrences_clr_bc.(current_state)(test_indices)';
            end

            % Perform TFCE analysis
            [TFCE_Obs, TFCE_Perm, P_Values, Info, ROF_Results] = test_rof_within_two_conditions(data_microstate, data_microstate_secondary, nshuffles);

            % Save results
            rof_file = fullfile(results_dir, sprintf('ROF_Results_%s_%s_%s_%s_%s.mat', ...
                selected_dataset, selected_visit, upper(stim_sites(1)), upper(stim_sites(2)), current_state));
            save(rof_file, 'ROF_Results');
        end

        % Plot microstate occurrences
        if save_figure
            figure_title = sprintf('%s %s %s vs. %s', selected_dataset, selected_visit, upper(stim_sites(1)), upper(stim_sites(2)));
            
            % Initialize a numerical p-values matrix filled with NaNs
            % Rows correspond to microstates, columns correspond to time points
            p_values_matrix = NaN(length(microstates), length(common_time_array));

            for ss = 1:length(microstates)
                % Extract current state
                current_state = microstates(ss);

                % Load p-values for the current state
                struct_name = fullfile(results_dir, sprintf('ROF_Results_%s_%s_%s_%s_%s.mat', ...
                    selected_dataset, selected_visit, upper(stim_sites(1)), upper(stim_sites(2)), current_state));
                loaded_data = load(struct_name);
                p_values = loaded_data.ROF_Results.P_Values;

                % Initialize a full p-value vector with NaNs
                full_p_values = NaN(1, length(common_time_array));

                % Assign p-values to the test period
                full_p_values(test_indices) = p_values;

                % Assign the full p-values vector to the matrix
                p_values_matrix(ss, :) = full_p_values;
            end

            plot_microstate_occurrence(data_struct_comparison, cellstr(stim_sites), ...
                p_values_matrix, microstates, common_time_array, ...
                figure_time_range, figure_title, save_figure);
        end
    end

    % Analyze RTF (Relative Transition Frequencies) if selected
    if any(strcmp(features2extract, 'RTF'))
        % Perform pairwise t-tests on transitions and compute Cohen's d
        RTF_Results = test_rtf_within_two_conditions(data_struct_comparison, stim_sites, apply_correction);

        % Save RTF results
        rtf_file = fullfile(results_dir, sprintf('RTF_Results_%s_%s_%s_%s.mat', ...
            selected_dataset, selected_visit, upper(stim_sites(1)), upper(stim_sites(2))));
        save(rtf_file, 'RTF_Results');

        % Plot transition T-scores as heatmaps
        plot_title = sprintf('Microstates Transitions T-scores %s %s %s vs. %s', ...
            selected_dataset, selected_visit, upper(stim_sites(1)), upper(stim_sites(2)));
        plot_transitions_t_scores_heatmap(microstates, RTF_Results, ...
            time_window_ranges.post_early, time_window_ranges.post_late, plot_title, save_figure);
    end
end

%% Compare Visit 1 and Visit 2 for specific stimulation site
if compare_visits
    % Verify that we have loaded the required data for both visits
    fprintf('Setting up visit comparison for site: %s\n', single_site_for_visit_comparison);
    
    % Check which visit is currently loaded and load the other
    if strcmp(selected_visit, "V1")
        % We have V1 loaded, need to load V2
        fprintf('Loading V2 data to compare with V1...\n');
        data_file_v2 = fullfile(parent_dir, 'data', sprintf('%s_V2.mat', selected_dataset));
        
        % Check if the V2 file exists
        if ~exist(data_file_v2, 'file')
            error('Cannot find V2 data file: %s', data_file_v2);
        end
        
        % Load V2 data
        temp_struct = load(data_file_v2);
        data_struct_v2 = temp_struct.data_struct;
        clear temp_struct;
        
        % V1 data is already loaded as data_struct
        data_struct_v1 = data_struct;
    else  % V2 is loaded
        % We already loaded V1 data during initialization
        if ~exist('data_struct_v1_matched', 'var')
            error('V1 data should have been loaded but is not available');
        end
        
        % V1 data is in data_struct_v1_matched, V2 is in data_struct
        data_struct_v1 = data_struct_v1_matched;
        data_struct_v2 = data_struct;
    end
    
    % Find matching subjects between V1 and V2
    subjects_v1 = {data_struct_v1.subject};
    subjects_v2 = {data_struct_v2.subject};
    
    % Initialize empty arrays for matching indices
    idx_v1 = [];
    idx_v2 = [];
    
    % For each subject in V1, check if it exists in V2
    for i = 1:length(subjects_v1)
        for j = 1:length(subjects_v2)
            % Compare subjects
            if isequal(subjects_v1{i}, subjects_v2{j})
                idx_v1 = [idx_v1, i];
                idx_v2 = [idx_v2, j];
                break;
            end
        end
    end
    
    fprintf('Found %d common subjects between V1 and V2.\n', length(idx_v1));
    
    % Filter to only include common subjects in both datasets
    data_struct_v1_matched = data_struct_v1(idx_v1);
    data_struct_v2_matched = data_struct_v2(idx_v2);
    
    % Extract site to compare (from variable single_site_for_visit_comparison)
    site_to_compare = single_site_for_visit_comparison;
    if ischar(site_to_compare) || isstring(site_to_compare)
        site_to_compare_char = char(site_to_compare);
    else
        error('single_site_for_visit_comparison must be a string or character array');
    end
    
    fprintf('Comparing V1 and V2 data for stimulation site: %s\n', upper(site_to_compare_char));
    
    % Extract relative frequencies of microstates for V1 and V2 data if not already done
    % Process V1 data
    if ~isfield(data_struct_v1_matched(1), site_to_compare_char) || ...
       ~isfield(data_struct_v1_matched(1).(site_to_compare_char), 'occurrences_clr_bc')
        fprintf('Extracting features for V1 data...\n');
        data_struct_v1_matched = extract_microstates_relative_frequencies(data_struct_v1_matched, microstates, baseline_indices, time_window_ranges, features2extract);
    end
    
    % Process V2 data
    if ~isfield(data_struct_v2_matched(1), site_to_compare_char) || ...
       ~isfield(data_struct_v2_matched(1).(site_to_compare_char), 'occurrences_clr_bc')
        fprintf('Extracting features for V2 data...\n');
        data_struct_v2_matched = extract_microstates_relative_frequencies(data_struct_v2_matched, microstates, baseline_indices, time_window_ranges, features2extract);
    end
    
    % Filter subjects that have data for this stimulation site in both V1 and V2
    valid_subjects_for_visit_comparison = false(length(data_struct_v1_matched), 1);
    
    for i = 1:length(data_struct_v1_matched)
        % Check if this subject has data for this site in both V1 and V2
        if isfield(data_struct_v1_matched(i), site_to_compare_char) && ...
           isfield(data_struct_v1_matched(i).(site_to_compare_char), 'occurrences_clr_bc') && ...
           isfield(data_struct_v2_matched(i), site_to_compare_char) && ...
           isfield(data_struct_v2_matched(i).(site_to_compare_char), 'occurrences_clr_bc')
            % Check if the first microstate exists in both datasets
            if isfield(data_struct_v1_matched(i).(site_to_compare_char).occurrences_clr_bc, microstates(1)) && ...
               isfield(data_struct_v2_matched(i).(site_to_compare_char).occurrences_clr_bc, microstates(1))
                valid_subjects_for_visit_comparison(i) = true;
            end
        end
    end
    
    % Create filtered data structures for visit comparison
    data_struct_v1_filtered = data_struct_v1_matched(valid_subjects_for_visit_comparison);
    data_struct_v2_filtered = data_struct_v2_matched(valid_subjects_for_visit_comparison);
    
    num_subjects_visit_comparison = sum(valid_subjects_for_visit_comparison);
    fprintf('Found %d subjects with data for site %s in both V1 and V2.\n', num_subjects_visit_comparison, upper(site_to_compare_char));
    
    if num_subjects_visit_comparison == 0
        error('No subjects found with data for site %s in both V1 and V2.', upper(site_to_compare_char));
    end
    
    % Analyze ROF (Relative Occurrence Frequencies) if selected
    if any(strcmp(features2extract, 'ROF'))
        for ss = 1:length(microstates)
            current_state = microstates(ss);
            
            % Prepare data for TFCE analysis
            data_microstate_v1 = zeros(num_subjects_visit_comparison, num_timepoints);
            data_microstate_v2 = zeros(num_subjects_visit_comparison, num_timepoints);
            
            for i = 1:num_subjects_visit_comparison
                % Extract data for V1 and V2
                data_microstate_v1(i, :) = data_struct_v1_filtered(i).(site_to_compare_char).occurrences_clr_bc.(current_state)(test_indices)';
                data_microstate_v2(i, :) = data_struct_v2_filtered(i).(site_to_compare_char).occurrences_clr_bc.(current_state)(test_indices)';
            end
            
            % Perform TFCE analysis comparing V1 vs V2
            [TFCE_Obs, TFCE_Perm, P_Values, Info, ROF_Results] = test_rof_within_two_conditions(data_microstate_v1, data_microstate_v2, nshuffles);
            
            % Save results
            rof_file = fullfile(results_dir, sprintf('ROF_Results_%s_V1_vs_V2_%s_%s.mat', ...
                selected_dataset, upper(site_to_compare_char), current_state));
            save(rof_file, 'ROF_Results');
        end
        
        % Plot microstate occurrences comparison
        if save_figure
            figure_title = sprintf('%s %s V1 vs. V2', selected_dataset, upper(site_to_compare_char));
            
            % Initialize a numerical p-values matrix filled with NaNs
            % Rows correspond to microstates, columns correspond to time points
            p_values_matrix = NaN(length(microstates), length(common_time_array));
            
            for ss = 1:length(microstates)
                % Extract current state
                current_state = microstates(ss);
                
                % Load p-values for the current state
                struct_name = fullfile(results_dir, sprintf('ROF_Results_%s_V1_vs_V2_%s_%s.mat', ...
                    selected_dataset, upper(site_to_compare_char), current_state));
                loaded_data = load(struct_name);
                p_values = loaded_data.ROF_Results.P_Values;
                
                % Initialize a full p-value vector with NaNs
                full_p_values = NaN(1, length(common_time_array));
                
                % Assign p-values to the test period
                full_p_values(test_indices) = p_values;
                
                % Assign the full p-values vector to the matrix
                p_values_matrix(ss, :) = full_p_values;
            end
            
            % Create cell array for visit comparison sites
            visit_comparison_sites = {"V1", "V2"};
            
            % Create a temporary combined data structure for visualization
            temp_data_struct = struct();
            
            % Copy data from V1 and V2 into the temporary structure
            for i = 1:num_subjects_visit_comparison
                temp_data_struct(i).subject = data_struct_v1_filtered(i).subject;
                temp_data_struct(i).V1 = data_struct_v1_filtered(i).(site_to_compare_char);
                temp_data_struct(i).V2 = data_struct_v2_filtered(i).(site_to_compare_char);
            end
            
            % Plot the comparison using the existing function
            plot_microstate_occurrence(temp_data_struct, visit_comparison_sites, ...
                p_values_matrix, microstates, common_time_array, ...
                figure_time_range, figure_title, save_figure);
        end
    end
    
    % Analyze RTF (Relative Transition Frequencies) if selected
    if any(strcmp(features2extract, 'RTF'))
        % Prepare data for RTF analysis
        % We need to create a modified data structure to use the existing functions
        rtf_data_struct = struct();
        
        % Copy relevant fields to the new structure
        for i = 1:num_subjects_visit_comparison
            rtf_data_struct(i).subject = data_struct_v1_filtered(i).subject;
            
            % Create visit fields and copy transition data
            rtf_data_struct(i).V1 = struct();
            rtf_data_struct(i).V2 = struct();
            
            % Check if transition_averages_bc exists
            if isfield(data_struct_v1_filtered(i).(site_to_compare_char), 'transition_averages_bc') && ...
               isfield(data_struct_v2_filtered(i).(site_to_compare_char), 'transition_averages_bc')
                % Copy transition averages for each visit
                rtf_data_struct(i).V1.transition_averages_bc = data_struct_v1_filtered(i).(site_to_compare_char).transition_averages_bc;
                rtf_data_struct(i).V2.transition_averages_bc = data_struct_v2_filtered(i).(site_to_compare_char).transition_averages_bc;
            end
        end
        
        % Perform RTF analysis comparing V1 vs V2
        rtf_visit_comparison_sites = ["V1", "V2"];
        RTF_Results = test_rtf_within_two_conditions(rtf_data_struct, rtf_visit_comparison_sites, apply_correction);
        
        % Save RTF results
        rtf_file = fullfile(results_dir, sprintf('RTF_Results_%s_V1_vs_V2_%s.mat', ...
            selected_dataset, upper(site_to_compare_char)));
        save(rtf_file, 'RTF_Results');
        
        % Plot transition T-scores as heatmaps
        plot_title = sprintf('Microstates Transitions T-scores %s %s V1 vs. V2', ...
            selected_dataset, upper(site_to_compare_char));
        plot_transitions_t_scores_heatmap(microstates, RTF_Results, ...
            time_window_ranges.post_early, time_window_ranges.post_late, plot_title, save_figure);
    end
end

%%
%% Compare Active vs Sham Stimulation
if compare_sham
    fprintf('Setting up active vs sham comparison...\n');
    
    % Identify active and sham sites
    active_sites = stim_sites;  % ["lpfc", "rpfc"]
    sham_sites = cellfun(@(x) [char(x) '_sham'], stim_sites, 'UniformOutput', false);
    
    % Find subjects with data for at least one active site AND at least one sham site
    valid_subjects_for_sham_comparison = false(num_subjects, 1);
    
    for i = 1:num_subjects
        % Check if subject has at least one active site
        has_active = false;
        for a = 1:length(active_sites)
            active_site_char = char(active_sites(a));
            if isfield(data_struct(i), active_site_char) && ~isempty(data_struct(i).(active_site_char))
                has_active = true;
                break;
            end
        end
        
        % Check if subject has at least one sham site
        has_sham = false;
        for s = 1:length(sham_sites)
            sham_site_char = sham_sites{s};
            if isfield(data_struct(i), sham_site_char) && ~isempty(data_struct(i).(sham_site_char))
                has_sham = true;
                break;
            end
        end
        
        % Subject is valid if they have both active and sham data
        if has_active && has_sham
            valid_subjects_for_sham_comparison(i) = true;
        end
    end
    
    % Create filtered data structure for sham comparison
    data_struct_sham_comparison = data_struct(valid_subjects_for_sham_comparison);
    num_subjects_sham = sum(valid_subjects_for_sham_comparison);
    fprintf('Found %d subjects with both active and sham data.\n', num_subjects_sham);
    
    if num_subjects_sham == 0
        error('No subjects found with both active and sham data.');
    end
    
    % Create combined active and sham data structures
    % We'll combine data from all available active sites into one condition
    % and all available sham sites into another condition
    for i = 1:num_subjects_sham
        % Initialize combined structures
        data_struct_sham_comparison(i).active = struct();
        data_struct_sham_comparison(i).sham = struct();
        
        % Find available active and sham data for this subject
        available_active_sites = {};
        available_sham_sites = {};
        
        for a = 1:length(active_sites)
            active_site_char = char(active_sites(a));
            if isfield(data_struct_sham_comparison(i), active_site_char) && ...
               ~isempty(data_struct_sham_comparison(i).(active_site_char))
                available_active_sites{end+1} = active_site_char;
            end
        end
        
        for s = 1:length(sham_sites)
            sham_site_char = sham_sites{s};
            if isfield(data_struct_sham_comparison(i), sham_site_char) && ...
               ~isempty(data_struct_sham_comparison(i).(sham_site_char))
                available_sham_sites{end+1} = sham_site_char;
            end
        end
        
        % Combine data from available sites (prioritize sites in order: lpfc, rpfc)
        % For active data
        if ~isempty(available_active_sites)
            % Take the first available active site (prioritizing lpfc over rpfc)
            primary_active_site = available_active_sites{1};
            
            % Copy data structure for occurrences and transitions
            if isfield(data_struct_sham_comparison(i).(primary_active_site), 'occurrences_clr_bc')
                data_struct_sham_comparison(i).active.occurrences_clr_bc = ...
                    data_struct_sham_comparison(i).(primary_active_site).occurrences_clr_bc;
            end
            
            if isfield(data_struct_sham_comparison(i).(primary_active_site), 'transition_averages_bc')
                data_struct_sham_comparison(i).active.transition_averages_bc = ...
                    data_struct_sham_comparison(i).(primary_active_site).transition_averages_bc;
            end
        end
        
        % For sham data
        if ~isempty(available_sham_sites)
            % Take the first available sham site
            primary_sham_site = available_sham_sites{1};
            
            % Copy data structure for occurrences and transitions
            if isfield(data_struct_sham_comparison(i).(primary_sham_site), 'occurrences_clr_bc')
                data_struct_sham_comparison(i).sham.occurrences_clr_bc = ...
                    data_struct_sham_comparison(i).(primary_sham_site).occurrences_clr_bc;
            end
            
            if isfield(data_struct_sham_comparison(i).(primary_sham_site), 'transition_averages_bc')
                data_struct_sham_comparison(i).sham.transition_averages_bc = ...
                    data_struct_sham_comparison(i).(primary_sham_site).transition_averages_bc;
            end
        end
    end
    
    % Analyze ROF (Relative Occurrence Frequencies) if selected
    if any(strcmp(features2extract, 'ROF'))
        for ss = 1:length(microstates)
            current_state = microstates(ss);
            
            % Prepare data for TFCE analysis
            data_microstate = zeros(num_subjects_sham, num_timepoints);
            data_microstate_secondary = zeros(num_subjects_sham, num_timepoints);
            
            for i = 1:num_subjects_sham
                if isfield(data_struct_sham_comparison(i).active, 'occurrences_clr_bc') && ...
                   isfield(data_struct_sham_comparison(i).active.occurrences_clr_bc, current_state) && ...
                   isfield(data_struct_sham_comparison(i).sham, 'occurrences_clr_bc') && ...
                   isfield(data_struct_sham_comparison(i).sham.occurrences_clr_bc, current_state)
                    
                    data_microstate(i, :) = data_struct_sham_comparison(i).active.occurrences_clr_bc.(current_state)(test_indices)';
                    data_microstate_secondary(i, :) = data_struct_sham_comparison(i).sham.occurrences_clr_bc.(current_state)(test_indices)';
                end
            end
            
            % Perform TFCE analysis
            [TFCE_Obs, TFCE_Perm, P_Values, Info, ROF_Results] = test_rof_within_two_conditions(data_microstate, data_microstate_secondary, nshuffles);
            
            % Save results
            rof_file = fullfile(results_dir, sprintf('ROF_Results_%s_%s_ACTIVE_vs_SHAM_%s.mat', ...
                selected_dataset, selected_visit, current_state));
            save(rof_file, 'ROF_Results');
        end
        
        % Plot microstate occurrences comparison
        if save_figure
            figure_title = sprintf('%s %s Active vs. Sham', selected_dataset, selected_visit);
            
            % Initialize a numerical p-values matrix filled with NaNs
            % Rows correspond to microstates, columns correspond to time points
            p_values_matrix = NaN(length(microstates), length(common_time_array));
            
            for ss = 1:length(microstates)
                % Extract current state
                current_state = microstates(ss);
                
                % Load p-values for the current state
                struct_name = fullfile(results_dir, sprintf('ROF_Results_%s_%s_ACTIVE_vs_SHAM_%s.mat', ...
                    selected_dataset, selected_visit, current_state));
                loaded_data = load(struct_name);
                p_values = loaded_data.ROF_Results.P_Values;
                
                % Initialize a full p-value vector with NaNs
                full_p_values = NaN(1, length(common_time_array));
                
                % Assign p-values to the test period
                full_p_values(test_indices) = p_values;
                
                % Assign the full p-values vector to the matrix
                p_values_matrix(ss, :) = full_p_values;
            end
            
            % Create cell array for comparison conditions
            comparison_conditions = {"active", "sham"};
            
            % Plot the comparison using the existing function
            plot_microstate_occurrence(data_struct_sham_comparison, comparison_conditions, ...
                p_values_matrix, microstates, common_time_array, ...
                figure_time_range, figure_title, save_figure);
        end
    end
    
    % Analyze RTF (Relative Transition Frequencies) if selected
    if any(strcmp(features2extract, 'RTF'))
        % Perform pairwise t-tests on transitions and compute Cohen's d
        rtf_comparison_conditions = ["active", "sham"];
        RTF_Results = test_rtf_within_two_conditions(data_struct_sham_comparison, rtf_comparison_conditions, apply_correction);
        
        % Save RTF results
        rtf_file = fullfile(results_dir, sprintf('RTF_Results_%s_%s_ACTIVE_vs_SHAM.mat', ...
            selected_dataset, selected_visit));
        save(rtf_file, 'RTF_Results');
        
        % Plot transition T-scores as heatmaps
        plot_title = sprintf('Microstates Transitions T-scores %s %s Active vs. Sham', ...
            selected_dataset, selected_visit);
        plot_transitions_t_scores_heatmap(microstates, RTF_Results, ...
            time_window_ranges.post_early, time_window_ranges.post_late, plot_title, save_figure);
    end
end