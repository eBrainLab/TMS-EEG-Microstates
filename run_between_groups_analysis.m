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
stim_sites = ["lpfc", "rpfc", "lmc", "rmc", "lpc", "rpc"];  % Stimulation sites ["lmc", "lpfc", "rmc", "rpfc"]
features2extract = {'ROF', 'RTF'}; % Metrics to extract

% Baseline and stimulus periods
baseline_indices = 1:991;  % Indices for baseline period
stim_time = 1001;  % Time index for TMS onset
post_stim = 400;  % Duration of post-stimulus period (in ms)
time_offset = 50;  % Offset to skip immediate post-stimulus artifacts
test_indices = stim_time + time_offset : stim_time + post_stim;  % Indices for testing

% Analysis parameters
nshuffles = 1000;  % Number of shuffles for permutation testing
cluster_percentile = 97.5;  % Percentile for cluster-based significance testing
apply_correction = 'bonferroni';  % Correction type ('fdr' or 'bonferroni')

% Transition intervals in milliseconds
time_window_ranges = struct(...
    'baseline', [-1000, -10], ...
    'post_early', [50, 400], ...
    'post_late', [500, 1000]);

% Flags for visualization and saving results
save_figure = true;
figure_time_range = [-100, 500];

%% Dataset Selection
healthy_dataset = "Dataset_HEALTHY";
mdd_dataset = "Dataset_MDD";

%% Data Loading and Preprocessing
fprintf('Loading healthy subjects data...\n');
data_file_healthy = fullfile(parent_dir, 'data', sprintf('%s.mat', healthy_dataset));
load(data_file_healthy);
data_struct_healthy = data_struct;  % Store the healthy subjects data
clear data_struct;  % Clear the loaded variable to avoid confusion

fprintf('Loading MDD subjects data...\n');
data_file_mdd = fullfile(parent_dir, 'data', sprintf('%s.mat', mdd_dataset));
load(data_file_mdd);
data_struct_mdd = data_struct;  % Store the MDD subjects data
clear data_struct;  % Clear the loaded variable to avoid confusion

% Extract frequencies of microstates for both groups
fprintf('Extracting features for healthy subjects...\n');
data_struct_healthy = extract_microstates_frequencies(data_struct_healthy, microstates, baseline_indices, time_window_ranges, features2extract);

fprintf('Extracting features for MDD subjects...\n');
data_struct_mdd = extract_microstates_frequencies(data_struct_mdd, microstates, baseline_indices, time_window_ranges, features2extract);

%% Results Directory
results_dir = fullfile(parent_dir, 'results');
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

%% Compare Between Groups for Each Stimulation Site
for s = 1:length(stim_sites)
    current_stim_site = stim_sites(s);
    current_site_char = char(current_stim_site);
    
    fprintf('\n====== Processing stimulation site: %s ======\n', upper(current_site_char));
    
    % Filter subjects with valid data for this stimulation site
    valid_healthy_subjects = [];
    for i = 1:length(data_struct_healthy)
        if isfield(data_struct_healthy(i), current_site_char)
            valid_healthy_subjects = [valid_healthy_subjects, i];
        end
    end
    
    valid_mdd_subjects = [];
    for i = 1:length(data_struct_mdd)
        if isfield(data_struct_mdd(i), current_site_char)
            valid_mdd_subjects = [valid_mdd_subjects, i];
        end
    end
    
    % Create filtered datasets with only valid subjects
    data_struct_healthy_filtered = data_struct_healthy(valid_healthy_subjects);
    data_struct_mdd_filtered = data_struct_mdd(valid_mdd_subjects);
    
    num_healthy_subjects = length(data_struct_healthy_filtered);
    num_mdd_subjects = length(data_struct_mdd_filtered);
    
    fprintf('Found %d valid healthy subjects and %d valid MDD subjects for site %s.\n', ...
            num_healthy_subjects, num_mdd_subjects, upper(current_site_char));
    
    if num_healthy_subjects == 0 || num_mdd_subjects == 0
        warning('Missing data for one or both groups for site %s. Skipping this site.', upper(current_site_char));
        continue;
    end
    
    % Process Occurrence Frequencies (AOF and ROF)
    % Process both absolute and relative occurrence frequencies
    occurrence_measures = {};
    if any(strcmp(features2extract, 'AOF'))
        occurrence_measures{end+1} = 'AOF';
    end
    if any(strcmp(features2extract, 'ROF'))
        occurrence_measures{end+1} = 'ROF';
    end
    
    for m = 1:length(occurrence_measures)
        current_measure = occurrence_measures{m};
        field_name = '';
        result_prefix = '';
        
        % Set appropriate field names based on the measure
        if strcmp(current_measure, 'AOF')
            field_name = 'occurrences_clr';
            result_prefix = 'AOF';
            fprintf('Starting AOF (Absolute Occurrence Frequencies) analysis for %s...\n', upper(current_site_char));
        elseif strcmp(current_measure, 'ROF')
            field_name = 'occurrences_clr_bc';
            result_prefix = 'ROF';
            fprintf('Starting ROF (Relative Occurrence Frequencies) analysis for %s...\n', upper(current_site_char));
        end
        
        % Check if the necessary fields exist
        valid_healthy_for_measure = [];
        for i = 1:num_healthy_subjects
            if isfield(data_struct_healthy_filtered(i).(current_site_char), field_name) && ...
               isfield(data_struct_healthy_filtered(i).(current_site_char).(field_name), microstates(1))
                valid_healthy_for_measure = [valid_healthy_for_measure, i];
            end
        end
        
        valid_mdd_for_measure = [];
        for i = 1:num_mdd_subjects
            if isfield(data_struct_mdd_filtered(i).(current_site_char), field_name) && ...
               isfield(data_struct_mdd_filtered(i).(current_site_char).(field_name), microstates(1))
                valid_mdd_for_measure = [valid_mdd_for_measure, i];
            end
        end
        
        fprintf('Found %d healthy subjects and %d MDD subjects with valid %s data for site %s.\n', ...
                length(valid_healthy_for_measure), length(valid_mdd_for_measure), current_measure, upper(current_site_char));
        
        if isempty(valid_healthy_for_measure) || isempty(valid_mdd_for_measure)
            warning('No valid %s data for one or both groups for site %s. Skipping %s analysis.', ...
                    current_measure, upper(current_site_char), current_measure);
            continue;
        end
        
        % Filter for subjects with valid data for this measure
        healthy_filtered_for_measure = data_struct_healthy_filtered(valid_healthy_for_measure);
        mdd_filtered_for_measure = data_struct_mdd_filtered(valid_mdd_for_measure);
        
        num_healthy_for_measure = length(healthy_filtered_for_measure);
        num_mdd_for_measure = length(mdd_filtered_for_measure);
        
        for ss = 1:length(microstates)
            current_state = microstates(ss);
            
            % Prepare data for TFCE analysis
            data_healthy = zeros(num_healthy_for_measure, length(test_indices));
            data_mdd = zeros(num_mdd_for_measure, length(test_indices));
            
            for i = 1:num_healthy_for_measure
                data_healthy(i, :) = healthy_filtered_for_measure(i).(current_site_char).(field_name).(current_state)(test_indices)';
            end
            
            for i = 1:num_mdd_for_measure
                data_mdd(i, :) = mdd_filtered_for_measure(i).(current_site_char).(field_name).(current_state)(test_indices)';
            end
            
            % Call appropriate testing function for between-groups comparison
            [TFCE_Obs, TFCE_Perm, P_Values, Info, Results] = test_rof_between_groups(data_mdd, data_healthy, nshuffles);
            
            % Save results with appropriate prefix
            results_file = fullfile(results_dir, sprintf('%s_Results_MDD_vs_HEALTHY_%s_%s.mat', ...
                                   result_prefix, upper(current_site_char), current_state));
            if strcmp(current_measure, 'AOF')
                AOF_Results = Results;
                save(results_file, 'AOF_Results');
            else  % ROF
                ROF_Results = Results;
                save(results_file, 'ROF_Results');
            end
            fprintf('Saved %s results for microstate %s\n', current_measure, current_state);
        end
        
        % Plot microstate occurrences comparison
        if save_figure
            figure_title = sprintf('MDD vs HEALTHY %s - %s', current_measure, upper(current_site_char));
            
            % Initialize a numerical p-values matrix filled with NaNs
            p_values_matrix = NaN(length(microstates), length(common_time_array));
            
            % Create temporary data_struct to be used by plot_microstate_occurrence_between_groups
            temp_data_struct_healthy = healthy_filtered_for_measure;
            temp_data_struct_mdd = mdd_filtered_for_measure;
            
            % Load the p-values for each microstate
            for ss = 1:length(microstates)
                current_state = microstates(ss);
                
                % Load p-values for the current state
                struct_name = fullfile(results_dir, sprintf('%s_Results_MDD_vs_HEALTHY_%s_%s.mat', ...
                           result_prefix, upper(current_site_char), current_state));
                loaded_data = load(struct_name);
                
                % Get the field name from the loaded struct
                if strcmp(current_measure, 'AOF')
                    p_values = loaded_data.AOF_Results.P_Values;
                else  % ROF
                    p_values = loaded_data.ROF_Results.P_Values;
                end
                
                % Initialize a full p-value vector with NaNs
                full_p_values = NaN(1, length(common_time_array));
                
                % Assign p-values to the test period
                full_p_values(test_indices) = p_values;
                
                % Assign the full p-values vector to the matrix
                p_values_matrix(ss, :) = full_p_values;
            end
            
            % Call the specialized between-groups plotting function
            % We need to specify which measure to plot (AOF or ROF)
            plot_microstate_occurrence_between_groups(...
                temp_data_struct_mdd, ...
                temp_data_struct_healthy, ...
                current_site_char, ...
                p_values_matrix, ...
                microstates, ...
                common_time_array, ...
                figure_time_range, ...
                figure_title, ...
                save_figure, ...
                field_name);  % Pass the field name to indicate which data to plot
        end
    end
    
    % Process Transition Frequencies (ATF and RTF)
    % Process both absolute and relative transition frequencies
    transition_measures = {};
    if any(strcmp(features2extract, 'ATF'))
        transition_measures{end+1} = 'ATF';
    end
    if any(strcmp(features2extract, 'RTF'))
        transition_measures{end+1} = 'RTF';
    end
    
    for m = 1:length(transition_measures)
        current_measure = transition_measures{m};
        field_name = '';
        result_prefix = '';
        
        % Set appropriate field names based on the measure
        if strcmp(current_measure, 'ATF')
            field_name = 'transition_averages';
            result_prefix = 'ATF';
            fprintf('Starting ATF (Absolute Transition Frequencies) analysis for %s...\n', upper(current_site_char));
        elseif strcmp(current_measure, 'RTF')
            field_name = 'transition_averages_bc';
            result_prefix = 'RTF';
            fprintf('Starting RTF (Relative Transition Frequencies) analysis for %s...\n', upper(current_site_char));
        end
        
        % Check if the necessary fields exist
        valid_healthy_for_measure = [];
        for i = 1:num_healthy_subjects
            if isfield(data_struct_healthy_filtered(i).(current_site_char), field_name) && ...
               isfield(data_struct_healthy_filtered(i).(current_site_char).(field_name), 'post_early') && ...
               isfield(data_struct_healthy_filtered(i).(current_site_char).(field_name), 'post_late')
                valid_healthy_for_measure = [valid_healthy_for_measure, i];
            end
        end
        
        valid_mdd_for_measure = [];
        for i = 1:num_mdd_subjects
            if isfield(data_struct_mdd_filtered(i).(current_site_char), field_name) && ...
               isfield(data_struct_mdd_filtered(i).(current_site_char).(field_name), 'post_early') && ...
               isfield(data_struct_mdd_filtered(i).(current_site_char).(field_name), 'post_late')
                valid_mdd_for_measure = [valid_mdd_for_measure, i];
            end
        end
        
        fprintf('Found %d healthy subjects and %d MDD subjects with valid %s data for site %s.\n', ...
                length(valid_healthy_for_measure), length(valid_mdd_for_measure), current_measure, upper(current_site_char));
        
        if isempty(valid_healthy_for_measure) || isempty(valid_mdd_for_measure)
            warning('No valid %s data for one or both groups for site %s. Skipping %s analysis.', ...
                    current_measure, upper(current_site_char), current_measure);
            continue;
        end
        
        % Filter for subjects with valid data for this measure
        healthy_filtered_for_measure = data_struct_healthy_filtered(valid_healthy_for_measure);
        mdd_filtered_for_measure = data_struct_mdd_filtered(valid_mdd_for_measure);
        
        % Call appropriate testing function for between-groups comparison
        % We need to pass the field_name to the test function to indicate which data to use
        Results = test_rtf_between_groups(mdd_filtered_for_measure, healthy_filtered_for_measure, ...
                                         current_site_char, apply_correction, field_name);
        
        % Save results with appropriate prefix
        results_file = fullfile(results_dir, sprintf('%s_Results_MDD_vs_HEALTHY_%s.mat', ...
                               result_prefix, upper(current_site_char)));
        
        if strcmp(current_measure, 'ATF')
            ATF_Results = Results;
            save(results_file, 'ATF_Results');
        else  % RTF
            RTF_Results = Results;
            save(results_file, 'RTF_Results');
        end
        
        % Plot transition T-scores as heatmaps
        plot_title = sprintf('Microstates %s T-scores MDD vs HEALTHY (%s)', ...
                             current_measure, upper(current_site_char));
        
        plot_transitions_t_scores_heatmap(microstates, Results, ...
                                         time_window_ranges.post_early, time_window_ranges.post_late, ...
                                         plot_title, save_figure);
    end
end

fprintf('Between-groups analysis for all metrics complete.\n');