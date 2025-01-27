%% Clear Workspace, Command Window, and Close All Figures
clear; clc; close all;

%% Constants and Directories
% Ensure that you are in the correct directory before running this script.
% Define the parent directory (update as needed)
parent_dir = "C:\Users\Amin\Documents\Pre_Post_Microstates_LPFC_LMC\Analysis_Scripts";
cd(parent_dir);

% Time arrays and labels
common_time_array = -1000:1:1000;  % Time range from -1000 ms to 1000 ms
microstates = ["A", "B", "C", "D", "E"];  % Microstate labels
stim_sites = ["lpfc", "lmc"];  % Stimulation sites
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
compare_sites = true;
selected_dataset = "Dataset1";
selected_visit = "Visit1";

%% Data Loading and Preprocessing
% Load data
data_file = fullfile(parent_dir, 'data', sprintf('%s_%s.mat', selected_dataset, selected_visit));
load(data_file);

% Extract relative frequencies of microstates
data_struct = extract_microstates_relative_frequencies(data_struct, microstates, baseline_indices, time_window_ranges, features2extract);

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

% Loop through stimulation sites
for s = 1:length(stim_sites)
    current_stim_site = stim_sites(s);

    % Analyze ROF (Relative Occurrence Frequencies) if selected
    if any(strcmp(features2extract, 'ROF'))
        for ss = 1:length(microstates)
            current_state = microstates(ss);

            % Prepare data for TFCE analysis
            for i = 1:num_subjects
                data_microstate(i, :) = data_struct(i).(current_stim_site).occurrences_clr_bc.(current_state)(test_indices)';
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

            plot_microstate_occurrence(data_struct, current_stim_site, ...
                p_values_matrix, microstates, common_time_array, ...
                figure_time_range, figure_title, save_figure);
        end
    end

    % Analyze RTF (Relative Transition Frequencies) if selected
    if any(strcmp(features2extract, 'RTF'))
        % Perform pairwise t-tests on transitions and compute Cohen's d
        RTF_Results = test_rtf(data_struct, current_stim_site, apply_correction);

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

if compare_sites

    % Loop through stimulation sites
    % Analyze ROF (Relative Occurrence Frequencies) if selected
    if any(strcmp(features2extract, 'ROF'))
        for ss = 1:length(microstates)
            current_state = microstates(ss);

            % Prepare data for TFCE analysis
            for i = 1:num_subjects
                data_microstate(i, :) = data_struct(i).(stim_sites(1)).occurrences_clr_bc.(current_state)(test_indices)';
                data_microstate_secondary(i, :) = data_struct(i).(stim_sites(2)).occurrences_clr_bc.(current_state)(test_indices)';
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

            plot_microstate_occurrence(data_struct, stim_sites, ...
                p_values_matrix, microstates, common_time_array, ...
                figure_time_range, figure_title, save_figure);
        end

    end

    % Analyze RTF (Relative Transition Frequencies) if selected
    if any(strcmp(features2extract, 'RTF'))
        % Perform pairwise t-tests on transitions and compute Cohen's d
        RTF_Results = test_rtf_within_two_conditions(data_struct, stim_sites, apply_correction);

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