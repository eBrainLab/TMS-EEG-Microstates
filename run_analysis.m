% -------------------------------------------------------------------------
% RUN_ANALYSIS - Main analysis script for TMS-EEG microstate dynamics
%
% This script orchestrates the complete analysis pipeline for investigating
% how TMS affects EEG microstate occurrence and transition patterns.
%
% Analysis steps:
%   1. Load preprocessed microstate data
%   2. Extract occurrence (ROF) and transition (RTF) frequencies
%   3. Perform statistical testing with TFCE
%   4. Generate publication-ready figures
%
% Outputs saved to results/ directory:
%   - ROF_Results_*.mat - Occurrence frequency statistics
%   - RTF_Results_*.mat - Transition frequency statistics
%   - Figures in PDF/PNG format
%
% See also: EXTRACT_MICROSTATES_FREQUENCIES, TEST_ROF_ONE_CONDITION, TEST_RTF
% -------------------------------------------------------------------------

%% Clear workspace and initialize
clear; clc; close all;


%% Define analysis parameters
% Time vectors and microstate labels
common_time_array = -1000:1:1000;               % Time range (ms)
microstates = ["A", "B", "C", "D", "E"];        % Microstate labels
stim_sites = ["lpfc", "rpfc", "lmc"];           % Stimulation sites
features2extract = {'ROF', 'RTF'};              % Features to extract

% Define analysis periods
baseline_indices = 1:991;                       % Baseline period indices
stim_time = 1001;                               % TMS onset index
post_stim = 1000;                               % Post-stim duration (ms)
time_offset = 20;                               % Skip immediate artifacts (ms)
test_indices = stim_time + time_offset : stim_time + post_stim;

% Statistical parameters
nshuffles = 5000;                               % Permutation count
cluster_percentile = 97.5;                      % Cluster threshold
apply_correction = 'bonferroni';                % Multiple comparison correction

% Time windows for analysis
time_window_ranges = struct(...
    'baseline', [-1000, -10], ...
    'post_tms', [20, 500]);

% Visualization parameters
save_figure = true;
figure_time_range = [-200, 600];

%% Select dataset and analysis type
% Dataset selection
selected_dataset = "Dataset1";
selected_visit = "Demo";

% Analysis flags
main_analysis = true;                           % Single-site analysis
compare_sites = true;                           % Between-site comparison

%% Load data
fprintf('Loading data...\n');
data_file_selected = fullfile('data', ...
    sprintf('%s_%s.mat', selected_dataset, selected_visit));
load(data_file_selected);

% Extract microstate frequencies
fprintf('Extracting microstate frequencies...\n');
data_struct = extract_microstates_frequencies(data_struct, microstates, ...
    baseline_indices, time_window_ranges, features2extract);

%% Create results directory
results_dir = fullfile('results');
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

%% Main analysis - Single stimulation sites
if main_analysis
    fprintf('\n=== Starting main analysis ===\n');
    
    % Process each stimulation site
    for s = 1:length(stim_sites)
        current_stim_site = stim_sites(s);
        current_site_char = char(current_stim_site);
        
        fprintf('\nAnalyzing site: %s\n', upper(current_site_char));
        
        % Get number of subjects
        num_subjects_site = length(data_struct);
        
        %% Analyze ROF (Relative Occurrence Frequencies)
        if any(strcmp(features2extract, 'ROF'))
            fprintf('Computing ROF statistics...\n');
            
            for ss = 1:length(microstates)
                current_state = microstates(ss);
                
                % Prepare data matrix
                data_microstate = zeros(num_subjects_site, length(test_indices));
                for i = 1:num_subjects_site
                    if isfield(data_struct(i), current_site_char) && ...
                       ~isempty(data_struct(i).(current_site_char))
                        data_microstate(i, :) = data_struct(i).(current_stim_site)...
                            .occurrences_clr_bc.(current_state)(test_indices)';
                    end
                end
                
                % Perform TFCE analysis
                [TFCE_Obs, TFCE_Perm, P_Values, Info, ROF_Results] = ...
                    test_rof_one_condition(data_microstate, nshuffles);
                
                % Save results
                rof_file = fullfile(results_dir, sprintf('ROF_Results_%s_%s_%s_%s.mat', ...
                    selected_dataset, selected_visit, upper(current_stim_site), current_state));
                save(rof_file, 'ROF_Results');
            end
            
            % Generate occurrence plot
            if save_figure
                figure_title = sprintf('%s %s %s', selected_dataset, ...
                    upper(current_stim_site), selected_visit);
                
                % Compile p-values matrix
                p_values_matrix = NaN(length(microstates), length(common_time_array));
                
                for ss = 1:length(microstates)
                    current_state = microstates(ss);
                    
                    % Load saved results
                    struct_name = fullfile(results_dir, ...
                        sprintf('ROF_Results_%s_%s_%s_%s.mat', ...
                        selected_dataset, selected_visit, ...
                        upper(current_stim_site), current_state));
                    loaded_data = load(struct_name);
                    p_values = loaded_data.ROF_Results.P_Values;
                    
                    % Create full p-value vector
                    full_p_values = NaN(1, length(common_time_array));
                    full_p_values(test_indices) = p_values;
                    p_values_matrix(ss, :) = full_p_values;
                end
                
                % Create plot
                plot_microstate_occurrence(data_struct, current_stim_site, ...
                    p_values_matrix, microstates, common_time_array, ...
                    figure_time_range, figure_title, save_figure);
            end
        end
        
        %% Analyze RTF (Relative Transition Frequencies)
        if any(strcmp(features2extract, 'RTF'))
            fprintf('Computing RTF statistics...\n');
            
            % Perform transition analysis
            RTF_Results = test_rtf(data_struct, current_stim_site, apply_correction);
            
            % Save results
            rtf_file = fullfile(results_dir, sprintf('RTF_Results_%s_%s_%s.mat', ...
                selected_dataset, selected_visit, upper(current_stim_site)));
            save(rtf_file, 'RTF_Results');
            
            % Generate heatmap
            plot_title = sprintf('Microstates Transitions T-scores %s %s %s', ...
                selected_dataset, selected_visit, upper(current_stim_site));
            plot_transitions_t_scores_heatmap(microstates, RTF_Results, ...
                time_window_ranges.post_tms, plot_title, save_figure);
        end
    end
end

%% Between-site comparison analysis
if compare_sites
    fprintf('\n=== Starting between-site comparison ===\n');
    
    % Define all site pairs to compare
    site_pairs = [1 2; 1 3; 2 3];  % lpfc vs rpfc, lpfc vs lmc, rpfc vs lmc
    
    % Process each site pair
    for pair_idx = 1:size(site_pairs, 1)
        site1_idx = site_pairs(pair_idx, 1);
        site2_idx = site_pairs(pair_idx, 2);
        site1_char = char(stim_sites(site1_idx));
        site2_char = char(stim_sites(site2_idx));
        
        fprintf('\n--- Comparing %s vs %s ---\n', upper(site1_char), upper(site2_char));
        
        % Find subjects with data for both sites
        num_subjects_total = length(data_struct);
        valid_subjects_for_comparison = false(num_subjects_total, 1);
        
        for i = 1:num_subjects_total
            if isfield(data_struct(i), site1_char) && ...
               ~isempty(data_struct(i).(site1_char)) && ...
               isfield(data_struct(i), site2_char) && ...
               ~isempty(data_struct(i).(site2_char))
                valid_subjects_for_comparison(i) = true;
            end
        end
        
        % Filter data
        data_struct_comparison = data_struct(valid_subjects_for_comparison);
        num_subjects_comparison = length(data_struct_comparison);
        
        fprintf('Found %d subjects with data for both %s and %s.\n', ...
                num_subjects_comparison, upper(site1_char), upper(site2_char));
        
        %% Analyze ROF comparison
        if any(strcmp(features2extract, 'ROF'))
            fprintf('Computing ROF comparison statistics...\n');
            
            for ss = 1:length(microstates)
                current_state = microstates(ss);
                
                % Prepare data matrices
                data_microstate = zeros(num_subjects_comparison, length(test_indices));
                data_microstate_secondary = zeros(num_subjects_comparison, length(test_indices));
                
                for i = 1:num_subjects_comparison
                    data_microstate(i, :) = data_struct_comparison(i).(site1_char)...
                        .occurrences_clr_bc.(current_state)(test_indices)';
                    data_microstate_secondary(i, :) = data_struct_comparison(i).(site2_char)...
                        .occurrences_clr_bc.(current_state)(test_indices)';
                end
                
                % Perform paired TFCE analysis
                [TFCE_Obs, TFCE_Perm, P_Values, Info, ROF_Results] = ...
                    test_rof_within_two_conditions(data_microstate, ...
                    data_microstate_secondary, nshuffles);
                
                % Save results
                rof_file = fullfile(results_dir, ...
                    sprintf('ROF_Results_%s_%s_%s_%s_%s.mat', ...
                    selected_dataset, selected_visit, upper(site1_char), ...
                    upper(site2_char), current_state));
                save(rof_file, 'ROF_Results');
            end
            
            % Generate comparison plot
            if save_figure
                figure_title = sprintf('%s %s %s vs. %s', selected_dataset, ...
                    selected_visit, upper(site1_char), upper(site2_char));
                
                % Compile p-values matrix
                p_values_matrix = NaN(length(microstates), length(common_time_array));
                
                for ss = 1:length(microstates)
                    current_state = microstates(ss);
                    
                    % Load saved results
                    struct_name = fullfile(results_dir, ...
                        sprintf('ROF_Results_%s_%s_%s_%s_%s.mat', ...
                        selected_dataset, selected_visit, upper(site1_char), ...
                        upper(site2_char), current_state));
                    loaded_data = load(struct_name);
                    p_values = loaded_data.ROF_Results.P_Values;
                    
                    % Create full p-value vector
                    full_p_values = NaN(1, length(common_time_array));
                    full_p_values(test_indices) = p_values;
                    p_values_matrix(ss, :) = full_p_values;
                end
                
                % Create plot
                plot_microstate_occurrence(data_struct_comparison, {site1_char, site2_char}, ...
                    p_values_matrix, microstates, common_time_array, ...
                    figure_time_range, figure_title, save_figure);
            end
        end
        
        %% Analyze RTF comparison
        if any(strcmp(features2extract, 'RTF'))
            fprintf('Computing RTF comparison statistics...\n');
            
            % Perform paired transition analysis
            RTF_Results = test_rtf_within_two_conditions(data_struct_comparison, ...
                [stim_sites(site1_idx), stim_sites(site2_idx)], apply_correction);
            
            % Save results
            rtf_file = fullfile(results_dir, sprintf('RTF_Results_%s_%s_%s_%s.mat', ...
                selected_dataset, selected_visit, upper(site1_char), ...
                upper(site2_char)));
            save(rtf_file, 'RTF_Results');
            
            % Generate heatmap
            plot_title = sprintf('Microstates Transitions T-scores %s %s %s vs. %s', ...
                selected_dataset, selected_visit, upper(site1_char), ...
                upper(site2_char));
            plot_transitions_t_scores_heatmap(microstates, RTF_Results, ...
                time_window_ranges.post_tms, plot_title, save_figure);
        end
    end  % End of site pairs loop
end

%% Complete
fprintf('\n=== Analysis complete ===\n');
fprintf('Results saved to: %s\n', results_dir);