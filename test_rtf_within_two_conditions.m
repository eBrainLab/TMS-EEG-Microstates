function results = test_rtf_within_two_conditions(data_struct, stim_sites, apply_correction)
% test_rtf_within_two_conditions performs statistical analysis to evaluate
% differences in transition dynamics between two stimulation conditions, 
% with optional p-value correction.
%
% Inputs:
%   - data_struct: Primary struct array containing transition matrices.
%   - apply_correction: String indicating the type of p-value correction
%                       ('fdr', 'bonferroni', or empty if no correction).
%
% Outputs:
%   - results: Struct containing hypothesis test results, p-values, t-statistics,
%              and Cohen's d for early and late transitions.

    % Initialize Variables
    % Preallocate difference matrices for efficiency
    num_states = size(data_struct(1).(stim_sites(1)).transition_averages_bc.post_early, 1);
    num_subjects = length(data_struct);
    difference_early = zeros(num_states, num_states, num_subjects);
    difference_late = zeros(num_states, num_states, num_subjects);

    % Compute differences for both post_early and post_late transitions
    for sub = 1:num_subjects
        % Subtract secondary from primary to get difference matrices
        difference_early(:, :, sub) = data_struct(sub).(stim_sites(1)).transition_averages_bc.post_early - ...
                                       data_struct(sub).(stim_sites(2)).transition_averages_bc.post_early;
        difference_late(:, :, sub) = data_struct(sub).(stim_sites(1)).transition_averages_bc.post_late - ...
                                      data_struct(sub).(stim_sites(2)).transition_averages_bc.post_late;
    end

    % Perform One-Sample T-Tests Against Zero
    % One-sample t-test: testing if the mean of the differences is zero
    [~, p_early, ~, stats_early] = ttest(difference_early, 0, 'Dim', 3, 'Tail', 'both');
    t_early = stats_early.tstat;  % Extract t-statistics

    [~, p_late, ~, stats_late] = ttest(difference_late, 0, 'Dim', 3, 'Tail', 'both');
    t_late = stats_late.tstat;  % Extract t-statistics

    % Exclude Self-Transitions by Setting Diagonal Elements to NaN
    for ii = 1:num_states
        p_early(ii, ii) = NaN;
        p_late(ii, ii) = NaN;
    end

    % Apply Multiple Comparisons Correction if Specified
    if nargin >= 3 && ~isempty(apply_correction)
        % Flatten p-values, excluding NaNs (diagonal elements)
        p_early_flat = p_early(~isnan(p_early));
        p_late_flat = p_late(~isnan(p_late));

        if strcmpi(apply_correction, "fdr")
            % Apply False Discovery Rate (FDR) correction
            p_early_corrected = mafdr(p_early_flat, 'BHFDR', true);
            p_late_corrected = mafdr(p_late_flat, 'BHFDR', true);
        elseif strcmpi(apply_correction, "bonferroni")
            % Apply Bonferroni correction (scaled by the number of tests)
            num_tests_early = numel(p_early_flat);
            num_tests_late = numel(p_late_flat);
            p_early_corrected = min(p_early_flat * num_tests_early, 1);
            p_late_corrected = min(p_late_flat * num_tests_late, 1);
        else
            error('Unsupported correction method. Use "fdr", "bonferroni", or leave empty for no correction.');
        end

        % Reshape the corrected p-values back to their original size, skipping diagonal elements
        p_early(~isnan(p_early)) = p_early_corrected;
        p_late(~isnan(p_late)) = p_late_corrected;
    end

    % Determine Significant Transitions Based on Corrected or Uncorrected p-values
    alpha = 0.05;  % Significance threshold
    h_early = p_early < alpha;
    h_late = p_late < alpha;
    
    % Compute Cohen's d for early and late transitions
    cohen_d_early = calculate_cohens_d(difference_early);
    cohen_d_late = calculate_cohens_d(difference_late);

    % Bundle all outputs into a struct
    results.h_early = h_early;
    results.p_early = p_early;
    results.t_early = t_early;
    results.cohen_d_early = cohen_d_early;
    results.h_late = h_late;
    results.p_late = p_late;
    results.t_late = t_late;
    results.cohen_d_late = cohen_d_late;
end

function cohen_d_matrix = calculate_cohens_d(combined_data)
% Helper function to calculate Cohen's d for a given transition matrix.
%
% Inputs:
%   - combined_data: 3D matrix of transition data aggregated across subjects.
%
% Output:
%   - cohen_d_matrix: Matrix containing Cohen's d for each transition.

    transition_size = size(combined_data, 1);  % Number of states
    cohen_d_matrix = zeros(transition_size, transition_size);

    % Calculate Cohen's d for each transition
    for row = 1:transition_size
        for col = 1:transition_size
            if row ~= col  % Exclude self-transitions if needed
                combined_series = squeeze(combined_data(row, col, :));  % Nx1 vector

                % Calculate mean and standard deviation
                M = mean(combined_series);
                SD = std(combined_series);

                % Handle cases where SD is zero to avoid division by zero
                if SD == 0
                    cohen_d = 0;
                else
                    cohen_d = M / SD;
                end

                % Store Cohen's d in the matrix
                cohen_d_matrix(row, col) = cohen_d;
            else
                % Optionally set self-transitions to NaN or zero
                cohen_d_matrix(row, col) = NaN;  % Or zero if preferred
            end
        end
    end
end