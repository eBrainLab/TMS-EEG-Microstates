function results = test_rtf(data_struct, stim_site, apply_correction)
% test_rtf performs statistical t-tests on the early and late post-TMS 
% transition matrices, computes Cohen's d, and applies optional p-value correction.
%
% Inputs:
%   - data_struct: Struct array containing transition matrices for early 
%                  and late post-stimulation periods.
%   - stim_site: String specifying the stimulation site field in the struct.
%   - apply_correction: String indicating the type of correction to apply 
%                       ('fdr', 'bonferroni', or empty if no correction).
%
% Output:
%   - results: Struct containing hypothesis test results, p-values, t-statistics,
%              and Cohen's d for early and late transitions.

    % Validate input data before proceeding
    valid_subjects = [];
    for kk = 1:length(data_struct)
        % Check if the subject has the required fields
        if isfield(data_struct(kk), stim_site) && ...
           isfield(data_struct(kk).(stim_site), 'transition_averages_bc') && ...
           isfield(data_struct(kk).(stim_site).transition_averages_bc, 'post_early') && ...
           isfield(data_struct(kk).(stim_site).transition_averages_bc, 'post_late')
            valid_subjects = [valid_subjects, kk];
        end
    end
    
    % Check if we found any valid subjects
    if isempty(valid_subjects)
        error(['No subjects found with transition_averages_bc data for stim_site "', stim_site, '".']);
    end
    
    % Report number of valid subjects found
    fprintf('Found %d valid subjects out of %d total for stim_site "%s".\n', ...
        length(valid_subjects), length(data_struct), stim_site);
    
    % Get dimensions from the first valid subject to preallocate arrays
    first_valid = valid_subjects(1);
    matrix_size = size(data_struct(first_valid).(stim_site).transition_averages_bc.post_early);
    num_valid_subjects = length(valid_subjects);
    
    % Preallocate arrays for transition matrices
    post_early_sample = zeros(matrix_size(1), matrix_size(2), num_valid_subjects);
    post_late_sample = zeros(matrix_size(1), matrix_size(2), num_valid_subjects);

    % Extract post-transition matrices for valid subjects only
    for idx = 1:num_valid_subjects
        kk = valid_subjects(idx);
        post_early_sample(:, :, idx) = data_struct(kk).(stim_site).transition_averages_bc.post_early;
        post_late_sample(:, :, idx) = data_struct(kk).(stim_site).transition_averages_bc.post_late;
    end

    % Perform t-test for early transitions
    [~, p_early, ~, stats] = ttest(...
        post_early_sample, 0, 'Dim', 3);
    t_early = stats.tstat;  % Extract t-statistics

    % Perform t-test for late transitions
    [~, p_late, ~, stats] = ttest(...
        post_late_sample, 0, 'Dim', 3);
    t_late = stats.tstat;  % Extract t-statistics

    % Temporarily set diagonal elements to NaN to exclude them from correction
    for ii = 1:size(p_early, 1)
        p_early(ii, ii) = NaN;
        p_late(ii, ii) = NaN;
    end

    % Apply p-value correction if specified
    if ~isempty(apply_correction)
        % Flatten p-values, excluding NaNs (diagonal elements)
        p_early_flat = p_early(~isnan(p_early));
        p_late_flat = p_late(~isnan(p_late));

        if strcmpi(apply_correction, "fdr")
            % Apply False Discovery Rate (FDR) correction
            p_early_corrected = mafdr(p_early_flat, 'BHFDR', true);
            p_late_corrected = mafdr(p_late_flat, 'BHFDR', true);
        elseif strcmpi(apply_correction, "bonferroni")
            % Apply Bonferroni correction (scaled by the number of tests)
            p_early_corrected = min(p_early_flat * numel(p_early_flat), 1);
            p_late_corrected = min(p_late_flat * numel(p_late_flat), 1);
        else
            warning('Unknown correction method "%s". No correction applied.', apply_correction);
            p_early_corrected = p_early_flat;
            p_late_corrected = p_late_flat;
        end

        % Reshape the corrected p-values back to their original size, skipping diagonal elements
        p_early(~isnan(p_early)) = p_early_corrected;
        p_late(~isnan(p_late)) = p_late_corrected;
    end

    % Update hypothesis test results based on the corrected or uncorrected p-values
    alpha = 0.05;  % Significance threshold
    h_early = p_early < alpha;
    h_late = p_late < alpha;

    % Compute Cohen's d for early and late transitions
    cohen_d_early = calculate_cohens_d(post_early_sample);
    cohen_d_late = calculate_cohens_d(post_late_sample);

    % Bundle all outputs into a struct
    results.h_early = h_early;
    results.p_early = p_early;
    results.t_early = t_early;
    results.cohen_d_early = cohen_d_early;
    results.h_late = h_late;
    results.p_late = p_late;
    results.t_late = t_late;
    results.cohen_d_late = cohen_d_late;
    results.valid_subjects = valid_subjects;
    results.num_valid_subjects = num_valid_subjects;
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