function results = test_rtf_between_groups(data_struct_healthy, data_struct_mdd, stim_site, apply_correction, field_name)
% test_rtf_between_groups performs statistical analysis to evaluate
% differences in transition dynamics between healthy and MDD groups, 
% with optional p-value correction.
%
% Inputs:
%   - data_struct_healthy: Struct array for healthy subjects.
%   - data_struct_mdd: Struct array for MDD subjects.
%   - stim_site: String specifying the stimulation site field in both structs.
%   - apply_correction: String for p-value correction ('fdr', 'bonferroni', or empty).
%   - field_name: String specifying which field to analyze ('transition_averages' for ATF,
%                'transition_averages_bc' for RTF). Default is 'transition_averages_bc'.
%
% Outputs:
%   - results: Struct containing hypothesis test results, p-values, t-statistics,
%              and Cohen's d for early and late transitions.

    % Set default field name if not provided
    if nargin < 5 || isempty(field_name)
        field_name = 'transition_averages_bc'; % Default to RTF for backward compatibility
    end

    % Determine metric type from field name for logging purposes
    if strcmp(field_name, 'transition_averages')
        metric_type = 'ATF';
    else
        metric_type = 'RTF';
    end
    
    % Get dimensions
    num_healthy = length(data_struct_healthy);
    num_mdd = length(data_struct_mdd);
    
    fprintf('Analyzing %s transitions for %d healthy subjects vs %d MDD subjects\n', ...
            metric_type, num_healthy, num_mdd);
    
    % Get matrix size from first healthy subject with valid data
    matrix_size = [];
    for i = 1:num_healthy
        if isfield(data_struct_healthy(i).(stim_site), field_name) && ...
           isfield(data_struct_healthy(i).(stim_site).(field_name), 'post_early')
            matrix_size = size(data_struct_healthy(i).(stim_site).(field_name).post_early);
            break;
        end
    end
    
    if isempty(matrix_size)
        error('Cannot find %s field in any healthy subjects data', field_name);
    end
    
    % Initialize arrays to store transition data
    % We'll collect early and late transition matrices for each group
    healthy_early = zeros(matrix_size(1), matrix_size(2), num_healthy);
    healthy_late = zeros(matrix_size(1), matrix_size(2), num_healthy);
    mdd_early = zeros(matrix_size(1), matrix_size(2), num_mdd);
    mdd_late = zeros(matrix_size(1), matrix_size(2), num_mdd);
    
    % Extract transition matrices for healthy subjects
    valid_healthy_count = 0;
    for i = 1:num_healthy
        if isfield(data_struct_healthy(i).(stim_site), field_name) && ...
           isfield(data_struct_healthy(i).(stim_site).(field_name), 'post_early') && ...
           isfield(data_struct_healthy(i).(stim_site).(field_name), 'post_late')
            valid_healthy_count = valid_healthy_count + 1;
            healthy_early(:, :, valid_healthy_count) = data_struct_healthy(i).(stim_site).(field_name).post_early;
            healthy_late(:, :, valid_healthy_count) = data_struct_healthy(i).(stim_site).(field_name).post_late;
        end
    end
    
    % Trim arrays if needed
    if valid_healthy_count < num_healthy
        fprintf('Note: Only %d out of %d healthy subjects had valid %s data.\n', ...
                valid_healthy_count, num_healthy, metric_type);
        healthy_early = healthy_early(:, :, 1:valid_healthy_count);
        healthy_late = healthy_late(:, :, 1:valid_healthy_count);
    end
    
    % Extract transition matrices for MDD subjects
    valid_mdd_count = 0;
    for i = 1:num_mdd
        if isfield(data_struct_mdd(i).(stim_site), field_name) && ...
           isfield(data_struct_mdd(i).(stim_site).(field_name), 'post_early') && ...
           isfield(data_struct_mdd(i).(stim_site).(field_name), 'post_late')
            valid_mdd_count = valid_mdd_count + 1;
            mdd_early(:, :, valid_mdd_count) = data_struct_mdd(i).(stim_site).(field_name).post_early;
            mdd_late(:, :, valid_mdd_count) = data_struct_mdd(i).(stim_site).(field_name).post_late;
        end
    end
    
    % Trim arrays if needed
    if valid_mdd_count < num_mdd
        fprintf('Note: Only %d out of %d MDD subjects had valid %s data.\n', ...
                valid_mdd_count, num_mdd, metric_type);
        mdd_early = mdd_early(:, :, 1:valid_mdd_count);
        mdd_late = mdd_late(:, :, 1:valid_mdd_count);
    end
    
    % Check if we have sufficient data
    if valid_healthy_count == 0 || valid_mdd_count == 0
        error('Insufficient data: %d healthy and %d MDD subjects with valid %s data.', ...
              valid_healthy_count, valid_mdd_count, metric_type);
    end
    
    % Perform independent samples t-test for early transitions
    [~, p_early, ~, stats_early] = ttest2(healthy_early, mdd_early, 'Dim', 3);
    t_early = stats_early.tstat;  % Extract t-statistics
    
    % Perform independent samples t-test for late transitions
    [~, p_late, ~, stats_late] = ttest2(healthy_late, mdd_late, 'Dim', 3);
    t_late = stats_late.tstat;  % Extract t-statistics
    
    % Set diagonal elements to NaN (self-transitions)
    for ii = 1:matrix_size(1)
        p_early(ii, ii) = NaN;
        p_late(ii, ii) = NaN;
        t_early(ii, ii) = NaN;
        t_late(ii, ii) = NaN;
    end
    
    % Apply multiple comparisons correction if specified
    if ~isempty(apply_correction)
        % Flatten p-values, excluding NaNs (diagonal elements)
        p_early_flat = p_early(~isnan(p_early));
        p_late_flat = p_late(~isnan(p_late));
        
        if strcmpi(apply_correction, "fdr")
            % Apply False Discovery Rate (FDR) correction
            p_early_corrected = mafdr(p_early_flat, 'BHFDR', true);
            p_late_corrected = mafdr(p_late_flat, 'BHFDR', true);
            fprintf('Applied FDR correction to p-values.\n');
        elseif strcmpi(apply_correction, "bonferroni")
            % Apply Bonferroni correction
            num_tests = numel(p_early_flat);
            p_early_corrected = min(p_early_flat * num_tests, 1);
            p_late_corrected = min(p_late_flat * num_tests, 1);
            fprintf('Applied Bonferroni correction to p-values (adjusted for %d tests).\n', num_tests);
        else
            warning('Unknown correction method "%s". No correction applied.', apply_correction);
            p_early_corrected = p_early_flat;
            p_late_corrected = p_late_flat;
        end

        % Reshape the corrected p-values back to matrix form
        p_early(~isnan(p_early)) = p_early_corrected;
        p_late(~isnan(p_late)) = p_late_corrected;
    end

    % Determine statistical significance based on corrected p-values
    alpha = 0.05;  % Significance level
    h_early = p_early < alpha;
    h_late = p_late < alpha;

    % Calculate effect sizes (Cohen's d) for each transition
    cohen_d_early = calculate_between_groups_cohens_d(healthy_early, mdd_early);
    cohen_d_late = calculate_between_groups_cohens_d(healthy_late, mdd_late);

    % Count significant transitions
    num_sig_early = sum(h_early(:));
    num_sig_late = sum(h_late(:));
    
    % Prepare correction text for display
    if ~isempty(apply_correction)
        correction_text = [', ' upper(apply_correction) ' corrected'];
    else
        correction_text = '';
    end
    
    fprintf('Found %d significant early transitions and %d significant late transitions (α=%.2f%s).\n', ...
            num_sig_early, num_sig_late, alpha, correction_text);

    % Bundle results into a struct
    results.h_early = h_early;
    results.p_early = p_early;
    results.t_early = t_early;
    results.cohen_d_early = cohen_d_early;
    results.h_late = h_late;
    results.p_late = p_late;
    results.t_late = t_late;
    results.cohen_d_late = cohen_d_late;
    results.num_healthy = valid_healthy_count;
    results.num_mdd = valid_mdd_count;
    results.metric_type = metric_type;
    results.field_name = field_name;
end

function cohen_d_matrix = calculate_between_groups_cohens_d(group1_data, group2_data)
% Calculate Cohen's d for between-groups effect size
%
% Inputs:
%   - group1_data: 3D matrix [states x states x subjects] for group 1 (healthy)
%   - group2_data: 3D matrix [states x states x subjects] for group 2 (MDD)
%
% Output:
%   - cohen_d_matrix: Matrix containing Cohen's d for each transition

    % Get dimensions
    transition_size = size(group1_data, 1);
    n1 = size(group1_data, 3);  % Number of subjects in group 1
    n2 = size(group2_data, 3);  % Number of subjects in group 2
    
    cohen_d_matrix = zeros(transition_size, transition_size);
    
    % Calculate Cohen's d for each transition
    for row = 1:transition_size
        for col = 1:transition_size
            if row ~= col  % Skip self-transitions
                % Extract data for this transition
                data1 = squeeze(group1_data(row, col, :));  % Group 1 data
                data2 = squeeze(group2_data(row, col, :));  % Group 2 data
                
                % Calculate means
                mean1 = mean(data1);
                mean2 = mean(data2);
                
                % Calculate pooled standard deviation
                pooled_sd = sqrt(((n1-1)*var(data1, 1) + (n2-1)*var(data2, 1)) / (n1+n2-2));
                
                % Calculate Cohen's d (avoiding division by zero)
                if pooled_sd > 0
                    cohen_d = (mean1 - mean2) / pooled_sd;
                else
                    cohen_d = 0;
                end
                
                % Store in matrix
                cohen_d_matrix(row, col) = cohen_d;
            else
                % Set diagonal elements to NaN
                cohen_d_matrix(row, col) = NaN;
            end
        end
    end
end