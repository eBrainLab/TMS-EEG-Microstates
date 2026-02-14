function results = test_rtf(data_struct, stim_site, apply_correction)
% -------------------------------------------------------------------------
% TEST_RTF - Statistical analysis of relative transition frequencies
%
% This function performs statistical t-tests on post-TMS transition matrices,
% computes Cohen's d effect sizes, and applies optional p-value correction.
%
% Syntax:
%   results = test_rtf(data_struct, stim_site, apply_correction)
%
% Inputs:
%   data_struct      - Structure array containing transition matrices
%   stim_site        - String specifying stimulation site field name
%   apply_correction - String: 'fdr', 'bonferroni', or empty (no correction)
%
% Outputs:
%   results - Structure containing:
%       .h_post_tms         - Binary hypothesis test results
%       .p_post_tms         - P-values for each transition
%       .t_post_tms         - T-statistics for each transition
%       .cohen_d_post_tms   - Cohen's d effect sizes
%       .valid_subjects     - Indices of valid subjects
%       .num_valid_subjects - Count of valid subjects
%
% Method:
%   1. Validate subjects with required data
%   2. Perform one-sample t-tests on transitions
%   3. Apply multiple comparison correction
%   4. Calculate effect sizes
%
% Example:
%   rtf_results = test_rtf(data_struct, 'lpfc', 'bonferroni');
%
% See also: CALCULATE_COHENS_D, TEST_RTF_WITHIN_TWO_CONDITIONS
% -------------------------------------------------------------------------

%% Validate input data
fprintf('Validating subjects for RTF analysis...\n');

valid_subjects = [];
for kk = 1:length(data_struct)
    % Check for required fields
    if isfield(data_struct(kk), stim_site) && ...
       isfield(data_struct(kk).(stim_site), 'transition_averages_bc')
        valid_subjects = [valid_subjects, kk];
    end
end

% Check if any valid subjects found
if isempty(valid_subjects)
    error(['No subjects found with transition_averages_bc data for stim_site "', ...
           stim_site, '".']);
end

fprintf('Found %d valid subjects out of %d total for stim_site "%s".\n', ...
    length(valid_subjects), length(data_struct), stim_site);

%% Extract transition data
% Get dimensions from first valid subject
first_valid = valid_subjects(1);
matrix_size = size(data_struct(first_valid).(stim_site).transition_averages_bc.post_tms);
num_valid_subjects = length(valid_subjects);

% Preallocate transition array
post_tms_sample = zeros(matrix_size(1), matrix_size(2), num_valid_subjects);

% Extract data for valid subjects
for idx = 1:num_valid_subjects
    kk = valid_subjects(idx);
    post_tms_sample(:, :, idx) = data_struct(kk).(stim_site).transition_averages_bc.post_tms;
end

%% Perform statistical testing
fprintf('Performing t-tests on transition matrices...\n');

% One-sample t-test against zero
if exist('ttest', 'file')
    [~, p_post_tms, ~, stats] = ttest(post_tms_sample, 0, 'Dim', 3);
    t_post_tms = stats.tstat;
else
    % Fallback: manual one-sample t-test without Statistics Toolbox
    fprintf('Note: Statistics Toolbox not found. Using built-in t-test implementation.\n');
    n_subj = size(post_tms_sample, 3);
    m = mean(post_tms_sample, 3);
    s = std(post_tms_sample, 0, 3);
    se = s ./ sqrt(n_subj);
    t_post_tms = m ./ se;
    df = n_subj - 1;
    % Two-tailed p-value via regularized incomplete beta function (base MATLAB)
    p_post_tms = betainc(df ./ (df + t_post_tms.^2), df/2, 0.5);
    p_post_tms(isnan(t_post_tms)) = NaN;
end

% Set diagonal (self-transitions) to NaN
for ii = 1:size(p_post_tms, 1)
    p_post_tms(ii, ii) = NaN;
end

%% Apply multiple comparison correction
if ~isempty(apply_correction)
    fprintf('Applying %s correction...\n', apply_correction);
    
    % Flatten p-values, excluding NaNs
    p_post_tms_flat = p_post_tms(~isnan(p_post_tms));
    
    if strcmpi(apply_correction, "fdr")
        % False Discovery Rate correction
        if exist('mafdr', 'file')
            p_post_tms_corrected = mafdr(p_post_tms_flat, 'BHFDR', true);
        else
            % Fallback: Benjamini-Hochberg FDR without Bioinformatics Toolbox
            fprintf('Note: Bioinformatics Toolbox not found. Using built-in BH-FDR correction.\n');
            m_tests = length(p_post_tms_flat);
            [p_sorted, sort_idx] = sort(p_post_tms_flat(:));
            adjusted = p_sorted .* m_tests ./ (1:m_tests)';
            % Enforce monotonicity (step-up procedure)
            for adj_i = m_tests-1:-1:1
                adjusted(adj_i) = min(adjusted(adj_i), adjusted(adj_i+1));
            end
            p_post_tms_corrected = min(adjusted, 1);
            % Restore original order
            p_temp = zeros(size(p_post_tms_flat));
            p_temp(sort_idx) = p_post_tms_corrected;
            p_post_tms_corrected = p_temp;
        end
        
    elseif strcmpi(apply_correction, "bonferroni")
        % Bonferroni correction
        num_tests = numel(p_post_tms_flat);
        p_post_tms_corrected = min(p_post_tms_flat * num_tests, 1);
        
    else
        warning('Unknown correction method "%s". No correction applied.', ...
                apply_correction);
        p_post_tms_corrected = p_post_tms_flat;
    end
    
    % Reshape corrected p-values back to matrix
    p_post_tms(~isnan(p_post_tms)) = p_post_tms_corrected;
end

%% Determine significant transitions
alpha = 0.05;
h_post_tms = p_post_tms < alpha;

%% Calculate effect sizes
fprintf('Computing Cohen''s d effect sizes...\n');
cohen_d_post_tms = calculate_cohens_d(post_tms_sample);

%% Compile results
results.h_post_tms = h_post_tms;
results.p_post_tms = p_post_tms;
results.t_post_tms = t_post_tms;
results.cohen_d_post_tms = cohen_d_post_tms;
results.valid_subjects = valid_subjects;
results.num_valid_subjects = num_valid_subjects;

fprintf('RTF analysis complete.\n');

end