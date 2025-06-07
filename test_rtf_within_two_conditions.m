function results = test_rtf_within_two_conditions(data_struct, stim_sites, apply_correction)
% -------------------------------------------------------------------------
% TEST_RTF_WITHIN_TWO_CONDITIONS - Paired comparison of transition frequencies
%
% This function performs statistical analysis to evaluate differences in
% transition dynamics between two stimulation conditions using paired tests.
%
% Syntax:
%   results = test_rtf_within_two_conditions(data_struct, stim_sites, apply_correction)
%
% Inputs:
%   data_struct      - Structure array containing transition matrices
%   stim_sites       - Array of two stimulation site names to compare
%   apply_correction - String: 'fdr', 'bonferroni', or empty (no correction)
%
% Outputs:
%   results - Structure containing:
%       .h_post_tms       - Binary hypothesis test results
%       .p_post_tms       - P-values for each transition
%       .t_post_tms       - T-statistics for each transition  
%       .cohen_d_post_tms - Cohen's d effect sizes
%
% Method:
%   1. Calculate differences between conditions
%   2. Perform one-sample t-tests on differences
%   3. Apply multiple comparison correction
%   4. Calculate effect sizes
%
% Example:
%   rtf_results = test_rtf_within_two_conditions(data_struct, ["lpfc","rpfc"], 'fdr');
%
% See also: TEST_RTF, CALCULATE_COHENS_D
% -------------------------------------------------------------------------

%% Initialize
fprintf('Starting paired RTF comparison between %s and %s...\n', ...
    upper(stim_sites(1)), upper(stim_sites(2)));

% Get dimensions
num_states = size(data_struct(1).(stim_sites(1)).transition_averages_bc.post_tms, 1);
num_subjects = length(data_struct);

% Preallocate difference matrices
difference_post_tms = zeros(num_states, num_states, num_subjects);

%% Calculate paired differences
fprintf('Computing paired differences...\n');

for sub = 1:num_subjects
    % Calculate difference: primary site - secondary site
    difference_post_tms(:, :, sub) = ...
        data_struct(sub).(stim_sites(1)).transition_averages_bc.post_tms - ...
        data_struct(sub).(stim_sites(2)).transition_averages_bc.post_tms;
end

%% Perform statistical testing
fprintf('Performing paired t-tests...\n');

% One-sample t-test on differences (testing if mean difference = 0)
[~, p_post_tms, ~, stats_post_tms] = ttest(difference_post_tms, 0, 'Dim', 3, 'Tail', 'both');
t_post_tms = stats_post_tms.tstat;

% Set diagonal (self-transitions) to NaN
for ii = 1:num_states
    p_post_tms(ii, ii) = NaN;
end

%% Apply multiple comparison correction
if nargin >= 3 && ~isempty(apply_correction)
    fprintf('Applying %s correction...\n', apply_correction);
    
    % Flatten p-values, excluding NaNs
    p_post_tms_flat = p_post_tms(~isnan(p_post_tms));
    
    if strcmpi(apply_correction, "fdr")
        % False Discovery Rate correction
        p_post_tms_corrected = mafdr(p_post_tms_flat, 'BHFDR', true);
        
    elseif strcmpi(apply_correction, "bonferroni")
        % Bonferroni correction
        num_tests_post_tms = numel(p_post_tms_flat);
        p_post_tms_corrected = min(p_post_tms_flat * num_tests_post_tms, 1);
        
    else
        error('Unsupported correction method. Use "fdr", "bonferroni", or leave empty.');
    end
    
    % Reshape corrected p-values back to matrix
    p_post_tms(~isnan(p_post_tms)) = p_post_tms_corrected;
end

%% Determine significant transitions
alpha = 0.05;
h_post_tms = p_post_tms < alpha;

%% Calculate effect sizes
fprintf('Computing Cohen''s d effect sizes...\n');
cohen_d_post_tms = calculate_cohens_d(difference_post_tms);

%% Compile results
results.h_post_tms = h_post_tms;
results.p_post_tms = p_post_tms;
results.t_post_tms = t_post_tms;
results.cohen_d_post_tms = cohen_d_post_tms;

fprintf('Paired RTF analysis complete.\n');

end