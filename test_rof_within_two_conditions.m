function [TFCE_Obs, TFCE_Perm, P_Values, Info, Results] = test_rof_within_two_conditions(...
    data_matrix1, data_matrix2, nPerm, H, E, progressFcn)
% -------------------------------------------------------------------------
% TEST_ROF_WITHIN_TWO_CONDITIONS - TFCE for paired within-subject design
%
% This function computes TFCE statistics for a two-condition paired design,
% using condition swapping for permutation testing.
%
% Syntax:
%   [TFCE_Obs, TFCE_Perm, P_Values, Info, Results] = ...
%       test_rof_within_two_conditions(data_matrix1, data_matrix2, nPerm, H, E)
%   [... ] = test_rof_within_two_conditions(..., progressFcn)
%
% Inputs:
%   data_matrix1 - Time series data [subjects x timepoints] for condition 1
%   data_matrix2 - Time series data [subjects x timepoints] for condition 2
%   nPerm        - Number of permutations (default: 5000)
%   H            - TFCE height parameter (default: 2)
%   E            - TFCE extent parameter (default: 0.5)
%   progressFcn  - (Optional) Function handle called as progressFcn(fraction, message)
%                  When provided, suppresses built-in waitbar and fprintf output.
%
% Outputs:
%   TFCE_Obs     - Observed TFCE values [1 x timepoints]
%   TFCE_Perm    - Max TFCE distribution from permutations [nPerm x 1]
%   P_Values     - P-values from permutation test [1 x timepoints]
%   Info         - Structure with analysis parameters
%   Results      - Structure with complete results including Cohen's d
%
% Method:
%   1. Calculate paired differences
%   2. Compute paired t-values
%   3. Apply TFCE integration
%   4. Generate null distribution via condition swapping
%   5. Compute p-values and effect sizes
%
% Example:
%   [TFCE_Obs, TFCE_Perm, P_Values, Info, Results] = ...
%       test_rof_within_two_conditions(data_site1, data_site2, 5000);
%
% See also: TEST_ROF_ONE_CONDITION, FIND_TEMPORAL_CLUSTER
% -------------------------------------------------------------------------

%% Set default parameters
if nargin < 3 || isempty(nPerm)
    nPerm = 5000;
end
if nargin < 4 || isempty(H)
    H = 2;
end
if nargin < 5 || isempty(E)
    E = 0.5;
end

% Determine output mode
useExternalProgress = nargin >= 6 && ~isempty(progressFcn) && isa(progressFcn, 'function_handle');
if useExternalProgress
    logfcn = @(varargin) [];  % suppress fprintf
else
    logfcn = @(varargin) fprintf(varargin{:});
end

%% Validate inputs
if size(data_matrix1,1) ~= size(data_matrix2,1) || ...
   size(data_matrix1,2) ~= size(data_matrix2,2)
    error('data_matrix1 and data_matrix2 must have the same dimensions.');
end

% Get dimensions
[nA, nTimepoints] = size(data_matrix1);

%% Calculate paired differences
logfcn('Calculating paired differences...\n');
DataDiff = data_matrix1 - data_matrix2;

%% Calculate observed paired t-values
T_Obs = mean(DataDiff, 1) ./ (std(DataDiff, [], 1) ./ sqrt(nA));

% Check for valid T-values
if max(abs(T_Obs(:))) < 1e-8
    error('All T-values are zero; no variability in differences.');
end

%% TFCE setup
max_t = max(abs(T_Obs));
start_threshold = 0;
threshold_step = 0.1;
thresholds = start_threshold : threshold_step : max_t;

%% TFCE integration on observed data
logfcn('Computing TFCE on observed data...\n');
TFCE_Obs = zeros(size(T_Obs));

for thresh = thresholds
    above_thresh = abs(T_Obs) > thresh;
    
    if any(above_thresh)
        for t = 1:nTimepoints
            if above_thresh(t)
                % Find temporal cluster
                cluster_mask = find_temporal_cluster(T_Obs, t);
                cluster_size = sum(cluster_mask);
                
                % Update TFCE value
                TFCE_Obs(t) = TFCE_Obs(t) + ...
                    (thresh^H * cluster_size^E * threshold_step);
            end
        end
    end
end

%% Permutation testing
logfcn('Starting permutation testing...\n');
maxTFCE = zeros(nPerm, 1);

% Create progress bar (only when no external progress callback)
if ~useExternalProgress
    hWaitbar = waitbar(0, 'Calculating Permutations...', ...
        'Name', 'TFCE Permutations', ...
        'CreateCancelBtn', 'setappdata(gcbf, ''canceling'', true)');
    setappdata(hWaitbar, 'canceling', false);
end

for i = 1:nPerm
    % Check for cancellation
    if ~useExternalProgress
        if getappdata(hWaitbar, 'canceling')
            delete(hWaitbar);
            error('Permutation test canceled by user.');
        end
    end
    
    % Randomly swap conditions (equivalent to sign flipping differences)
    flip_signs = (rand(nA,1) > 0.5)*2 - 1;  % Vector of +1 or -1
    flip_signs = repmat(flip_signs, [1 nTimepoints]);
    
    % Apply sign flip to differences
    perm_data = flip_signs .* DataDiff;
    
    % Calculate permuted t-values
    T_Perm = mean(perm_data, 1) ./ (std(perm_data, [], 1) ./ sqrt(nA));
    
    % TFCE for permuted data
    temp_TFCE = zeros(size(T_Perm));
    for thresh = thresholds
        above_thresh = abs(T_Perm) > thresh;
        
        if any(above_thresh)
            for t = 1:nTimepoints
                if above_thresh(t)
                    % Find temporal cluster
                    cluster_mask = find_temporal_cluster(T_Perm, t);
                    cluster_size = sum(cluster_mask);
                    
                    % Update TFCE value
                    temp_TFCE(t) = temp_TFCE(t) + ...
                        (thresh^H * cluster_size^E * threshold_step);
                end
            end
        end
    end
    
    % Store maximum TFCE value
    maxTFCE(i) = max(abs(temp_TFCE));
    
    % Update progress
    if useExternalProgress
        progressFcn(i / nPerm, sprintf('Permutations: %d/%d', i, nPerm));
    else
        waitbar(i / nPerm, hWaitbar, sprintf('Permutations: %d/%d', i, nPerm));
    end
end

if ~useExternalProgress
    delete(hWaitbar);
end
logfcn('Permutation testing complete.\n');

%% Calculate p-values
logfcn('Calculating p-values and effect sizes...\n');

if exist('histc', 'file') || exist('histc', 'builtin')
    % Use histogram method for p-values
    edges = [maxTFCE; max(abs(TFCE_Obs(:)))];
    [~, bin] = histc(abs(TFCE_Obs), sort(edges));
    P_Values = 1 - bin ./ (nPerm + 2);
else
    % Fallback: direct permutation p-value computation (histc removed)
    logfcn('Note: histc not available. Using direct p-value computation.\n');
    P_Values = zeros(size(TFCE_Obs));
    for t_idx = 1:length(TFCE_Obs)
        P_Values(t_idx) = (sum(maxTFCE >= abs(TFCE_Obs(t_idx))) + 1) / (nPerm + 1);
    end
end

%% Identify significant clusters and compute effect sizes
% Find significant time points
sig_mask = (P_Values < 0.05);

% Initialize effect size arrays
cohens_d = zeros(size(T_Obs));
processed_timepoints = false(size(T_Obs));
cluster_d_values = [];

% Process each significant cluster
for t = 1:nTimepoints
    if sig_mask(t) && ~processed_timepoints(t)
        % Find cluster restricted to significant points
        cluster_mask = find_temporal_cluster(T_Obs, t, sig_mask);
        
        % Mark as processed
        processed_timepoints(cluster_mask) = true;
        
        % Calculate cluster-level Cohen's d
        cluster_data = mean(DataDiff(:, cluster_mask), 2);
        d_value = mean(cluster_data) / std(cluster_data);
        
        % Assign to all points in cluster
        cohens_d(cluster_mask) = d_value;
        
        % Store for averaging
        cluster_d_values(end+1) = d_value;
    end
end

% Calculate average effect size across clusters
if ~isempty(cluster_d_values)
    avg_cohens_d = mean(cluster_d_values);
else
    avg_cohens_d = NaN;
end

% Sort permutation distribution
TFCE_Perm = sort(maxTFCE);

%% Compile results
% Analysis information
c = clock;
Info.Comments = sprintf('TFCE paired t-test analysis conducted at %02d:%02d on %s', ...
    c(4), c(5), date);
Info.Parameters.nPerm = nPerm;
Info.Parameters.nTimepoints = nTimepoints;
Info.Parameters.GroupSize = nA;
Info.Parameters.H = H;
Info.Parameters.E = E;
Info.Parameters.threshold_step = threshold_step;
Info.Parameters.start_threshold = start_threshold;

% Results structure
Results.Obs = T_Obs;
Results.TFCE_Obs = TFCE_Obs;
Results.maxTFCE = TFCE_Perm;
Results.P_Values = P_Values;
Results.Cohens_d = cohens_d;
Results.Avg_Cohens_d = avg_cohens_d;

logfcn('Analysis complete.\n');

end