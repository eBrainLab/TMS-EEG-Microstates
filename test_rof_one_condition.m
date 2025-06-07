function [TFCE_Obs, TFCE_Perm, P_Values, Info, Results] = test_rof_one_condition(...
    data_matrix, nPerm, H, E)
% -------------------------------------------------------------------------
% TEST_ROF_ONE_CONDITION - TFCE analysis for single-condition time series
%
% This function computes Threshold-Free Cluster Enhancement (TFCE) statistics
% for one-condition time series data using sign-flipping permutations.
%
% Syntax:
%   [TFCE_Obs, TFCE_Perm, P_Values, Info, Results] = test_rof_one_condition(...
%       data_matrix, nPerm, H, E)
%
% Inputs:
%   data_matrix - Time series data [subjects x timepoints]
%   nPerm       - Number of permutations (default: 5000)
%   H           - TFCE height parameter (default: 2)
%   E           - TFCE extent parameter (default: 0.5)
%
% Outputs:
%   TFCE_Obs    - Observed TFCE values [1 x timepoints]
%   TFCE_Perm   - Max TFCE distribution from permutations [nPerm x 1]
%   P_Values    - P-values from permutation test [1 x timepoints]
%   Info        - Structure with analysis parameters
%   Results     - Structure with complete results including Cohen's d
%
% Method:
%   1. Calculate one-sample t-values
%   2. Apply TFCE integration
%   3. Generate null distribution via sign-flipping
%   4. Compute p-values and effect sizes
%
% Example:
%   [TFCE_Obs, TFCE_Perm, P_Values, Info, Results] = ...
%       test_rof_one_condition(data, 5000);
%
% See also: FIND_TEMPORAL_CLUSTER, CALCULATE_COHENS_D
% -------------------------------------------------------------------------

%% Set default parameters
if nargin < 2 || isempty(nPerm)
    nPerm = 5000;
end
if nargin < 3 || isempty(H)
    H = 2;
end
if nargin < 4 || isempty(E)
    E = 0.5;
end

%% Initialize
% Store data in cell for compatibility
Data{1} = data_matrix;
[nA, nTimepoints] = size(Data{1});

% Preallocate
maxTFCE = zeros(nPerm, 1);

%% Calculate observed T-values
fprintf('Calculating observed T-values...\n');
T_Obs = mean(Data{1}, 1) ./ (std(Data{1}) ./ sqrt(nA));

% Check for valid T-values
if max(abs(T_Obs(:))) < 1e-5
    error('T-values are all zero. Check input data.');
end

%% TFCE integration on observed data
% Set up thresholds
max_t = max(abs(T_Obs));
start_threshold = 0;
threshold_step = 0.1;
thresholds = start_threshold : threshold_step : max_t;

% Initialize TFCE values
TFCE_Obs = zeros(size(T_Obs));

% Integrate over thresholds
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
fprintf('TFCE integration complete.\n');

%% Permutation testing
fprintf('Starting permutation testing...\n');

% Create progress bar
hWaitbar = waitbar(0, 'Calculating Permutations...', ...
    'Name', 'TFCE Permutations', ...
    'CreateCancelBtn', 'setappdata(gcbf, ''canceling'', true)');
setappdata(hWaitbar, 'canceling', false);

for i = 1:nPerm
    % Check for cancellation
    if getappdata(hWaitbar, 'canceling')
        delete(hWaitbar);
        error('Permutation test canceled by user.');
    end
    
    % Random sign flipping
    Signs = [-1, 1];
    SignSwitch = randsample(Signs, nA, 'true')';
    SignSwitch = repmat(SignSwitch, [1, nTimepoints]);
    
    % Apply sign flip
    nData = SignSwitch .* Data{1};
    
    % Calculate permuted T-values
    T_Perm = mean(nData, 1) ./ (std(nData) ./ sqrt(nA));
    
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
    waitbar(i / nPerm, hWaitbar, sprintf('Permutations: %d/%d', i, nPerm));
end

delete(hWaitbar);
fprintf('Permutation testing complete.\n');

%% Calculate p-values
fprintf('Calculating p-values and effect sizes...\n');

% Use histogram method for p-values
edges = [maxTFCE; max(abs(TFCE_Obs(:)))];
[~, bin] = histc(abs(TFCE_Obs), sort(edges));
P_Values = 1 - bin ./ (nPerm + 2);

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
        this_cluster = find_temporal_cluster(TFCE_Obs, t, sig_mask);
        
        % Mark as processed
        processed_timepoints(this_cluster) = true;
        
        % Calculate cluster-level Cohen's d
        cluster_data = mean(Data{1}(:, this_cluster), 2);
        d_value = mean(cluster_data) / std(cluster_data);
        
        % Assign to all points in cluster
        cohens_d(this_cluster) = d_value;
        
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
Info.Comments = sprintf('TFCE analysis conducted at %02d:%02d on %s', ...
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
Results.TFCE_Perm = TFCE_Perm;
Results.P_Values = P_Values;
Results.Cohens_d = cohens_d;
Results.Avg_Cohens_d = avg_cohens_d;

fprintf('Analysis complete.\n');

end