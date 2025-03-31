function [TFCE_Obs, TFCE_Perm, P_Values, Info, Results] = test_rof_within_two_conditions(data_matrix1, data_matrix2, nPerm, H, E)
% Computes TFCE statistics for a two-condition (paired) within-subject design
% using swapping of conditions for permutations.
%
% Inputs:
% - data_matrix1: Time series data matrix [subjects x timepoints] for condition 1
% - data_matrix2: Time series data matrix [subjects x timepoints] for condition 2
% - nPerm: Number of permutations (default 5000)
% - H: TFCE height parameter (default 2)
% - E: TFCE extent parameter (default 0.5)
%
% Outputs:
% - TFCE_Obs: Observed TFCE values
% - TFCE_Perm: Distribution of max TFCE values across permutations (useful for thresholding)
% - P_Values: P-values (per timepoint) from permutation test
% - Info: Structure with analysis parameters
% - Results: Structure with T-values, Cohen’s d, etc.
%
% -------------------------------------------------------------------------
% Author:   Amin Kabir
% Date:     Jan 2025
% -------------------------------------------------------------------------

%% 1) Set defaults if not specified
if nargin < 3 || isempty(nPerm)
    nPerm = 5000;
end
if nargin < 4 || isempty(H)
    H = 2;
end
if nargin < 5 || isempty(E)
    E = 0.5;
end

%% 2) Basic checks
if size(data_matrix1,1) ~= size(data_matrix2,1) || ...
   size(data_matrix1,2) ~= size(data_matrix2,2)
    error('data_matrix1 and data_matrix2 must be [subjects x timepoints] of the same size.');
end

% Number of subjects and timepoints
[nA, nTimepoints] = size(data_matrix1);

%% 3) Compute the within-subject differences
DataDiff = data_matrix1 - data_matrix2; % [subjects x timepoints]

%% 4) Observed T-values for paired t-test
disp('Calculating Observed Differences...');
T_Obs = mean(DataDiff, 1) ./ (std(DataDiff, [], 1) ./ sqrt(nA));

if max(abs(T_Obs(:))) < 1e-8
    error('All T-values are 0; there is no variability.');
end

%% 5) TFCE setup
max_t         = max(abs(T_Obs));
start_threshold = 0;
threshold_step  = 0.1;
thresholds      = start_threshold : threshold_step : max_t;

%% 6) TFCE on observed data
TFCE_Obs = zeros(size(T_Obs));
for thresh = thresholds
    above_thresh = abs(T_Obs) > thresh;
    if any(above_thresh)
        for t = 1:nTimepoints
            if above_thresh(t)
                % Find contiguous timepoints of same sign
                cluster_mask = find_temporal_cluster(T_Obs, t);
                cluster_size = sum(cluster_mask);
                TFCE_Obs(t) = TFCE_Obs(t) + ...
                              (thresh^H * cluster_size^E * threshold_step);
            end
        end
    end
end

%% 7) Permutation testing
disp('Calculating Permutations...');
maxTFCE = zeros(nPerm,1);

hWaitbar = waitbar(0, 'Calculating Permutations...', ...
    'Name', 'TFCE Permutations', ...
    'CreateCancelBtn', 'setappdata(gcbf, ''canceling'', true)');
setappdata(hWaitbar, 'canceling', false);

for i = 1:nPerm
    % Check for cancel button press
    if getappdata(hWaitbar, 'canceling')
        delete(hWaitbar);
        error('Permutation test canceled by user.');
    end

    % Randomly swap conditions (equivalent to sign flipping for the differences)
    flip_signs = (rand(nA,1) > 0.5)*2 - 1; % Vector of +1 or -1
    flip_signs = repmat(flip_signs, [1 nTimepoints]);
    
    % Apply sign flip
    perm_data = flip_signs .* DataDiff;  % [subjects x timepoints]
    
    % Permuted T-values
    T_Perm = mean(perm_data,1) ./ (std(perm_data,[],1) ./ sqrt(nA));
    
    % Compute TFCE for this permutation
    temp_TFCE = zeros(size(T_Perm));
    for thresh = thresholds
        above_thresh = abs(T_Perm) > thresh;
        if any(above_thresh)
            for t = 1:nTimepoints
                if above_thresh(t)
                    cluster_mask = find_temporal_cluster(T_Perm, t);
                    cluster_size = sum(cluster_mask);
                    temp_TFCE(t) = temp_TFCE(t) + ...
                                   (thresh^H * cluster_size^E * threshold_step);
                end
            end
        end
    end
    
    % Store the maximum TFCE value
    maxTFCE(i) = max(abs(temp_TFCE));
    
    waitbar(i / nPerm, hWaitbar, sprintf('Permutations: %d/%d', i, nPerm));
end

delete(hWaitbar);

%% 8) P-values from permutation distribution
disp('Calculating P-Values...');
edges      = [maxTFCE; max(abs(TFCE_Obs(:)))];
[~, bin]   = histc(abs(TFCE_Obs), sort(edges));
P_Values   = 1 - bin ./ (nPerm + 2);

%% 9) Effect sizes (Cohen's d) for significant clusters
sig_mask             = (P_Values < 0.05);
cohens_d             = zeros(size(T_Obs));
processed_timepoints = false(size(T_Obs));
cluster_d_values     = [];

for t = 1:nTimepoints
    if sig_mask(t) && ~processed_timepoints(t)
        % Find contiguous cluster of same sign
        cluster_mask = find_temporal_cluster(T_Obs, t);
        % Restrict to only significant timepoints
        cluster_mask = cluster_mask & sig_mask;
        
        processed_timepoints(cluster_mask) = true;
        
        % Extract data and compute cluster-level effect size
        cluster_data = mean(DataDiff(:, cluster_mask), 2);
        d_value      = mean(cluster_data) / std(cluster_data);
        
        % Assign that same effect size to all points in the cluster
        cohens_d(cluster_mask) = d_value;
        
        % Keep track for averaging (if desired)
        cluster_d_values(end+1) = d_value;
    end
end

% Compute average Cohen's d across all significant clusters
if ~isempty(cluster_d_values)
    avg_cohens_d = mean(cluster_d_values);
else
    avg_cohens_d = NaN;
end

%% 10) Assign output for TFCE_Perm
% This is the distribution of maximum TFCE values across all permutations.
TFCE_Perm = sort(maxTFCE);

%% 11) Save test information
c = clock;
Info.Comments = sprintf('TFCE paired t-test analysis run at %02d:%02d on %s', ...
                        c(4), c(5), date);
Info.Parameters.nPerm           = nPerm;
Info.Parameters.nTimepoints     = nTimepoints;
Info.Parameters.GroupSize       = nA;
Info.Parameters.H               = H;
Info.Parameters.E               = E;
Info.Parameters.threshold_step  = threshold_step;
Info.Parameters.start_threshold = start_threshold;

%% 12) Organize results
Results.Obs          = T_Obs;
Results.TFCE_Obs     = TFCE_Obs;
Results.maxTFCE      = TFCE_Perm;  % sorted distribution of max TFCE
Results.P_Values     = P_Values;
Results.Cohens_d     = cohens_d;   % one Cohen's d value per timepoint
Results.Avg_Cohens_d = avg_cohens_d;

disp('Done!');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        HELPER FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cluster_mask = find_temporal_cluster(T_values, time_idx)
% Finds contiguous timepoints around 'time_idx' that share the same sign.
% Returns a logical mask with true at cluster points, false otherwise.

cluster_mask         = false(size(T_values));
cluster_mask(time_idx) = true;

% Expand forward
t = time_idx;
while t < length(T_values) && sign(T_values(t+1)) == sign(T_values(t))
    cluster_mask(t+1) = true;
    t = t + 1;
end

% Expand backward
t = time_idx;
while t > 1 && sign(T_values(t-1)) == sign(T_values(t))
    cluster_mask(t-1) = true;
    t = t - 1;
end
end
