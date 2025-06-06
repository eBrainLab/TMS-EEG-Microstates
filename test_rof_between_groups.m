function [TFCE_Obs, TFCE_Perm, P_Values, Info, Results] = test_rof_between_groups(data_matrix1, data_matrix2, nPerm, H, E)
% Computes TFCE statistics for a two-group (unpaired) between-subjects design
% using permutation testing.
%
% Inputs:
% - data_matrix1: Time series data matrix [subjects x timepoints] for group 1 (healthy)
% - data_matrix2: Time series data matrix [subjects x timepoints] for group 2 (MDD)
% - nPerm: Number of permutations (default 5000)
% - H: TFCE height parameter (default 2)
% - E: TFCE extent parameter (default 0.5)
%
% Outputs:
% - TFCE_Obs: Observed TFCE values
% - TFCE_Perm: Distribution of max TFCE values across permutations (useful for thresholding)
% - P_Values: P-values (per timepoint) from permutation test
% - Info: Structure with analysis parameters
% - Results: Structure with T-values, Cohen's d, etc.
%
% -------------------------------------------------------------------------
% Author: Based on code by Amin Kabir, adapted for between-groups analysis
% Date: April 2025
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
if size(data_matrix1,2) ~= size(data_matrix2,2)
    error('data_matrix1 and data_matrix2 must have the same number of timepoints.');
end

% Number of subjects in each group and timepoints
[nA, nTimepoints] = size(data_matrix1);
[nB, ~] = size(data_matrix2);

fprintf('Group 1: %d subjects, Group 2: %d subjects\n', nA, nB);

%% 3) Compute observed independent samples t-test
disp('Calculating Observed Between-Group Differences...');

% Calculate means for each group
mean_1 = mean(data_matrix1, 1);
mean_2 = mean(data_matrix2, 1);

% Calculate standard deviations for each group
std_1 = std(data_matrix1, [], 1);
std_2 = std(data_matrix2, [], 1);

% Calculate pooled standard deviation
s_pooled = sqrt(((nA-1).*std_1.^2 + (nB-1).*std_2.^2) ./ (nA + nB - 2));

% Calculate t-statistic for independent samples t-test
T_Obs = (mean_1 - mean_2) ./ (s_pooled .* sqrt(1/nA + 1/nB));

if max(abs(T_Obs(:))) < 1e-8
    error('All T-values are 0; there is no variability or difference between groups.');
end

%% 4) TFCE setup
max_t = max(abs(T_Obs));
start_threshold = 0;
threshold_step = 0.1;
thresholds = start_threshold : threshold_step : max_t;

%% 5) TFCE on observed data
TFCE_Obs = zeros(size(T_Obs));
for thresh = thresholds
    above_thresh = abs(T_Obs) > thresh;
    if any(above_thresh)
        for t = 1:nTimepoints
            if above_thresh(t)
                % Find contiguous timepoints of same sign
                cluster_mask = find_temporal_cluster(T_Obs, t);
                cluster_size = sum(cluster_mask);
                TFCE_Obs(t) = TFCE_Obs(t) + (thresh^H * cluster_size^E * threshold_step);
            end
        end
    end
end

%% 6) Permutation testing
disp('Calculating Permutations...');
maxTFCE = zeros(nPerm, 1);

hWaitbar = waitbar(0, 'Calculating Permutations...', ...
    'Name', 'TFCE Permutations', ...
    'CreateCancelBtn', 'setappdata(gcbf, ''canceling'', true)');
setappdata(hWaitbar, 'canceling', false);

% Combine all subjects from both groups for permutation
all_data = [data_matrix1; data_matrix2];
total_subjects = nA + nB;

for i = 1:nPerm
    % Check for cancel button press
    if getappdata(hWaitbar, 'canceling')
        delete(hWaitbar);
        error('Permutation test canceled by user.');
    end
    
    % Random permutation of subject indices
    perm_indices = randperm(total_subjects);
    
    % Assign random subjects to group 1 and group 2 (preserving group sizes)
    perm_group1 = all_data(perm_indices(1:nA), :);
    perm_group2 = all_data(perm_indices(nA+1:end), :);
    
    % Calculate permuted statistics
    mean_perm1 = mean(perm_group1, 1);
    mean_perm2 = mean(perm_group2, 1);
    
    std_perm1 = std(perm_group1, [], 1);
    std_perm2 = std(perm_group2, [], 1);
    
    s_pooled_perm = sqrt(((nA-1).*std_perm1.^2 + (nB-1).*std_perm2.^2) ./ (nA + nB - 2));
    
    % Calculate t-statistic for permuted groups
    T_Perm = (mean_perm1 - mean_perm2) ./ (s_pooled_perm .* sqrt(1/nA + 1/nB));
    
    % Compute TFCE for this permutation
    temp_TFCE = zeros(size(T_Perm));
    for thresh = thresholds
        above_thresh = abs(T_Perm) > thresh;
        if any(above_thresh)
            for t = 1:nTimepoints
                if above_thresh(t)
                    cluster_mask = find_temporal_cluster(T_Perm, t);
                    cluster_size = sum(cluster_mask);
                    temp_TFCE(t) = temp_TFCE(t) + (thresh^H * cluster_size^E * threshold_step);
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
edges = [maxTFCE; max(abs(TFCE_Obs(:)))];
[~, bin] = histc(abs(TFCE_Obs), sort(edges));
P_Values = 1 - bin ./ (nPerm + 2);

%% 9) Effect sizes (Cohen's d) for significant clusters
sig_mask = (P_Values < 0.05);
cohens_d = zeros(size(T_Obs));
processed_timepoints = false(size(T_Obs));
cluster_d_values = [];

for t = 1:nTimepoints
    if sig_mask(t) && ~processed_timepoints(t)
        % Find the contiguous cluster around this timepoint
        cluster_mask = find_temporal_cluster_significant(T_Obs, sig_mask, t);
        processed_timepoints(cluster_mask) = true;
        
        % Calculate mean values for these timepoints for each group
        mean_cluster_1 = mean(mean(data_matrix1(:, cluster_mask), 2));
        mean_cluster_2 = mean(mean(data_matrix2(:, cluster_mask), 2));
        
        % Calculate pooled standard deviation for this cluster
        pooled_std = sqrt(((nA-1)*var(mean(data_matrix1(:, cluster_mask), 2)) + ...
                          (nB-1)*var(mean(data_matrix2(:, cluster_mask), 2))) / ...
                          (nA + nB - 2));
                      
        % Calculate Cohen's d for between-groups comparison
        if pooled_std > 0
            d_value = (mean_cluster_1 - mean_cluster_2) / pooled_std;
        else
            d_value = 0;
        end
        
        % Assign to all timepoints in this cluster
        cohens_d(cluster_mask) = d_value;
        
        % Keep track of cluster-level effect sizes
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
Info.Comments = sprintf('TFCE independent t-test analysis run at %02d:%02d on %s', ...
                       c(4), c(5), date);
Info.Parameters.nPerm = nPerm;
Info.Parameters.nTimepoints = nTimepoints;
Info.Parameters.Group1Size = nA;
Info.Parameters.Group2Size = nB;
Info.Parameters.H = H;
Info.Parameters.E = E;
Info.Parameters.threshold_step = threshold_step;
Info.Parameters.start_threshold = start_threshold;

%% 12) Organize results
Results.Obs = T_Obs;
Results.TFCE_Obs = TFCE_Obs;
Results.maxTFCE = TFCE_Perm;  % sorted distribution of max TFCE
Results.P_Values = P_Values;
Results.Cohens_d = cohens_d;   % one Cohen's d value per timepoint
Results.Avg_Cohens_d = avg_cohens_d;

disp('TFCE analysis complete!');
end