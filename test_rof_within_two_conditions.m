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
% - TFCE_Perm: Permuted TFCE values (the distribution of maximum TFCE per permutation)
% - P_Values: P-values from permutation test
% - Info: Structure with analysis parameters
% - Results: Structure with test results (T_Obs, Cohens_d, etc.)
%
% -------------------------------------------------------------------------
% Author:   Amin Kabir
% Date:     Jan 2025
% -------------------------------------------------------------------------

%% 1. Set defaults if not specified
if nargin < 3 || isempty(nPerm)
    nPerm = 5000;
end
if nargin < 4 || isempty(H)
    H = 2;
end
if nargin < 5 || isempty(E)
    E = 0.5;
end

%% 2. Basic checks
if size(data_matrix1,1) ~= size(data_matrix2,1) || size(data_matrix1,2) ~= size(data_matrix2,2)
    error('data_matrix1 and data_matrix2 must be the same size: [subjects x timepoints].');
end

% Number of subjects and timepoints
[nA, nTimepoints] = size(data_matrix1);

%% 3. Compute the difference for each subject
DataDiff = data_matrix1 - data_matrix2;    % [subjects x timepoints]

%% 4. Compute observed T-values (paired t-test)
%    T = mean(differences)/ (std(differences)/ sqrt(nSubjects))
disp('Calculating Observed Differences...');
T_Obs = mean(DataDiff,1) ./ (std(DataDiff,[],1) ./ sqrt(nA));

% Safety check
if max(abs(T_Obs(:))) < 1e-8
    error('All T-values are 0; there is no variability.');
end

%% 5. Set up thresholds for TFCE integration
max_t = max(abs(T_Obs));
start_threshold = 0;
threshold_step = 0.1;
thresholds = start_threshold : threshold_step : max_t;

%% 6. Calculate TFCE values for observed data
TFCE_Obs = zeros(size(T_Obs));
for thresh = thresholds
    above_thresh = abs(T_Obs) > thresh;
    if any(above_thresh)
        for t = 1:nTimepoints
            if above_thresh(t)
                cluster_mask = find_temporal_cluster(T_Obs, t);
                cluster_size = sum(cluster_mask);
                TFCE_Obs(t) = TFCE_Obs(t) + (thresh^H * cluster_size^E * threshold_step);
            end
        end
    end
end

%% 7. Permutation testing
disp('Calculating Permutations...');
maxTFCE = zeros(nPerm,1);

% Initialize waitbar
hWaitbar = waitbar(0, 'Calculating Permutations...', 'Name', 'TFCE Permutations', ...
                   'CreateCancelBtn', 'setappdata(gcbf, ''canceling'', true)');
setappdata(hWaitbar, 'canceling', false);

for i = 1:nPerm
    % Check for cancel button press
    if getappdata(hWaitbar, 'canceling')
        delete(hWaitbar);
        error('Permutation test canceled by user.');
    end

    % For each subject, randomly swap conditions with 50% probability
    % This is equivalent to sign-flipping the subject's difference.
    flip_signs = (rand(nA,1) > 0.5)*2 - 1;  % Vector of +1 or -1
    flip_signs = repmat(flip_signs, [1 nTimepoints]);
    
    % Apply sign flip
    perm_data = flip_signs .* DataDiff;  % [subjects x timepoints]
    
    % Compute permuted T-values
    T_Perm = mean(perm_data,1) ./ (std(perm_data,[],1) ./ sqrt(nA));
    
    % Compute TFCE for this permutation
    TFCE_Perm = zeros(size(T_Perm));
    for thresh = thresholds
        above_thresh = abs(T_Perm) > thresh;
        if any(above_thresh)
            for t = 1:nTimepoints
                if above_thresh(t)
                    cluster_mask = find_temporal_cluster(T_Perm, t);
                    cluster_size = sum(cluster_mask);
                    TFCE_Perm(t) = TFCE_Perm(t) + (thresh^H * cluster_size^E * threshold_step);
                end
            end
        end
    end
    
    % Store the maximum TFCE value for this permutation
    maxTFCE(i) = max(abs(TFCE_Perm));
    
    % Update progress bar
    waitbar(i / nPerm, hWaitbar, sprintf('Permutations: %d/%d', i, nPerm));
end

% Close waitbar
delete(hWaitbar);

%% 8. Calculate p-values
disp('Calculating P-Values...');
% For each timepoint, find where TFCE_Obs lies in distribution of maxTFCE
edges = [maxTFCE; max(abs(TFCE_Obs(:)))];
[~, bin] = histc(abs(TFCE_Obs), sort(edges));
P_Values = 1 - bin ./ (nPerm + 2);

%% 9. Calculate effect sizes (Cohen's d) for significant clusters (p < 0.05)
sig_clusters = P_Values < 0.05;
cohens_d = zeros(size(T_Obs));
processed = false(size(T_Obs));

for t = find(sig_clusters)
    if ~processed(t)
        cluster_mask = find_temporal_cluster(T_Obs, t);
        processed(cluster_mask) = true;
        
        % Compute effect size in this cluster
        cluster_data = mean(DataDiff(:, cluster_mask), 2);  % average difference in cluster
        d_value = mean(cluster_data) / std(cluster_data);
        cohens_d(cluster_mask) = d_value;
    end
end

%% 10. Save test information
c = clock;
Info.Comments = sprintf('TFCE paired t-test analysis run at %d:%d on %s', ...
                        c(4), c(5), date);
Info.Parameters.nPerm         = nPerm;
Info.Parameters.nTimepoints   = nTimepoints;
Info.Parameters.GroupSize     = nA;
Info.Parameters.H             = H;
Info.Parameters.E             = E;
Info.Parameters.threshold_step= threshold_step;
Info.Parameters.start_threshold= start_threshold;

%% 11. Organize results
Results.Obs         = T_Obs;
Results.TFCE_Obs    = TFCE_Obs;
Results.maxTFCE     = sort(maxTFCE);
Results.P_Values    = P_Values;
Results.Cohens_d    = cohens_d;

disp('Done!');

end
