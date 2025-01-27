function [TFCE_Obs, TFCE_Perm, P_Values, Info, Results] = test_rof_one_condition(data_matrix, nPerm, H, E)
% Computes TFCE statistics for one condition time series data using sign-flipping permutations
%
% Inputs:
% - data_matrix: Time series data matrix [subjects x timepoints]
% - nPerm: Number of permutations (default 5000)
% - H: TFCE height parameter (default 2)
% - E: TFCE extent parameter (default 0.5)
%
% Outputs:
% - TFCE_Obs: Observed TFCE values
% - TFCE_Perm: Permuted TFCE values
% - P_Values: P-values from permutation test
% - Info: Structure with analysis parameters
% - Results: Structure with test results

% Set defaults if not specified
if nargin < 2 || isempty(nPerm)
    nPerm = 5000;
end

if nargin < 3 || isempty(H)
    H = 2;
end

if nargin < 4 || isempty(E)
    E = 0.5;
end

% Store the matrix in the cell array
Data{1} = data_matrix;

% Get data dimensions
[nA, nTimepoints] = size(Data{1});

% Initialize variables
maxTFCE = zeros(nPerm,1);

% Calculate observed T-values
display('Calculating Actual Differences...')

% One-sample/dependent t-test
T_Obs = mean(Data{1},1)./(std(Data{1})./sqrt(nA));

% Check for non-zero T-values
if max(abs(T_Obs(:))) < 0.00001
    error('T-values were all 0')
end

% Set up thresholds for TFCE integration
max_t = max(abs(T_Obs));
start_threshold = 0;
threshold_step = 0.1;
thresholds = start_threshold:threshold_step:max_t;

% Calculate TFCE values for observed data
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

display('Done')

% Permutation testing
display('Calculating Permutations...')

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
    
    % Sign-flipping for one-sample/dependent test
    Signs = [-1,1];
    SignSwitch = randsample(Signs, nA, 'true')';
    SignSwitch = repmat(SignSwitch, [1, nTimepoints]);
    
    nData = SignSwitch .* Data{1};
    T_Perm = mean(nData, 1) ./ (std(nData) ./ sqrt(size(nData, 1)));
    
    % Calculate TFCE values for permuted data
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
    
    maxTFCE(i) = max(abs(TFCE_Perm(:)));
    
    % Update progress bar
    waitbar(i / nPerm, hWaitbar, sprintf('Permutations: %d/%d', i, nPerm));
end

% Close waitbar
delete(hWaitbar);

display('Done');

% Calculate p-values and effect sizes
display('Calculating P-Values and Effect Sizes...')

edges = [maxTFCE; max(abs(TFCE_Obs(:)))];
[~,bin] = histc(abs(TFCE_Obs),sort(edges));
P_Values = 1-bin./(nPerm+2);

% Calculate Cohen's d for significant clusters
sig_clusters = P_Values < 0.05;
cohens_d = zeros(size(T_Obs));
processed_timepoints = false(size(T_Obs));

for t = find(sig_clusters)
    if ~processed_timepoints(t)
        cluster_mask = find_temporal_cluster(T_Obs, t);
        processed_timepoints(cluster_mask) = true;
        cluster_data = mean(Data{1}(:,cluster_mask), 2);
        d_value = mean(cluster_data) / std(cluster_data);
        cohens_d(cluster_mask) = d_value;
    end
end

% Save test information
c = clock;
Info.Comments = ['TFCE analysis conducted at ', num2str(c(4)), ':', num2str(c(5)), ' on ', date];

Info.Parameters.nPerm = nPerm;
Info.Parameters.nTimepoints = nTimepoints;
Info.Parameters.GroupSize = nA;
Info.Parameters.H = H;
Info.Parameters.E = E;
Info.Parameters.threshold_step = threshold_step;
Info.Parameters.start_threshold = start_threshold;

Results.Obs = T_Obs;
Results.TFCE_Obs = TFCE_Obs;
Results.maxTFCE = sort(maxTFCE);
Results.P_Values = P_Values;
Results.Cohens_d = cohens_d;

end
