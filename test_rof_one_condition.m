function [TFCE_Obs, TFCE_Perm, P_Values, Info, Results] = test_rof_one_condition(data_matrix, nPerm, H, E)
% Computes TFCE statistics for one-condition time series data using sign-flipping permutations.
%
% Inputs:
% - data_matrix: Time series data matrix [subjects x timepoints]
% - nPerm: Number of permutations (default 5000)
% - H: TFCE height parameter (default 2)
% - E: TFCE extent parameter (default 0.5)
%
% Outputs:
% - TFCE_Obs: Observed TFCE values
% - TFCE_Perm: Permuted TFCE values (max across timepoints for each permutation)
% - P_Values: P-values from permutation test
% - Info: Structure with analysis parameters
% - Results: Structure with test results (includes T-values, TFCE, effect sizes, etc.)

    % ---------------------------------------------------------------------
    % 1) Set defaults if not specified
    % ---------------------------------------------------------------------
    if nargin < 2 || isempty(nPerm)
        nPerm = 5000;
    end
    if nargin < 3 || isempty(H)
        H = 2;
    end
    if nargin < 4 || isempty(E)
        E = 0.5;
    end

    % ---------------------------------------------------------------------
    % 2) Prepare data, compute T-values
    % ---------------------------------------------------------------------
    Data{1} = data_matrix;  % keep consistent with your original usage
    [nA, nTimepoints] = size(Data{1});

    % Pre-allocate array to store max TFCE values across permutations
    maxTFCE = zeros(nPerm,1);

    % One-sample T-values
    fprintf('Calculating Actual Differences...\n');
    T_Obs = mean(Data{1}, 1) ./ (std(Data{1}) ./ sqrt(nA));

    if max(abs(T_Obs(:))) < 1e-5
        error('T-values were all 0!');
    end

    % ---------------------------------------------------------------------
    % 3) TFCE integration on observed T-values
    % ---------------------------------------------------------------------
    max_t = max(abs(T_Obs));
    start_threshold = 0;
    threshold_step = 0.1;
    thresholds = start_threshold : threshold_step : max_t;

    TFCE_Obs = zeros(size(T_Obs));
    for thresh = thresholds
        above_thresh = abs(T_Obs) > thresh; 
        if any(above_thresh)
            for t = 1 : nTimepoints
                if above_thresh(t)
                    % Use the old sign-based function
                    cluster_mask = find_temporal_cluster(T_Obs, t);
                    cluster_size = sum(cluster_mask);
                    TFCE_Obs(t) = TFCE_Obs(t) + (thresh^H * cluster_size^E * threshold_step);
                end
            end
        end
    end
    fprintf('Done.\n');

    % ---------------------------------------------------------------------
    % 4) Permutation testing (sign flipping) for TFCE
    % ---------------------------------------------------------------------
    fprintf('Calculating Permutations...\n');
    hWaitbar = waitbar(0, 'Calculating Permutations...', ...
        'Name', 'TFCE Permutations', ...
        'CreateCancelBtn', 'setappdata(gcbf, ''canceling'', true)');
    setappdata(hWaitbar, 'canceling', false);

    for i = 1 : nPerm
        % Check for cancel button press
        if getappdata(hWaitbar, 'canceling')
            delete(hWaitbar);
            error('Permutation test canceled by user.');
        end

        % Randomly sign-flip for one-sample T
        Signs = [-1,1];
        SignSwitch = randsample(Signs, nA, 'true')';
        SignSwitch = repmat(SignSwitch, [1, nTimepoints]);

        nData = SignSwitch .* Data{1};
        T_Perm = mean(nData, 1) ./ (std(nData) ./ sqrt(nA));

        % Compute TFCE for permuted data
        temp_TFCE = zeros(size(T_Perm));
        for thresh = thresholds
            above_thresh = abs(T_Perm) > thresh;
            if any(above_thresh)
                for t = 1 : nTimepoints
                    if above_thresh(t)
                        % Use the old sign-based function
                        cluster_mask = find_temporal_cluster(T_Perm, t);
                        cluster_size = sum(cluster_mask);
                        temp_TFCE(t) = temp_TFCE(t) + (thresh^H * cluster_size^E * threshold_step);
                    end
                end
            end
        end

        % Keep track of the maximum TFCE value
        maxTFCE(i) = max(abs(temp_TFCE));
        waitbar(i / nPerm, hWaitbar, sprintf('Permutations: %d/%d', i, nPerm));
    end

    delete(hWaitbar);
    fprintf('Done.\n');

    % ---------------------------------------------------------------------
    % 5) Compute p-values from permutation distribution
    % ---------------------------------------------------------------------
    fprintf('Calculating P-Values and Effect Sizes...\n');
    % Use histogram/bin method for final p-values
    edges = [maxTFCE; max(abs(TFCE_Obs(:)))];
    [~, bin] = histc(abs(TFCE_Obs), sort(edges));
    P_Values = 1 - bin ./ (nPerm + 2);

    % ---------------------------------------------------------------------
    % 6) Identify significant clusters post-TFCE & compute effect sizes
    % ---------------------------------------------------------------------
    % sig_mask is true only where p-values are < 0.05
    sig_mask = (P_Values < 0.05);

    cohens_d = zeros(size(T_Obs));      % store cluster-level d in each timepoint
    processed_timepoints = false(size(T_Obs));
    cluster_d_values = [];              % store the d for each distinct cluster

    for t = 1 : nTimepoints
        % Only check timepoints that are (a) significant and (b) not already processed
        if sig_mask(t) && ~processed_timepoints(t)
            % Identify the contiguous cluster around t, using TFCE_Obs
            % and restricting to the significance mask
            this_cluster = find_temporal_cluster_significant(TFCE_Obs, sig_mask, t);

            % Mark them so we don't re-calculate
            processed_timepoints(this_cluster) = true;

            % Extract data for these timepoints; compute cluster-level effect size
            cluster_data = mean(Data{1}(:, this_cluster), 2);  
            d_value = mean(cluster_data) / std(cluster_data);

            % Assign the same cluster-level Cohen's d to each point in that cluster
            cohens_d(this_cluster) = d_value;

            % Keep track if you want to average across clusters
            cluster_d_values(end+1) = d_value;
        end
    end

    avg_cohens_d = mean(cluster_d_values);
    TFCE_Perm = sort(maxTFCE);

    % ---------------------------------------------------------------------
    % 7) Build output structures
    % ---------------------------------------------------------------------
    c = clock;
    Info.Comments = sprintf('TFCE analysis conducted at %02d:%02d on %s', c(4), c(5), date);
    Info.Parameters.nPerm = nPerm;
    Info.Parameters.nTimepoints = nTimepoints;
    Info.Parameters.GroupSize = nA;
    Info.Parameters.H = H;
    Info.Parameters.E = E;
    Info.Parameters.threshold_step = threshold_step;
    Info.Parameters.start_threshold = start_threshold;

    Results.Obs          = T_Obs;
    Results.TFCE_Obs     = TFCE_Obs;
    Results.TFCE_Perm      = TFCE_Perm;
    Results.P_Values     = P_Values;
    Results.Cohens_d     = cohens_d;
    Results.Avg_Cohens_d = avg_cohens_d;

end
