function cluster_mask = find_temporal_cluster(T_values, time_idx, sig_mask)
% -------------------------------------------------------------------------
% FIND_TEMPORAL_CLUSTER - Find contiguous temporal clusters
%
% This function identifies contiguous time points around a starting index
% that share the same sign in T-values. Used for both TFCE integration
% and post-TFCE cluster analysis.
%
% Syntax:
%   cluster_mask = find_temporal_cluster(T_values, time_idx)
%   cluster_mask = find_temporal_cluster(T_values, time_idx, sig_mask)
%
% Inputs:
%   T_values - Vector of T-statistics or TFCE values
%   time_idx - Starting time index for cluster search
%   sig_mask - (Optional) Logical mask indicating significant time points
%              If not provided, clusters based only on sign consistency
%
% Outputs:
%   cluster_mask - Logical mask with true at cluster points, false otherwise
%
% Method:
%   1. Start from time_idx
%   2. Expand forward/backward while sign remains consistent
%   3. If sig_mask provided, also require significance
%
% Example:
%   T_vals = [1.2, 1.5, -0.3, -0.8, -0.5, 1.1];
%   cluster = find_temporal_cluster(T_vals, 2);  % Returns [1 1 0 0 0 0]
%
% See also: TEST_ROF_ONE_CONDITION, TEST_ROF_WITHIN_TWO_CONDITIONS
% -------------------------------------------------------------------------

%% Initialize
% Create output mask
cluster_mask = false(size(T_values));
cluster_mask(time_idx) = true;

% Get reference sign from starting point
ref_sign = sign(T_values(time_idx));

% Check if significance constraint is provided
check_significance = nargin >= 3 && ~isempty(sig_mask);

%% Expand cluster forward
t = time_idx;
while t < length(T_values)
    % Check if next point meets criteria
    next_point_valid = sign(T_values(t+1)) == ref_sign;
    
    % Apply significance constraint if provided
    if check_significance
        next_point_valid = next_point_valid && sig_mask(t+1);
    end
    
    % Include in cluster or stop expansion
    if next_point_valid
        cluster_mask(t+1) = true;
        t = t + 1;
    else
        break;
    end
end

%% Expand cluster backward
t = time_idx;
while t > 1
    % Check if previous point meets criteria
    prev_point_valid = sign(T_values(t-1)) == ref_sign;
    
    % Apply significance constraint if provided
    if check_significance
        prev_point_valid = prev_point_valid && sig_mask(t-1);
    end
    
    % Include in cluster or stop expansion
    if prev_point_valid
        cluster_mask(t-1) = true;
        t = t - 1;
    else
        break;
    end
end

end