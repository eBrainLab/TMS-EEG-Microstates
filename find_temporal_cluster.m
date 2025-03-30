function cluster_mask = find_temporal_cluster_for_TFCE(T_values, time_idx)
% -------------------------------------------------------------------------
% CLUSTER FINDING for TFCE integration
%
% This subfunction simply finds contiguous time points around "time_idx"
% that share the same sign in T_values. It does NOT check significance,
% because TFCE is performed by integrating across thresholds, not p-values.
% -------------------------------------------------------------------------
    cluster_mask = false(size(T_values));
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