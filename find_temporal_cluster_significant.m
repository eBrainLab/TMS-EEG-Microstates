function cluster_mask = find_temporal_cluster_significant(T_values, sig_mask, start_idx)
% -------------------------------------------------------------------------
% CLUSTER FINDING for POST-TFCE SIGNIFICANCE
%
% This subfunction finds contiguous time points around "start_idx" that:
%   (1) pass sig_mask (i.e., p<0.05),
%   (2) optionally share the same sign in T_values (if you want that).
%
% This ensures that we only group truly significant time points into a single
% cluster for final effect-size calculations.
% -------------------------------------------------------------------------
    cluster_mask = false(size(T_values));
    cluster_mask(start_idx) = true;

    % (Optional) If you DO want to restrict clusters to points of the same sign:
    ref_sign = sign(T_values(start_idx));

    % Expand forward
    t = start_idx;
    while (t < length(T_values)) ...
          && sig_mask(t+1) ...
          && sign(T_values(t+1)) == ref_sign
        cluster_mask(t+1) = true;
        t = t + 1;
    end

    % Expand backward
    t = start_idx;
    while (t > 1) ...
          && sig_mask(t-1) ...
          && sign(T_values(t-1)) == ref_sign
        cluster_mask(t-1) = true;
        t = t - 1;
    end
end