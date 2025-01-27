function cluster_mask = find_temporal_cluster(T_values, timepoint)
    % Helper function to find contiguous temporal clusters
    % Returns logical mask of cluster around given timepoint
    
    cluster_mask = false(size(T_values));
    cluster_mask(timepoint) = true;
    
    % Expand cluster forward in time
    t = timepoint;
    while t < length(T_values) && sign(T_values(t)) == sign(T_values(t+1))
        cluster_mask(t+1) = true;
        t = t + 1;
    end
    
    % Expand cluster backward in time
    t = timepoint;
    while t > 1 && sign(T_values(t)) == sign(T_values(t-1))
        cluster_mask(t-1) = true;
        t = t - 1;
    end
end