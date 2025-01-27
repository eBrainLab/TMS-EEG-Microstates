function y = centered_log_ratio(x, delta)
% Applies centered log-ratio (CLR) transformation to the input vector `x`.
% Efficiently handles zeros using the multiplicative replacement method.
%
% Parameters:
%   x     - A compositional vector (1 x D) containing the proportions of
%           microstate occurrences.
%   delta - A scaling factor that determines the proportion of the smallest
%           non-zero component used for replacing zeros.
%
% Returns:
%   y     - The CLR-transformed vector.

% Number of components
D = length(x);

% Identify zero elements
zero_idx = (x == 0);
num_zeros = sum(zero_idx);

if num_zeros > 0
    % Total sum of the non-zero components
    x_nonzero_sum = sum(x(~zero_idx));
    
    % Replace zeros using the multiplicative replacement method
    x_replaced = x;
    epsilon = delta * min(x(x > 0)) / D;
    x_replaced(zero_idx) = epsilon;
    
    % Adjust the non-zero components to compensate for the added mass
    x_replaced(~zero_idx) = x(~zero_idx) - (epsilon * num_zeros) * (x(~zero_idx) / x_nonzero_sum);
    
    % Ensure the data is still compositional (sums to the original total)
    x_replaced = x_replaced / sum(x_replaced) * sum(x);
else
    x_replaced = x;
end

% Calculate the geometric mean of the replaced input vector
geo_mean = geomean(x_replaced);

% Perform the CLR transformation
y = log(x_replaced / geo_mean);
end
