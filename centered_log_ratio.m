function y = centered_log_ratio(x, delta)
% -------------------------------------------------------------------------
% CENTERED_LOG_RATIO - Apply centered log-ratio transformation
%
% This function applies the centered log-ratio (CLR) transformation to a
% compositional vector, efficiently handling zeros using the multiplicative
% replacement method.
%
% Syntax:
%   y = centered_log_ratio(x, delta)
%
% Inputs:
%   x     - Compositional vector (1 x D) containing proportions of
%           microstate occurrences
%   delta - Scaling factor that determines the proportion of the smallest
%           non-zero component used for replacing zeros
%
% Outputs:
%   y     - CLR-transformed vector (1 x D)
%
% Method:
%   1. Replace zeros using multiplicative replacement
%   2. Calculate geometric mean of replaced values
%   3. Apply CLR transformation: log(x / geometric_mean)
%
% Example:
%   proportions = [0.3, 0, 0.5, 0.2];
%   transformed = centered_log_ratio(proportions, 0.5);
%
% See also: EXTRACT_MICROSTATES_FREQUENCIES
% -------------------------------------------------------------------------

%% Initialize
% Get number of components
D = length(x);

% Identify zero elements
zero_idx = (x == 0);
num_zeros = sum(zero_idx);

%% Handle zeros using multiplicative replacement
if num_zeros > 0
    % Total sum of non-zero components
    x_nonzero_sum = sum(x(~zero_idx));
    
    % Initialize replacement vector
    x_replaced = x;
    
    % Calculate replacement value
    epsilon = delta * min(x(x > 0)) / D;
    
    % Replace zeros with epsilon
    x_replaced(zero_idx) = epsilon;
    
    % Adjust non-zero components to maintain compositional constraint
    x_replaced(~zero_idx) = x(~zero_idx) - ...
        (epsilon * num_zeros) * (x(~zero_idx) / x_nonzero_sum);
    
    % Ensure data remains compositional (sums to original total)
    x_replaced = x_replaced / sum(x_replaced) * sum(x);
else
    % No zeros to replace
    x_replaced = x;
end

%% Apply CLR transformation
% Calculate geometric mean of replaced vector
geo_mean = geomean(x_replaced);

% Perform CLR transformation
y = log(x_replaced / geo_mean);

end