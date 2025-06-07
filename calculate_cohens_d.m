function cohen_d_matrix = calculate_cohens_d(combined_data)
% -------------------------------------------------------------------------
% CALCULATE_COHENS_D - Calculate Cohen's d effect size for transition matrices
%
% This function computes Cohen's d for each transition in a microstate
% transition matrix, measuring the standardized effect size.
%
% Syntax:
%   cohen_d_matrix = calculate_cohens_d(combined_data)
%
% Inputs:
%   combined_data - 3D matrix of transition data [states x states x subjects]
%                   Contains transition frequencies for multiple subjects
%
% Outputs:
%   cohen_d_matrix - 2D matrix [states x states] containing Cohen's d values
%                    Diagonal elements (self-transitions) are set to NaN
%
% Formula:
%   Cohen's d = Mean / Standard Deviation
%   This is appropriate for one-sample tests (testing against zero)
%
% Example:
%   transition_data = randn(5, 5, 20); % 5 states, 20 subjects
%   cohen_d = calculate_cohens_d(transition_data);
%
% See also: TEST_RTF, TEST_RTF_WITHIN_TWO_CONDITIONS
% -------------------------------------------------------------------------

%% Initialize
% Get dimensions
transition_size = size(combined_data, 1);

% Preallocate output matrix
cohen_d_matrix = zeros(transition_size, transition_size);

%% Calculate Cohen's d for each transition
for row = 1:transition_size
    for col = 1:transition_size
        if row ~= col  % Exclude self-transitions
            % Extract data for this specific transition across all subjects
            combined_series = squeeze(combined_data(row, col, :));  % Nx1 vector
            
            % Calculate mean and standard deviation
            M = mean(combined_series);
            SD = std(combined_series);
            
            % Calculate Cohen's d
            if SD == 0
                % Handle division by zero
                cohen_d = 0;
            else
                cohen_d = M / SD;
            end
            
            % Store in output matrix
            cohen_d_matrix(row, col) = cohen_d;
        else
            % Set diagonal elements (self-transitions) to NaN
            cohen_d_matrix(row, col) = NaN;
        end
    end
end

end