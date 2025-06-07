function plot_transitions_t_scores_heatmap(microstates, rtf_results, time_interval, ...
    figure_title, save_figure)
% -------------------------------------------------------------------------
% PLOT_TRANSITIONS_T_SCORES_HEATMAP - Visualize transition T-scores as heatmap
%
% This function creates a heatmap visualization of T-scores from microstate
% transition analysis, showing significant transitions between states.
%
% Syntax:
%   plot_transitions_t_scores_heatmap(microstates, rtf_results, time_interval, ...
%       figure_title, save_figure)
%
% Inputs:
%   microstates   - Cell array of microstate labels for axes
%   rtf_results   - Structure containing T-scores and hypothesis test results
%                   Fields: t_post_tms, h_post_tms
%   time_interval - Time interval vector (in ms) for subtitle
%   figure_title  - String for main figure title
%   save_figure   - Boolean to save figure as PNG
%
% Features:
%   - Non-significant transitions set to zero
%   - Self-transitions (diagonal) set to NaN
%   - Jet colormap for clear visualization
%   - High-resolution export (500 DPI)
%
% Example:
%   plot_transitions_t_scores_heatmap(["A","B","C","D","E"], rtf_results, ...
%       [20,500], 'LPFC Transitions', true);
%
% See also: TEST_RTF, TEST_RTF_WITHIN_TWO_CONDITIONS
% -------------------------------------------------------------------------

%% Prepare data
% Set non-significant T-scores to zero
rtf_results.t_post_tms(rtf_results.h_post_tms == 0) = 0;

% Set diagonal elements (self-transitions) to NaN
for ii = 1:size(rtf_results.t_post_tms, 1)
    rtf_results.t_post_tms(ii, ii) = NaN;
end

% Create categorical labels for heatmap axes
plot_ax = categorical(microstates);

%% Create figure
figure;
set(gcf, 'Position', [100, 200, 1600, 1000]);
title(figure_title, 'FontName', 'Times New Roman');

% Create heatmap
h1 = heatmap(plot_ax, plot_ax, rtf_results.t_post_tms);
h1.Title = sprintf("Post-TMS: %d ms to %d ms vs. Baseline", ...
    time_interval(1), time_interval(end));
h1.FontName = 'Times New Roman';

% Apply jet colormap
colormap jet;

%% Save figure
if save_figure
    safe_figure_title = matlab.lang.makeValidName(figure_title);
    exportgraphics(gcf, [safe_figure_title, '.png'], 'Resolution', 500);
    close all;
end

end