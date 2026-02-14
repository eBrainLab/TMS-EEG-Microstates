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

% Create heatmap
if exist('heatmap', 'file')
    h1 = heatmap(plot_ax, plot_ax, rtf_results.t_post_tms);
    h1.Title = sprintf("Post-TMS: %d ms to %d ms vs. Baseline", ...
        time_interval(1), time_interval(end));
    h1.FontName = 'Times New Roman';
    colormap jet;
else
    % Fallback for MATLAB versions prior to R2017a
    fprintf('Note: heatmap function not available (requires R2017a+). Using imagesc.\n');
    imagesc(rtf_results.t_post_tms);
    set(gca, 'XTick', 1:length(microstates), 'XTickLabel', cellstr(microstates));
    set(gca, 'YTick', 1:length(microstates), 'YTickLabel', cellstr(microstates));
    xlabel('To State', 'FontName', 'Times New Roman');
    ylabel('From State', 'FontName', 'Times New Roman');
    title(sprintf("Post-TMS: %d ms to %d ms vs. Baseline", ...
        time_interval(1), time_interval(end)), 'FontName', 'Times New Roman');
    colorbar;
    colormap jet;
    axis equal tight;
    % Add text annotations for cell values
    data_mat = rtf_results.t_post_tms;
    for r = 1:size(data_mat, 1)
        for c = 1:size(data_mat, 2)
            if ~isnan(data_mat(r, c))
                text(c, r, sprintf('%.2f', data_mat(r, c)), ...
                    'HorizontalAlignment', 'center', 'FontName', 'Times New Roman');
            end
        end
    end
end

% Main figure title (sgtitle requires R2018b+)
if exist('sgtitle', 'file')
    sgtitle(figure_title, 'FontName', 'Times New Roman');
end

%% Save figure
if save_figure
    safe_figure_title = matlab.lang.makeValidName(figure_title);
    if exist('exportgraphics', 'file')
        exportgraphics(gcf, [safe_figure_title, '.png'], 'Resolution', 500);
    else
        % Fallback for MATLAB versions prior to R2020a
        fprintf('Note: exportgraphics not available (requires R2020a+). Using print instead.\n');
        print(gcf, safe_figure_title, '-dpng', '-r500');
    end
    close all;
end

end