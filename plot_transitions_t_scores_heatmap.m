function plot_transitions_t_scores_heatmap(microstates, rtf_results, early_interval, late_interval, figure_title, save_figure)
    % Visualizes the T-scores of Markov transition matrices
    % for early and late responses as heatmaps.
    %
    % Inputs:
    %   - microstates: List of microstate labels for the axes of the heatmap.
    %   - rtf_results: Struct containing T-scores and hypothesis test results.
    %   - early_interval: Time interval (in ms) for early response.
    %   - late_interval: Time interval (in ms) for late response.
    %   - figure_title: Title for the figure.
    %   - save_figure: Boolean to determine whether to save the figure.

    % Set insignificant T-scores to zero
    rtf_results.t_early(rtf_results.h_early == 0) = 0;
    rtf_results.t_late(rtf_results.h_late == 0) = 0;

    % Set diagonal elements to NaN to avoid self-transition bias
    for ii = 1:size(rtf_results.t_early, 1)
        rtf_results.t_early(ii, ii) = NaN;
        rtf_results.t_late(ii, ii) = NaN;
    end

    % Create categorical labels for the heatmap axes
    plot_ax = categorical(microstates);

    % Create a new figure with a title spanning both subplots
    figure;
    set(gcf, 'Position', [100, 200, 1600, 1000]); % [left, bottom, width, height]
    sgtitle(figure_title, 'FontName', 'Times New Roman');

    % Plot the heatmap for early response transitions
    subplot(1, 2, 1);
    h1 = heatmap(plot_ax, plot_ax, rtf_results.t_early);
    h1.Title = sprintf("Early response: %d ms to %d ms vs. pre", early_interval(1), early_interval(end));
    h1.FontName = 'Times New Roman';

    % Plot the heatmap for late response transitions
    subplot(1, 2, 2);
    h2 = heatmap(plot_ax, plot_ax, rtf_results.t_late);
    h2.Title = sprintf("Late response: %d ms to %d ms vs. pre", late_interval(1), late_interval(end));
    h2.FontName = 'Times New Roman';

    % Apply the 'jet' colormap to both heatmaps
    colormap jet;

    % Save the figure as a PNG file with specified resolution if requested
    if save_figure
        safe_figure_title = matlab.lang.makeValidName(figure_title);
        exportgraphics(gcf, [safe_figure_title, '.png'], 'Resolution', 500);
        close all;
    end
end
