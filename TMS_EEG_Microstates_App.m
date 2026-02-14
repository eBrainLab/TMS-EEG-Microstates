classdef TMS_EEG_Microstates_App < matlab.apps.AppBase
% TMS_EEG_MICROSTATES_APP  Interactive GUI for TMS-EEG microstate analysis.
%
%   TMS_EEG_Microstates_App   launches the application.
%
%   The app provides four tabs:
%     1. Data & Settings   - load .mat data and configure parameters
%     2. Visualization     - browse microstate sequences per trial
%     3. Analysis          - configure and run ROF / RTF analysis
%     4. Results           - view figures and export
%
%   See also: extract_microstates_frequencies, test_rof_one_condition,
%             test_rof_within_two_conditions, test_rtf,
%             test_rtf_within_two_conditions

    %% ===== UI Component Properties ======================================
    properties (Access = public)
        UIFigure                    matlab.ui.Figure
        TabGroup                    matlab.ui.container.TabGroup

        % --- Tab 1: Data & Settings ---
        DataTab                     matlab.ui.container.Tab
        LoadDataButton              matlab.ui.control.Button
        FileNameLabel               matlab.ui.control.Label
        InfoPanel                   matlab.ui.container.Panel
        NumSubjectsLabel            matlab.ui.control.Label
        StimSitesLabel              matlab.ui.control.Label
        MicrostatesDetectedLabel    matlab.ui.control.Label
        NumTrialsLabel              matlab.ui.control.Label
        SettingsPanel               matlab.ui.container.Panel
        TMSIndexSpinner             matlab.ui.control.Spinner
        SampleRateSpinner           matlab.ui.control.Spinner
        MicrostateLabelsField       matlab.ui.control.EditField

        % --- Tab 2: Visualization ---
        VizTab                      matlab.ui.container.Tab
        SubjectDropDown             matlab.ui.control.DropDown
        SiteDropDown                matlab.ui.control.DropDown
        TrialDropDown               matlab.ui.control.DropDown
        RibbonAxes                  matlab.ui.control.UIAxes
        SequenceAxes                matlab.ui.control.UIAxes

        % --- Tab 3: Analysis ---
        AnalysisTab                 matlab.ui.container.Tab
        AnalysisTypeGroup           matlab.ui.container.ButtonGroup
        SingleSiteRadio             matlab.ui.control.RadioButton
        BetweenSiteRadio            matlab.ui.control.RadioButton
        AnalysisSite1DropDown       matlab.ui.control.DropDown
        AnalysisSite2DropDown       matlab.ui.control.DropDown
        AnalysisSite2Label          matlab.ui.control.Label
        ROFCheckBox                 matlab.ui.control.CheckBox
        RTFCheckBox                 matlab.ui.control.CheckBox
        BaselineStartField          matlab.ui.control.NumericEditField
        BaselineEndField            matlab.ui.control.NumericEditField
        PostStimStartField          matlab.ui.control.NumericEditField
        PostStimEndROFField         matlab.ui.control.NumericEditField
        PostStimEndRTFField         matlab.ui.control.NumericEditField
        % Ignore window is derived: BaselineEnd .. PostStimStart
        NPermSpinner                matlab.ui.control.Spinner
        CorrectionDropDown          matlab.ui.control.DropDown
        RunAnalysisButton           matlab.ui.control.Button
        AnalysisStatusLabel         matlab.ui.control.Label

        % --- Tab 4: Results ---
        ResultsTab                  matlab.ui.container.Tab
        ResultDropDown              matlab.ui.control.DropDown
        ResultAxes                  matlab.ui.control.UIAxes
        ExportFigureButton          matlab.ui.control.Button
        ExportResultsButton         matlab.ui.control.Button
    end

    %% ===== Data Properties ==============================================
    properties (Access = private)
        DataStruct                  % Raw loaded structure array
        DataStructProcessed         % Processed structure (with frequencies)
        Microstates                 % String array of microstate labels
        StimSites                   % String array of stimulation site names
        NumSubjects     double      % Number of subjects
        CommonTimeArray double      % Time vector in ms
        TMSIndex        double      % TMS onset sample index
        SampleRate      double      % Sampling rate in Hz
        AnalysisResults struct      % Array of stored results
        DataLoaded      logical     % True after successful load
        DataFilePath    char        % Path to loaded .mat file
        VizEnabled      logical     % True when Visualization tab is ready
        AnalysisEnabled logical     % True when Analysis tab is ready
        ResultsEnabled  logical     % True when Results tab is ready
    end

    %% ===== Callback Methods =============================================
    methods (Access = private)

        % ----- Tab navigation guard --------------------------------------

        function TabGroupSelectionChanged(app, event)
            newTab = event.NewValue;
            if (newTab == app.VizTab && ~app.VizEnabled) || ...
               (newTab == app.AnalysisTab && ~app.AnalysisEnabled) || ...
               (newTab == app.ResultsTab && ~app.ResultsEnabled)
                app.TabGroup.SelectedTab = event.OldValue;
                uialert(app.UIFigure, ...
                    'Please load data first.', 'Tab Not Available');
            end
        end

        % ----- Tab 1 callbacks -------------------------------------------

        function LoadDataButtonPushed(app, ~)
            [file, path] = uigetfile('*.mat', 'Select Microstate Data File');
            if isequal(file, 0), return; end

            filepath = fullfile(path, file);
            app.DataFilePath = filepath;

            d = uiprogressdlg(app.UIFigure, 'Title', 'Loading', ...
                'Message', 'Loading data file...', 'Indeterminate', 'on');

            try
                loaded = load(filepath);

                % Locate variable that contains microstate_data
                vars = fieldnames(loaded);
                found = false;
                for v = 1:length(vars)
                    if isstruct(loaded.(vars{v})) && ...
                            isfield(loaded.(vars{v}), 'microstate_data')
                        app.DataStruct = loaded.(vars{v});
                        found = true;
                        break;
                    end
                end

                if ~found
                    close(d);
                    uialert(app.UIFigure, ...
                        'No structure with a "microstate_data" field found.', ...
                        'Invalid Data');
                    return;
                end

                app.NumSubjects = length(app.DataStruct);

                % Detect stimulation sites
                ms_data = app.DataStruct(1).microstate_data;
                app.StimSites = string(fieldnames(ms_data))';

                % Detect microstate labels from first trial of first site
                first_site = app.StimSites(1);
                first_trial = ms_data.(first_site){1};
                if iscell(first_trial)
                    unique_labels = unique(string(first_trial));
                else
                    unique_labels = unique(first_trial);
                end
                app.Microstates = sort(unique_labels(:))';

                % Read settings and build time array
                app.TMSIndex   = app.TMSIndexSpinner.Value;
                app.SampleRate = app.SampleRateSpinner.Value;
                nSamples = length(first_trial);
                app.CommonTimeArray = ...
                    ((1:nSamples) - app.TMSIndex) * (1000 / app.SampleRate);

                % Update info labels
                app.FileNameLabel.Text = sprintf('File: %s', file);
                app.NumSubjectsLabel.Text = ...
                    sprintf('Subjects: %d', app.NumSubjects);
                app.StimSitesLabel.Text = ...
                    sprintf('Stimulation Sites: %s', strjoin(upper(app.StimSites), ', '));
                app.MicrostatesDetectedLabel.Text = ...
                    sprintf('Microstates: %s', strjoin(app.Microstates, ', '));
                nTrials = length(ms_data.(first_site));
                app.NumTrialsLabel.Text = ...
                    sprintf('Trials (Subj 1, %s): %d', upper(char(first_site)), nTrials);
                app.MicrostateLabelsField.Value = strjoin(app.Microstates, ', ');

                % Enable tabs and populate dropdowns
                app.VizEnabled      = true;
                app.AnalysisEnabled = true;
                app.populateVizDropdowns();
                app.populateAnalysisDropdowns();
                app.DataLoaded = true;

                close(d);
                app.updateVisualization();

            catch ME
                close(d);
                uialert(app.UIFigure, ...
                    sprintf('Error loading file:\n%s', ME.message), 'Load Error');
            end
        end

        function TMSIndexChanged(app, ~)
            if ~app.DataLoaded, return; end
            app.TMSIndex = app.TMSIndexSpinner.Value;
            app.rebuildTimeArray();
            app.updateVisualization();
        end

        function SampleRateChanged(app, ~)
            if ~app.DataLoaded, return; end
            app.SampleRate = app.SampleRateSpinner.Value;
            app.rebuildTimeArray();
            app.updateVisualization();
        end

        function MicrostateLabelsChanged(app, ~)
            if ~app.DataLoaded, return; end
            raw = app.MicrostateLabelsField.Value;
            app.Microstates = strtrim(string(strsplit(raw, ',')));
            app.populateAnalysisDropdowns();
        end

        % ----- Tab 2 callbacks -------------------------------------------

        function VizSelectionChanged(app, ~)
            if ~app.DataLoaded, return; end
            app.updateTrialDropdown();
            app.updateVisualization();
        end

        function TrialChanged(app, ~)
            if ~app.DataLoaded, return; end
            app.updateVisualization();
        end

        % ----- Tab 3 callbacks -------------------------------------------

        function AnalysisTypeChanged(app, ~)
            isBetween = app.BetweenSiteRadio.Value;
            if isBetween
                app.AnalysisSite2DropDown.Enable = 'on';
                app.AnalysisSite2Label.Enable    = 'on';
            else
                app.AnalysisSite2DropDown.Enable = 'off';
                app.AnalysisSite2Label.Enable    = 'off';
            end
        end

        function RunAnalysisButtonPushed(app, ~)
            app.runAnalysis();
        end

        % ----- Tab 4 callbacks -------------------------------------------

        function ResultDropDownChanged(app, ~)
            app.displaySelectedResult();
        end

        function ExportFigureButtonPushed(app, ~)
            [file, path] = uiputfile( ...
                {'*.png'; '*.pdf'; '*.fig'}, 'Save Figure As');
            if isequal(file, 0), return; end
            filepath = fullfile(path, file);
            try
                if endsWith(filepath, '.fig')
                    fig = figure('Visible', 'off');
                    copyobj(app.ResultAxes, fig);
                    savefig(fig, filepath);
                    close(fig);
                elseif exist('exportgraphics', 'file')
                    exportgraphics(app.ResultAxes, filepath, 'Resolution', 300);
                else
                    fig = figure('Visible', 'off');
                    newAx = copyobj(app.ResultAxes, fig);
                    set(newAx, 'Units', 'normalized', 'Position', [.1 .1 .8 .8]);
                    print(fig, filepath, '-dpng', '-r300');
                    close(fig);
                end
                uialert(app.UIFigure, ...
                    sprintf('Figure saved to:\n%s', filepath), ...
                    'Export Successful', 'Icon', 'success');
            catch ME
                uialert(app.UIFigure, ...
                    sprintf('Export failed: %s', ME.message), 'Export Error');
            end
        end

        function ExportResultsButtonPushed(app, ~)
            [file, path] = uiputfile('*.mat', 'Save Results As');
            if isequal(file, 0), return; end
            filepath = fullfile(path, file);
            try
                AnalysisResults = app.AnalysisResults; %#ok<PROPLC>
                save(filepath, 'AnalysisResults');
                uialert(app.UIFigure, ...
                    sprintf('Results saved to:\n%s', filepath), ...
                    'Export Successful', 'Icon', 'success');
            catch ME
                uialert(app.UIFigure, ...
                    sprintf('Export failed: %s', ME.message), 'Export Error');
            end
        end
    end

    %% ===== Helper Methods ===============================================
    methods (Access = private)

        function updateProgress(~, dlg, fraction, message)
            % Update the uiprogressdlg from a progress callback
            if isvalid(dlg)
                dlg.Value   = min(max(fraction, 0), 1);
                dlg.Message = message;
                drawnow;
            end
        end

        function rebuildTimeArray(app)
            first_site = app.StimSites(1);
            nSamples = length( ...
                app.DataStruct(1).microstate_data.(first_site){1});
            app.CommonTimeArray = ...
                ((1:nSamples) - app.TMSIndex) * (1000 / app.SampleRate);
        end

        % ----- Dropdown helpers ------------------------------------------

        function populateVizDropdowns(app)
            subj_items = arrayfun(@(x) sprintf('Subject %d', x), ...
                1:app.NumSubjects, 'UniformOutput', false);
            app.SubjectDropDown.Items     = subj_items;
            app.SubjectDropDown.ItemsData = 1:app.NumSubjects;
            app.SubjectDropDown.Value     = 1;

            site_items = cellstr(upper(app.StimSites));
            app.SiteDropDown.Items     = site_items;
            app.SiteDropDown.ItemsData = cellstr(app.StimSites);
            app.SiteDropDown.Value     = char(app.StimSites(1));

            app.updateTrialDropdown();
        end

        function updateTrialDropdown(app)
            sub_idx = app.SubjectDropDown.Value;
            site    = string(app.SiteDropDown.Value);

            if isfield(app.DataStruct(sub_idx).microstate_data, site)
                nTrials = length(app.DataStruct(sub_idx).microstate_data.(site));
            else
                nTrials = 0;
            end

            if nTrials > 0
                trial_items = arrayfun(@(x) sprintf('Trial %d', x), ...
                    1:nTrials, 'UniformOutput', false);
                app.TrialDropDown.Items     = trial_items;
                app.TrialDropDown.ItemsData = 1:nTrials;
                app.TrialDropDown.Value     = 1;
            else
                app.TrialDropDown.Items     = {'No trials'};
                app.TrialDropDown.ItemsData = 0;
                app.TrialDropDown.Value     = 0;
            end
        end

        function populateAnalysisDropdowns(app)
            site_items = cellstr(upper(app.StimSites));
            app.AnalysisSite1DropDown.Items     = site_items;
            app.AnalysisSite1DropDown.ItemsData = cellstr(app.StimSites);
            app.AnalysisSite1DropDown.Value     = char(app.StimSites(1));

            app.AnalysisSite2DropDown.Items     = site_items;
            app.AnalysisSite2DropDown.ItemsData = cellstr(app.StimSites);
            if length(app.StimSites) > 1
                app.AnalysisSite2DropDown.Value = char(app.StimSites(2));
            end
        end

        % ----- Visualization (Tab 2) ------------------------------------

        function updateVisualization(app)
            sub_idx   = app.SubjectDropDown.Value;
            site      = string(app.SiteDropDown.Value);
            trial_idx = app.TrialDropDown.Value;
            if isempty(sub_idx) || isequal(trial_idx, 0), return; end

            if ~isfield(app.DataStruct(sub_idx).microstate_data, site)
                return;
            end
            ms_data = app.DataStruct(sub_idx).microstate_data.(site);
            if trial_idx > length(ms_data), return; end

            sequence = ms_data{trial_idx};
            if iscell(sequence), sequence = string(sequence); end

            nSamples  = length(sequence);
            nStates   = length(app.Microstates);
            colors    = app.getStateColors(nStates);
            time_ms   = app.CommonTimeArray;
            if length(time_ms) ~= nSamples
                time_ms = ((1:nSamples) - app.TMSIndex) * (1000/app.SampleRate);
            end

            % Convert labels to numeric codes
            numeric_seq = zeros(1, nSamples);
            for s = 1:nStates
                numeric_seq(sequence == app.Microstates(s)) = s;
            end

            % Common x-axis range
            x_range = [time_ms(1), time_ms(end)];

            % ---- Ribbon plot (top) ----
            delete(findall(app.RibbonAxes, 'Type', 'Legend'));
            cla(app.RibbonAxes, 'reset');
            imagesc(app.RibbonAxes, time_ms, 1, numeric_seq);
            colormap(app.RibbonAxes, colors);
            app.RibbonAxes.CLim = [0.5, nStates + 0.5];
            app.RibbonAxes.YTick = [];
            app.RibbonAxes.XLim  = x_range;
            app.RibbonAxes.FontSize = 10;
            hold(app.RibbonAxes, 'on');
            plot(app.RibbonAxes, [0 0], [0.5 1.5], 'w--', 'LineWidth', 2);
            hold(app.RibbonAxes, 'off');
            title(app.RibbonAxes, sprintf('Subject %d  |  %s  |  Trial %d', ...
                sub_idx, upper(char(site)), trial_idx), 'FontSize', 13);

            % ---- Step / segment plot (bottom) ----
            delete(findall(app.SequenceAxes, 'Type', 'Legend'));
            cla(app.SequenceAxes, 'reset');
            hold(app.SequenceAxes, 'on');

            % Find contiguous segments for efficiency
            changes = [1, find(diff(numeric_seq) ~= 0) + 1, nSamples + 1];
            legend_handles = gobjects(nStates, 1);
            legend_drawn   = false(nStates, 1);

            for seg = 1:length(changes) - 1
                seg_start = changes(seg);
                seg_end   = changes(seg + 1) - 1;
                state     = numeric_seq(seg_start);
                if state < 1 || state > nStates, continue; end

                t1 = time_ms(seg_start);
                if seg_end < nSamples
                    t2 = time_ms(seg_end + 1);
                else
                    t2 = time_ms(seg_end) + (1000 / app.SampleRate);
                end

                h = patch(app.SequenceAxes, ...
                    [t1, t2, t2, t1], ...
                    [state-0.4, state-0.4, state+0.4, state+0.4], ...
                    colors(state, :), 'EdgeColor', 'none', 'FaceAlpha', 0.8);

                if ~legend_drawn(state)
                    legend_handles(state) = h;
                    legend_drawn(state)   = true;
                else
                    set(h, 'HandleVisibility', 'off');
                end
            end

            % TMS marker
            plot(app.SequenceAxes, [0 0], [0.2, nStates + 0.8], ...
                'r--', 'LineWidth', 2, 'HandleVisibility', 'off');

            % Format axes - match x-range to ribbon, A on top
            app.SequenceAxes.XLim       = x_range;
            app.SequenceAxes.YTick      = 1:nStates;
            app.SequenceAxes.YTickLabel = cellstr(app.Microstates);
            app.SequenceAxes.YLim       = [0.2, nStates + 0.8];
            app.SequenceAxes.YDir       = 'reverse';
            app.SequenceAxes.FontSize   = 11;
            xlabel(app.SequenceAxes, 'Time (ms)', 'FontSize', 12);
            ylabel(app.SequenceAxes, 'Microstate',  'FontSize', 12);
            grid(app.SequenceAxes, 'on');

            % Link x-axes so zoom/pan stays synchronized
            linkaxes([app.RibbonAxes, app.SequenceAxes], 'x');

            % Legend (only drawn states)
            valid = legend_drawn;
            if any(valid)
                legend(app.SequenceAxes, legend_handles(valid), ...
                    cellstr(app.Microstates(valid)), ...
                    'Location', 'northeast', 'FontSize', 10);
            end

            hold(app.SequenceAxes, 'off');
        end

        % ----- Analysis (Tab 3) -----------------------------------------

        function runAnalysis(app)
            if ~app.DataLoaded
                uialert(app.UIFigure, 'Please load data first.', 'No Data');
                return;
            end

            microstates = app.Microstates;
            do_rof = app.ROFCheckBox.Value;
            do_rtf = app.RTFCheckBox.Value;
            if ~do_rof && ~do_rtf
                uialert(app.UIFigure, ...
                    'Select at least one feature (ROF or RTF).', 'No Features');
                return;
            end

            % Read parameters
            baseline_start  = app.BaselineStartField.Value;
            baseline_end    = app.BaselineEndField.Value;
            post_stim_start = app.PostStimStartField.Value;
            post_stim_end_rof = app.PostStimEndROFField.Value;
            post_stim_end_rtf = app.PostStimEndRTFField.Value;
            nshuffles = app.NPermSpinner.Value;

            switch app.CorrectionDropDown.Value
                case 'Bonferroni', apply_correction = 'bonferroni';
                case 'FDR',        apply_correction = 'fdr';
                otherwise,         apply_correction = '';
            end

            is_single  = app.SingleSiteRadio.Value;
            site1      = string(app.AnalysisSite1DropDown.Value);
            site2      = string(app.AnalysisSite2DropDown.Value);

            % Compute indices from time array
            ta = app.CommonTimeArray;
            baseline_idx  = find(ta >= baseline_start & ta <= baseline_end);
            test_idx_rof  = find(ta >= post_stim_start & ta <= post_stim_end_rof);

            time_window_ranges = struct( ...
                'baseline', [baseline_start, baseline_end], ...
                'post_tms', [post_stim_start, post_stim_end_rtf]);

            features2extract = {};
            if do_rof, features2extract{end+1} = 'ROF'; end
            if do_rtf, features2extract{end+1} = 'RTF'; end

            % Reset results
            app.AnalysisResults = struct( ...
                'name', {}, 'type', {}, 'data', {}, ...
                'site1', {}, 'site2', {}, ...
                'p_values_matrix', {}, 'data_struct_used', {});

            d = uiprogressdlg(app.UIFigure, ...
                'Title', 'Analysis in Progress', ...
                'Message', 'Extracting microstate frequencies...', ...
                'Cancelable', 'on');

            try
                % --- Step 1: extract frequencies -------------------------
                d.Value   = 0.05;
                d.Message = 'Extracting microstate frequencies...';
                app.AnalysisStatusLabel.Text = 'Status: Extracting frequencies...';
                drawnow;

                evalc_output = evalc('data_proc = extract_microstates_frequencies(app.DataStruct, microstates, baseline_idx, time_window_ranges, features2extract);'); %#ok<NASGU>
                app.DataStructProcessed = data_proc;

                if d.CancelRequested
                    app.AnalysisStatusLabel.Text = 'Status: Cancelled.';
                    close(d); return;
                end

                % --- Step 2: run tests -----------------------------------
                if is_single
                    app.runSingleSiteAnalysis(d, data_proc, ...
                        char(site1), microstates, test_idx_rof, ...
                        nshuffles, apply_correction, do_rof, do_rtf, ta);
                else
                    app.runBetweenSiteAnalysis(d, data_proc, ...
                        char(site1), char(site2), microstates, ...
                        test_idx_rof, nshuffles, apply_correction, ...
                        do_rof, do_rtf, ta);
                end

                % --- Step 3: populate results dropdown -------------------
                disp_items   = {};
                disp_indices = [];
                for r = 1:length(app.AnalysisResults)
                    tp = app.AnalysisResults(r).type;
                    if strcmp(tp, 'ROF') || strcmp(tp, 'RTF')
                        disp_items{end+1}   = app.AnalysisResults(r).name; %#ok<AGROW>
                        disp_indices(end+1) = r; %#ok<AGROW>
                    end
                end

                app.ResultDropDown.Items     = disp_items;
                app.ResultDropDown.ItemsData = disp_indices;

                if ~isempty(disp_items)
                    app.ResultDropDown.Value = disp_indices(1);
                    app.ResultsEnabled = true;
                    app.displaySelectedResult();
                end

                d.Value   = 1;
                d.Message = 'Analysis complete!';
                app.AnalysisStatusLabel.Text = 'Status: Analysis complete.';
                pause(0.5);
                close(d);

                app.TabGroup.SelectedTab = app.ResultsTab;

            catch ME
                if isvalid(d), close(d); end
                app.AnalysisStatusLabel.Text = ...
                    sprintf('Status: Error - %s', ME.message);
                uialert(app.UIFigure, ...
                    sprintf('Analysis error:\n%s\n\nIn: %s (line %d)', ...
                    ME.message, ME.stack(1).name, ME.stack(1).line), ...
                    'Analysis Error');
            end
        end

        function runSingleSiteAnalysis(app, d, data_proc, current_site, ...
                microstates, test_idx_rof, nshuffles, apply_correction, ...
                do_rof, do_rtf, time_array)

            total_steps = 0;
            if do_rof, total_steps = total_steps + length(microstates); end
            if do_rtf, total_steps = total_steps + 1; end
            step = 0;
            num_subjects = length(data_proc);

            if do_rof
                p_values_matrix = NaN(length(microstates), length(time_array));

                for ss = 1:length(microstates)
                    if d.CancelRequested
                        close(d);
                        app.AnalysisStatusLabel.Text = 'Status: Cancelled.';
                        return;
                    end
                    step = step + 1;
                    current_state = microstates(ss);
                    d.Value   = 0.1 + 0.8 * step / total_steps;
                    d.Message = sprintf('ROF: %s - Microstate %s (%d/%d)', ...
                        upper(current_site), current_state, ss, length(microstates));
                    app.AnalysisStatusLabel.Text = ['Status: ', d.Message];
                    drawnow;

                    data_mat = zeros(num_subjects, length(test_idx_rof));
                    for i = 1:num_subjects
                        if isfield(data_proc(i), current_site) && ...
                                ~isempty(data_proc(i).(current_site))
                            data_mat(i, :) = data_proc(i).(current_site) ...
                                .occurrences_clr_bc.(current_state)(test_idx_rof)';
                        end
                    end

                    % Progress callback for permutations
                    baseFrac = 0.1 + 0.8 * (step - 1) / total_steps;
                    stepFrac = 0.8 / total_steps;
                    pfcn = @(frac, msg) app.updateProgress(d, ...
                        baseFrac + stepFrac * frac, ...
                        sprintf('ROF %s - %s: %s', upper(current_site), current_state, msg));

                    [~, ~, P_Values, ~, ROF_Result] = ...
                        test_rof_one_condition(data_mat, nshuffles, [], [], pfcn);

                    full_p = NaN(1, length(time_array));
                    full_p(test_idx_rof) = P_Values;
                    p_values_matrix(ss, :) = full_p;

                    app.appendResult(sprintf('ROF - %s - %s', ...
                        upper(current_site), current_state), ...
                        'ROF_individual', ROF_Result, current_site, '', ...
                        [], []);
                end

                app.appendResult(sprintf('ROF Overview - %s', ...
                    upper(current_site)), 'ROF', [], current_site, '', ...
                    p_values_matrix, data_proc);
            end

            if do_rtf
                step = step + 1;
                d.Value   = 0.1 + 0.8 * step / total_steps;
                d.Message = sprintf('RTF: %s', upper(current_site));
                app.AnalysisStatusLabel.Text = ['Status: ', d.Message];
                drawnow;

                evalc_output = evalc('RTF_Result = test_rtf(data_proc, current_site, apply_correction);'); %#ok<NASGU>

                app.appendResult(sprintf('RTF - %s', upper(current_site)), ...
                    'RTF', RTF_Result, current_site, '', [], []);
            end
        end

        function runBetweenSiteAnalysis(app, d, data_proc, site1, site2, ...
                microstates, test_idx_rof, nshuffles, apply_correction, ...
                do_rof, do_rtf, time_array)

            % Filter subjects that have both sites
            num_total  = length(data_proc);
            valid_mask = false(num_total, 1);
            for i = 1:num_total
                if isfield(data_proc(i), site1) && ...
                        ~isempty(data_proc(i).(site1)) && ...
                        isfield(data_proc(i), site2) && ...
                        ~isempty(data_proc(i).(site2))
                    valid_mask(i) = true;
                end
            end
            data_comp = data_proc(valid_mask);

            if isempty(data_comp)
                close(d);
                uialert(app.UIFigure, ...
                    'No subjects with data for both sites.', ...
                    'No Valid Subjects');
                return;
            end

            total_steps = 0;
            if do_rof, total_steps = total_steps + length(microstates); end
            if do_rtf, total_steps = total_steps + 1; end
            step = 0;
            num_comp = length(data_comp);

            if do_rof
                p_values_matrix = NaN(length(microstates), length(time_array));

                for ss = 1:length(microstates)
                    if d.CancelRequested
                        close(d);
                        app.AnalysisStatusLabel.Text = 'Status: Cancelled.';
                        return;
                    end
                    step = step + 1;
                    current_state = microstates(ss);
                    d.Value   = 0.1 + 0.8 * step / total_steps;
                    d.Message = sprintf('ROF: %s vs %s - %s (%d/%d)', ...
                        upper(site1), upper(site2), current_state, ...
                        ss, length(microstates));
                    app.AnalysisStatusLabel.Text = ['Status: ', d.Message];
                    drawnow;

                    dm1 = zeros(num_comp, length(test_idx_rof));
                    dm2 = zeros(num_comp, length(test_idx_rof));
                    for i = 1:num_comp
                        dm1(i, :) = data_comp(i).(site1) ...
                            .occurrences_clr_bc.(current_state)(test_idx_rof)';
                        dm2(i, :) = data_comp(i).(site2) ...
                            .occurrences_clr_bc.(current_state)(test_idx_rof)';
                    end

                    % Progress callback for permutations
                    baseFrac = 0.1 + 0.8 * (step - 1) / total_steps;
                    stepFrac = 0.8 / total_steps;
                    pfcn = @(frac, msg) app.updateProgress(d, ...
                        baseFrac + stepFrac * frac, ...
                        sprintf('ROF %s vs %s - %s: %s', upper(site1), upper(site2), current_state, msg));

                    [~, ~, P_Values, ~, ROF_Result] = ...
                        test_rof_within_two_conditions(dm1, dm2, nshuffles, [], [], pfcn);

                    full_p = NaN(1, length(time_array));
                    full_p(test_idx_rof) = P_Values;
                    p_values_matrix(ss, :) = full_p;

                    app.appendResult(sprintf('ROF - %s vs %s - %s', ...
                        upper(site1), upper(site2), current_state), ...
                        'ROF_individual', ROF_Result, site1, site2, [], []);
                end

                app.appendResult(sprintf('ROF Overview - %s vs %s', ...
                    upper(site1), upper(site2)), 'ROF', [], ...
                    site1, site2, p_values_matrix, data_comp);
            end

            if do_rtf
                step = step + 1;
                d.Value   = 0.1 + 0.8 * step / total_steps;
                d.Message = sprintf('RTF: %s vs %s', upper(site1), upper(site2));
                app.AnalysisStatusLabel.Text = ['Status: ', d.Message];
                drawnow;

                evalc_output = evalc('RTF_Result = test_rtf_within_two_conditions(data_comp, [string(site1), string(site2)], apply_correction);'); %#ok<NASGU>

                app.appendResult(sprintf('RTF - %s vs %s', ...
                    upper(site1), upper(site2)), 'RTF', ...
                    RTF_Result, site1, site2, [], []);
            end
        end

        function appendResult(app, name, type, data, site1, site2, pmat, ds)
            idx = length(app.AnalysisResults) + 1;
            app.AnalysisResults(idx).name             = name;
            app.AnalysisResults(idx).type             = type;
            app.AnalysisResults(idx).data             = data;
            app.AnalysisResults(idx).site1            = site1;
            app.AnalysisResults(idx).site2            = site2;
            app.AnalysisResults(idx).p_values_matrix  = pmat;
            app.AnalysisResults(idx).data_struct_used = ds;
        end

        % ----- Results display (Tab 4) -----------------------------------

        function displaySelectedResult(app)
            if isempty(app.AnalysisResults), return; end
            idx    = app.ResultDropDown.Value;
            result = app.AnalysisResults(idx);
            delete(findall(app.ResultAxes, 'Type', 'Legend'));
            delete(findall(app.ResultAxes, 'Type', 'ColorBar'));
            cla(app.ResultAxes, 'reset');

            if strcmp(result.type, 'ROF')
                app.plotROFOverview(result, app.ResultAxes);
            elseif strcmp(result.type, 'RTF')
                app.plotRTFHeatmap(result, app.ResultAxes);
            end
        end

        % ----- Plotting methods (adapted from standalone functions) ------

        function plotROFOverview(app, result, ax)
            % Adapted from plot_microstate_occurrence.m  (renders into UIAxes)

            ds   = result.data_struct_used;
            pmat = result.p_values_matrix;
            s1   = result.site1;
            s2   = result.site2;
            is_comp = ~isempty(s2);

            ms  = app.Microstates;
            ta  = app.CommonTimeArray;
            nS  = length(ms);
            nT  = length(ta);
            nSub = length(ds);

            % Display range
            fig_range = [max(-200, ta(1)), min(600, ta(end))];
            ti = find(ta >= fig_range(1) & ta <= fig_range(2));

            % Collect per-microstate data across subjects
            all_data = struct();
            for m = 1:nS
                st = char(ms(m));
                all_data.(st) = NaN(nT, nSub);
            end

            for sub = 1:nSub
                for m = 1:nS
                    st = char(ms(m));
                    if is_comp
                        if isfield(ds(sub), s1) && isfield(ds(sub), s2) && ...
                                isfield(ds(sub).(s1), 'occurrences_clr_bc') && ...
                                isfield(ds(sub).(s2), 'occurrences_clr_bc')
                            all_data.(st)(:, sub) = ...
                                ds(sub).(s1).occurrences_clr_bc.(ms(m)) - ...
                                ds(sub).(s2).occurrences_clr_bc.(ms(m));
                        end
                    else
                        if isfield(ds(sub), s1) && ...
                                isfield(ds(sub).(s1), 'occurrences_clr_bc')
                            all_data.(st)(:, sub) = ...
                                ds(sub).(s1).occurrences_clr_bc.(ms(m));
                        end
                    end
                end
            end

            % Statistics
            df = nSub - 1;
            if exist('tinv', 'file')
                tc = tinv(1 - 0.025, df);
            elseif exist('betaincinv', 'file')
                x  = betaincinv(0.05, df/2, 0.5);
                tc = sqrt(df * (1 - x) / x);
            else
                z  = sqrt(2) * erfinv(0.95);
                tc = z + (z^3 + z) / (4 * df);
            end

            avgs = struct();  cis = struct();
            for m = 1:nS
                st  = char(ms(m));
                avg = mean(all_data.(st), 2, 'omitnan');
                se  = std(all_data.(st), 0, 2, 'omitnan') ./ ...
                      sqrt(sum(~isnan(all_data.(st)), 2));
                ci  = tc * se;

                % Smooth
                if exist('smooth', 'file')
                    avg = smooth(avg, 0.02, 'loess');
                    ci  = smooth(ci,  0.02, 'loess');
                else
                    ws = max(3, round(0.02 * length(avg)));
                    if mod(ws,2)==0, ws = ws + 1; end
                    k = ones(ws,1)/ws;
                    avg = conv(avg, k, 'same');
                    ci  = conv(ci,  k, 'same');
                end

                % Exclude ignore window
                excl = ta >= app.BaselineEndField.Value & ...
                       ta <= app.PostStimStartField.Value;
                avg(excl) = NaN;
                ci(excl)  = NaN;

                avgs.(st) = avg;
                cis.(st)  = ci;
            end

            colors = app.getStateColors(nS);
            hold(ax, 'on');

            % Shaded CIs and curves
            ign_s = app.BaselineEndField.Value;
            ign_e = app.PostStimStartField.Value;
            for m = 1:nS
                st  = char(ms(m));
                avg = avgs.(st)(ti)';
                ci  = cis.(st)(ti)';

                for region = 1:2
                    if region == 1
                        ri = ta(ti) < ign_s;
                    else
                        ri = ta(ti) > ign_e;
                    end
                    ok = ri & ~isnan(avg) & ~isnan(ci);
                    if any(ok)
                        xf = [ta(ti(ok)), fliplr(ta(ti(ok)))];
                        yf = [avg(ok)+ci(ok), fliplr(avg(ok)-ci(ok))];
                        fill(ax, xf, yf, colors(m,:), ...
                            'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
                            'HandleVisibility', 'off');
                    end
                end

                plot(ax, ta(ti), avg, 'Color', colors(m,:), ...
                    'LineWidth', 2, 'DisplayName', char(ms(m)));
            end

            % Compute Cohen's d per microstate per timepoint: d = mean/std
            cohens_d_all = struct();
            for m = 1:nS
                st  = char(ms(m));
                d_arr = all_data.(st);          % nT x nSub
                mu  = mean(d_arr, 2, 'omitnan');
                sd  = std(d_arr, 0, 2, 'omitnan');
                sd(sd == 0) = NaN;
                cohens_d_all.(st) = mu ./ sd;   % nT x 1
            end

            % Significance bars
            avg_mat = zeros(nS, nT);
            for m = 1:nS, avg_mat(m,:) = avgs.(char(ms(m)))'; end

            g_min = inf; g_max = -inf;
            for m = 1:nS
                v = avgs.(char(ms(m)))(ti);
                g_min = min(g_min, min(v,[],'omitnan'));
                g_max = max(g_max, max(v,[],'omitnan'));
            end
            bh = 0.15;

            pos_s = []; neg_s = [];
            for s = 1:nS
                si = find(pmat(s,:) < 0.05);
                if ~isempty(si)
                    if mean(avg_mat(s, si), 'omitnan') > 0
                        pos_s(end+1) = s; %#ok<AGROW>
                    else
                        neg_s(end+1) = s; %#ok<AGROW>
                    end
                end
            end

            for p = 1:length(pos_s)
                s  = pos_s(p);
                si = find(pmat(s,:) < 0.05);
                yp = g_max + p * bh;
                app.drawSigBarsInteractive(ax, ta, si, yp, colors(s,:), ...
                    char(ms(s)), cohens_d_all.(char(ms(s))));
            end
            for n = 1:length(neg_s)
                s  = neg_s(n);
                si = find(pmat(s,:) < 0.05);
                yp = g_min - n * bh;
                app.drawSigBarsInteractive(ax, ta, si, yp, colors(s,:), ...
                    char(ms(s)), cohens_d_all.(char(ms(s))));
            end

            % Reference lines
            line(ax, [0 0], [-3 3], 'Color', 'r', 'LineStyle', '--', ...
                'LineWidth', 2, 'HandleVisibility', 'off');
            line(ax, fig_range, [0 0], 'Color', 'k', 'LineStyle', '--', ...
                'HandleVisibility', 'off');

            ax.XLim = fig_range;
            ax.YLim = [-3, 3];
            xlabel(ax, 'Latency (ms)', 'FontSize', 12);
            ylabel(ax, 'Relative Occurrence Frequency', 'FontSize', 12);
            if is_comp
                title(ax, sprintf('ROF: %s vs %s', upper(s1), upper(s2)), ...
                    'FontSize', 14);
            else
                title(ax, sprintf('ROF: %s', upper(s1)), 'FontSize', 14);
            end
            legend(ax, 'Location', 'northeast', 'FontSize', 10);
            grid(ax, 'on');
            hold(ax, 'off');
        end

        function drawSigBarsInteractive(~, ax, ta, sig_idx, y_pos, col, stateName, cohens_d_vec)
            % Draw clickable significance bars. Clicking shows microstate
            % name, time range, and average Cohen's d for that cluster.
            if isempty(sig_idx), return; end
            dd = diff(sig_idx);
            sp = find(dd > 1);
            starts = [sig_idx(1), sig_idx(sp + 1)];
            ends   = [sig_idx(sp),   sig_idx(end)];
            for seg = 1:length(starts)
                seg_idx = starts(seg):ends(seg);
                avg_d   = mean(abs(cohens_d_vec(seg_idx)), 'omitnan');
                t_start = ta(starts(seg));
                t_end   = ta(ends(seg));

                h = line(ax, [t_start, t_end], ...
                    [y_pos, y_pos], 'Color', col, 'LineWidth', 5, ...
                    'HandleVisibility', 'off', 'HitTest', 'on', ...
                    'PickableParts', 'all');

                % Store info for click callback
                tipStr = sprintf('Microstate %s\n%d – %d ms\nCohen''s d = %.3f', ...
                    stateName, round(t_start), round(t_end), avg_d);
                h.UserData = struct('tipStr', tipStr, ...
                    'x', (t_start + t_end) / 2, 'y', y_pos);
                h.ButtonDownFcn = @(src, ~) showSigBarTip(ax, src);
            end

            function showSigBarTip(axObj, src)
                % Remove any previous annotation from this axes
                old = findobj(axObj.Parent, 'Tag', 'SigBarAnnotation');
                delete(old);
                % Create a text annotation at the bar position
                text(axObj, src.UserData.x, src.UserData.y, ...
                    src.UserData.tipStr, ...
                    'FontSize', 11, 'FontWeight', 'bold', ...
                    'BackgroundColor', [1 1 0.85], ...
                    'EdgeColor', [0.3 0.3 0.3], ...
                    'Margin', 5, ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'bottom', ...
                    'Tag', 'SigBarAnnotation', ...
                    'HitTest', 'on', 'PickableParts', 'all', ...
                    'ButtonDownFcn', @(s,~) delete(s));
            end
        end

        function plotRTFHeatmap(app, result, ax)
            % Adapted from plot_transitions_t_scores_heatmap.m (into UIAxes)

            rtf = result.data;
            ms  = app.Microstates;
            nS  = length(ms);

            t_data = rtf.t_post_tms;
            h_data = rtf.h_post_tms;

            % Get Cohen's d if available
            if isfield(rtf, 'cohen_d_post_tms')
                d_data = rtf.cohen_d_post_tms;
            else
                d_data = NaN(nS, nS);
            end

            t_data(h_data == 0) = 0;
            for ii = 1:nS, t_data(ii, ii) = NaN; end

            imagesc(ax, t_data, 'AlphaData', ~isnan(t_data), ...
                'HitTest', 'off');
            colormap(ax, jet);
            cb = colorbar(ax);
            cb.Label.String = 'T-score';

            ax.XTick      = 1:nS;
            ax.XTickLabel = cellstr(ms);
            ax.YTick      = 1:nS;
            ax.YTickLabel = cellstr(ms);
            xlabel(ax, 'To State',   'FontSize', 12);
            ylabel(ax, 'From State', 'FontSize', 12);
            ax.FontSize = 11;

            % Interactive text annotations — click to show details
            for r = 1:nS
                for c = 1:nS
                    if ~isnan(t_data(r, c))
                        ht = text(ax, c, r, sprintf('%.2f', t_data(r, c)), ...
                            'HorizontalAlignment', 'center', ...
                            'FontSize', 10, 'Color', 'k', ...
                            'HitTest', 'on', 'PickableParts', 'all');

                        % Build tooltip info
                        fromSt = char(ms(r));
                        toSt   = char(ms(c));
                        d_val  = d_data(r, c);
                        isSig  = h_data(r, c) ~= 0;
                        if isSig
                            tipStr = sprintf('%s → %s\nT = %.3f\nCohen''s d = %.3f\np < 0.05 (significant)', ...
                                fromSt, toSt, t_data(r,c), d_val);
                        else
                            tipStr = sprintf('%s → %s\nT = %.3f\nCohen''s d = %.3f\nn.s.', ...
                                fromSt, toSt, rtf.t_post_tms(r,c), d_val);
                        end
                        ht.UserData = struct('tipStr', tipStr, 'col', c, 'row', r);
                        ht.ButtonDownFcn = @(src, ~) showRTFTip(ax, src);
                    end
                end
            end

            ps = app.PostStimStartField.Value;
            pe = app.PostStimEndRTFField.Value;
            if ~isempty(result.site2)
                title(ax, sprintf('RTF T-scores: %s vs %s  (Post-TMS: %d-%d ms)', ...
                    upper(result.site1), upper(result.site2), ps, pe), ...
                    'FontSize', 14);
            else
                title(ax, sprintf('RTF T-scores: %s  (Post-TMS: %d-%d ms)', ...
                    upper(result.site1), ps, pe), 'FontSize', 14);
            end

            axis(ax, 'equal', 'tight');

            function showRTFTip(axObj, src)
                old = findobj(axObj.Parent, 'Tag', 'RTFAnnotation');
                delete(old);
                text(axObj, src.UserData.col, src.UserData.row - 0.35, ...
                    src.UserData.tipStr, ...
                    'FontSize', 11, 'FontWeight', 'bold', ...
                    'BackgroundColor', [1 1 0.85], ...
                    'EdgeColor', [0.3 0.3 0.3], ...
                    'Margin', 5, ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'bottom', ...
                    'Tag', 'RTFAnnotation', ...
                    'HitTest', 'on', 'PickableParts', 'all', ...
                    'ButtonDownFcn', @(s,~) delete(s));
            end
        end

        % ----- Color utility ---------------------------------------------

        function colors = getStateColors(~, nStates)
            % Use MATLAB's built-in 'lines' colormap to match the original
            % plot_microstate_occurrence.m color scheme exactly.
            colors = lines(nStates);
        end
    end

    %% ===== UI Construction ==============================================
    methods (Access = private)

        function createComponents(app)

            % --- Main figure ---------------------------------------------
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [80, 80, 1280, 820];
            app.UIFigure.Name = 'TMS-EEG Microstate Analysis';
            app.UIFigure.Resize = 'on';
            app.UIFigure.AutoResizeChildren = 'off'; % grids handle resizing

            % --- Tab group (inside a grid layout for auto-resize) -----------
            mainGrid = uigridlayout(app.UIFigure, [1, 1]);
            mainGrid.Padding = [0 0 0 0];
            app.TabGroup = uitabgroup(mainGrid);

            app.createDataTab();
            app.createVizTab();
            app.createAnalysisTab();
            app.createResultsTab();

            % Tab navigation guard
            app.TabGroup.SelectionChangedFcn = ...
                createCallbackFcn(app, @TabGroupSelectionChanged, true);
        end

        % -----------------------------------------------------------------
        function createDataTab(app)
            app.DataTab       = uitab(app.TabGroup);
            app.DataTab.Title = '  Data & Settings  ';

            gl = uigridlayout(app.DataTab, [4, 2]);
            gl.RowHeight    = {'fit', 'fit', 'fit', '1x'};
            gl.ColumnWidth  = {'fit', '1x'};
            gl.Padding      = [20 20 20 20];
            gl.RowSpacing   = 15;

            % Load button
            app.LoadDataButton = uibutton(gl, 'push');
            app.LoadDataButton.Text      = '  Load Data File (.mat)  ';
            app.LoadDataButton.FontSize  = 14;
            app.LoadDataButton.FontWeight = 'bold';
            app.LoadDataButton.Layout.Row    = 1;
            app.LoadDataButton.Layout.Column = 1;
            app.LoadDataButton.ButtonPushedFcn = ...
                createCallbackFcn(app, @LoadDataButtonPushed, true);

            app.FileNameLabel = uilabel(gl);
            app.FileNameLabel.Text     = 'No file loaded';
            app.FileNameLabel.FontSize = 13;
            app.FileNameLabel.FontAngle = 'italic';
            app.FileNameLabel.Layout.Row    = 1;
            app.FileNameLabel.Layout.Column = 2;

            % --- Info panel ----------------------------------------------
            app.InfoPanel = uipanel(gl);
            app.InfoPanel.Title      = 'Data Information';
            app.InfoPanel.FontSize   = 13;
            app.InfoPanel.FontWeight = 'bold';
            app.InfoPanel.Layout.Row    = 2;
            app.InfoPanel.Layout.Column = [1, 2];

            ig = uigridlayout(app.InfoPanel, [4, 1]);
            ig.RowHeight = {'fit', 'fit', 'fit', 'fit'};
            ig.Padding   = [12 8 12 8];

            app.NumSubjectsLabel           = uilabel(ig);
            app.NumSubjectsLabel.Text      = 'Subjects: --';
            app.NumSubjectsLabel.FontSize  = 12;

            app.StimSitesLabel             = uilabel(ig);
            app.StimSitesLabel.Text        = 'Stimulation Sites: --';
            app.StimSitesLabel.FontSize    = 12;

            app.MicrostatesDetectedLabel         = uilabel(ig);
            app.MicrostatesDetectedLabel.Text    = 'Microstates: --';
            app.MicrostatesDetectedLabel.FontSize = 12;

            app.NumTrialsLabel             = uilabel(ig);
            app.NumTrialsLabel.Text        = 'Trials: --';
            app.NumTrialsLabel.FontSize    = 12;

            % --- Settings panel ------------------------------------------
            app.SettingsPanel = uipanel(gl);
            app.SettingsPanel.Title      = 'Settings';
            app.SettingsPanel.FontSize   = 13;
            app.SettingsPanel.FontWeight = 'bold';
            app.SettingsPanel.Layout.Row    = 3;
            app.SettingsPanel.Layout.Column = [1, 2];

            sg = uigridlayout(app.SettingsPanel, [3, 2]);
            sg.RowHeight    = {'fit', 'fit', 'fit'};
            sg.ColumnWidth  = {'fit', '1x'};
            sg.Padding      = [12 8 12 8];

            lbl = uilabel(sg); lbl.Text = 'TMS Onset Index:';   lbl.FontSize = 12;
            app.TMSIndexSpinner = uispinner(sg);
            app.TMSIndexSpinner.Value  = 1001;
            app.TMSIndexSpinner.Limits = [1, Inf];
            app.TMSIndexSpinner.Step   = 1;
            app.TMSIndexSpinner.FontSize = 12;
            app.TMSIndexSpinner.ValueChangedFcn = ...
                createCallbackFcn(app, @TMSIndexChanged, true);

            lbl2 = uilabel(sg); lbl2.Text = 'Sample Rate (Hz):'; lbl2.FontSize = 12;
            app.SampleRateSpinner = uispinner(sg);
            app.SampleRateSpinner.Value  = 1000;
            app.SampleRateSpinner.Limits = [1, Inf];
            app.SampleRateSpinner.Step   = 1;
            app.SampleRateSpinner.FontSize = 12;
            app.SampleRateSpinner.ValueChangedFcn = ...
                createCallbackFcn(app, @SampleRateChanged, true);

            lbl3 = uilabel(sg); lbl3.Text = 'Microstate Labels (comma-separated):'; lbl3.FontSize = 12;
            app.MicrostateLabelsField = uieditfield(sg, 'text');
            app.MicrostateLabelsField.Value    = 'A, B, C, D, E';
            app.MicrostateLabelsField.FontSize = 12;
            app.MicrostateLabelsField.ValueChangedFcn = ...
                createCallbackFcn(app, @MicrostateLabelsChanged, true);
        end

        % -----------------------------------------------------------------
        function createVizTab(app)
            app.VizTab       = uitab(app.TabGroup);
            app.VizTab.Title = '  Visualization  ';

            gl = uigridlayout(app.VizTab, [3, 6]);
            gl.RowHeight    = {'fit', '0.15x', '1x'};
            gl.ColumnWidth  = {'fit', '1x', 'fit', '1x', 'fit', '1x'};
            gl.Padding      = [15 12 15 12];
            gl.RowSpacing   = 8;

            % Subject
            lbl = uilabel(gl); lbl.Text = 'Subject:'; lbl.FontSize = 12;
            lbl.Layout.Row = 1; lbl.Layout.Column = 1;
            app.SubjectDropDown = uidropdown(gl);
            app.SubjectDropDown.Items    = {'--'};
            app.SubjectDropDown.FontSize = 12;
            app.SubjectDropDown.Layout.Row = 1; app.SubjectDropDown.Layout.Column = 2;
            app.SubjectDropDown.ValueChangedFcn = ...
                createCallbackFcn(app, @VizSelectionChanged, true);

            % Site
            lbl2 = uilabel(gl); lbl2.Text = 'Stim. Site:'; lbl2.FontSize = 12;
            lbl2.Layout.Row = 1; lbl2.Layout.Column = 3;
            app.SiteDropDown = uidropdown(gl);
            app.SiteDropDown.Items    = {'--'};
            app.SiteDropDown.FontSize = 12;
            app.SiteDropDown.Layout.Row = 1; app.SiteDropDown.Layout.Column = 4;
            app.SiteDropDown.ValueChangedFcn = ...
                createCallbackFcn(app, @VizSelectionChanged, true);

            % Trial
            lbl3 = uilabel(gl); lbl3.Text = 'Trial:'; lbl3.FontSize = 12;
            lbl3.Layout.Row = 1; lbl3.Layout.Column = 5;
            app.TrialDropDown = uidropdown(gl);
            app.TrialDropDown.Items    = {'--'};
            app.TrialDropDown.FontSize = 12;
            app.TrialDropDown.Layout.Row = 1; app.TrialDropDown.Layout.Column = 6;
            app.TrialDropDown.ValueChangedFcn = ...
                createCallbackFcn(app, @TrialChanged, true);

            % Ribbon axes
            app.RibbonAxes = uiaxes(gl);
            app.RibbonAxes.Layout.Row    = 2;
            app.RibbonAxes.Layout.Column = [1, 6];
            title(app.RibbonAxes, '');
            app.RibbonAxes.YTick = [];

            % Main sequence axes
            app.SequenceAxes = uiaxes(gl);
            app.SequenceAxes.Layout.Row    = 3;
            app.SequenceAxes.Layout.Column = [1, 6];
            xlabel(app.SequenceAxes, 'Time (ms)');
            ylabel(app.SequenceAxes, 'Microstate');
            title(app.SequenceAxes, 'Microstate Sequence');
        end

        % -----------------------------------------------------------------
        function createAnalysisTab(app)
            app.AnalysisTab       = uitab(app.TabGroup);
            app.AnalysisTab.Title = '  Analysis  ';

            % Main layout: 5 rows x 2 cols
            %   Row 1: Analysis Type  | Site Selection
            %   Row 2: Features       | Statistics
            %   Row 3: Time Parameters (span 2 cols)
            %   Row 4: Note about ignore window
            %   Row 5: spacer
            %   Row 6: Run button + status (span 2 cols)
            gl = uigridlayout(app.AnalysisTab, [6, 2]);
            gl.RowHeight   = {85, 'fit', 'fit', 'fit', '1x', 'fit'};
            gl.ColumnWidth = {'1x', '1x'};
            gl.Padding     = [20 15 20 15];
            gl.RowSpacing  = 8;
            gl.ColumnSpacing = 14;

            % ============================================================
            % Row 1, Col 1: Analysis Type
            % ============================================================
            app.AnalysisTypeGroup = uibuttongroup(gl);
            app.AnalysisTypeGroup.Title      = 'Analysis Type';
            app.AnalysisTypeGroup.FontSize   = 14;
            app.AnalysisTypeGroup.FontWeight = 'bold';
            app.AnalysisTypeGroup.Layout.Row    = 1;
            app.AnalysisTypeGroup.Layout.Column = 1;
            app.AnalysisTypeGroup.SelectionChangedFcn = ...
                createCallbackFcn(app, @AnalysisTypeChanged, true);

            app.SingleSiteRadio = uiradiobutton(app.AnalysisTypeGroup);
            app.SingleSiteRadio.Text     = 'Single Stimulation Site';
            app.SingleSiteRadio.FontSize = 13;
            app.SingleSiteRadio.Position = [12, 22, 300, 22];
            app.SingleSiteRadio.Value    = true;

            app.BetweenSiteRadio = uiradiobutton(app.AnalysisTypeGroup);
            app.BetweenSiteRadio.Text     = 'Between Two Stimulation Sites';
            app.BetweenSiteRadio.FontSize = 13;
            app.BetweenSiteRadio.Position = [12, 2, 300, 22];

            % ============================================================
            % Row 1, Col 2: Site Selection
            % ============================================================
            sitePanel = uipanel(gl);
            sitePanel.Title      = 'Site Selection';
            sitePanel.FontSize   = 14;
            sitePanel.FontWeight = 'bold';
            sitePanel.Layout.Row    = 1;
            sitePanel.Layout.Column = 2;

            sg = uigridlayout(sitePanel, [2, 2]);
            sg.RowHeight   = {'1x', '1x'};
            sg.ColumnWidth = {'fit', '1x'};
            sg.Padding     = [10 6 10 6];

            uilabel(sg, 'Text', 'Site 1:', 'FontSize', 13);
            app.AnalysisSite1DropDown = uidropdown(sg);
            app.AnalysisSite1DropDown.Items    = {'--'};
            app.AnalysisSite1DropDown.FontSize = 13;

            app.AnalysisSite2Label = uilabel(sg);
            app.AnalysisSite2Label.Text     = 'Site 2:';
            app.AnalysisSite2Label.FontSize = 13;
            app.AnalysisSite2Label.Enable   = 'off';

            app.AnalysisSite2DropDown = uidropdown(sg);
            app.AnalysisSite2DropDown.Items    = {'--'};
            app.AnalysisSite2DropDown.FontSize = 13;
            app.AnalysisSite2DropDown.Enable   = 'off';

            % ============================================================
            % Row 2, Col 1: Features to Analyze
            % ============================================================
            featPanel = uipanel(gl);
            featPanel.Title      = 'Features to Analyze';
            featPanel.FontSize   = 14;
            featPanel.FontWeight = 'bold';
            featPanel.Layout.Row    = 2;
            featPanel.Layout.Column = 1;

            fg = uigridlayout(featPanel, [2, 1]);
            fg.RowHeight = {'fit', 'fit'};
            fg.Padding   = [10 8 10 8];

            app.ROFCheckBox = uicheckbox(fg);
            app.ROFCheckBox.Text     = 'ROF (Relative Occurrence Frequencies)';
            app.ROFCheckBox.Value    = true;
            app.ROFCheckBox.FontSize = 13;

            app.RTFCheckBox = uicheckbox(fg);
            app.RTFCheckBox.Text     = 'RTF (Relative Transition Frequencies)';
            app.RTFCheckBox.Value    = true;
            app.RTFCheckBox.FontSize = 13;

            % ============================================================
            % Row 2, Col 2: Statistical Parameters
            % ============================================================
            statPanel = uipanel(gl);
            statPanel.Title      = 'Statistical Parameters';
            statPanel.FontSize   = 14;
            statPanel.FontWeight = 'bold';
            statPanel.Layout.Row    = 2;
            statPanel.Layout.Column = 2;

            stg = uigridlayout(statPanel, [2, 2]);
            stg.RowHeight   = {'fit', 'fit'};
            stg.ColumnWidth = {'fit', '1x'};
            stg.Padding     = [10 8 10 8];

            uilabel(stg, 'Text', 'Permutations:', 'FontSize', 13);
            app.NPermSpinner = uispinner(stg);
            app.NPermSpinner.Value    = 5000;
            app.NPermSpinner.Limits   = [100, 100000];
            app.NPermSpinner.Step     = 500;
            app.NPermSpinner.FontSize = 13;

            uilabel(stg, 'Text', 'Correction:', 'FontSize', 13);
            app.CorrectionDropDown = uidropdown(stg);
            app.CorrectionDropDown.Items    = {'Bonferroni', 'FDR', 'None'};
            app.CorrectionDropDown.Value    = 'Bonferroni';
            app.CorrectionDropDown.FontSize = 13;

            % ============================================================
            % Row 3: Time Parameters (full width, spans 2 cols)
            % ============================================================
            timePanel = uipanel(gl);
            timePanel.Title      = 'Time Parameters (ms, relative to TMS)';
            timePanel.FontSize   = 14;
            timePanel.FontWeight = 'bold';
            timePanel.Layout.Row    = 3;
            timePanel.Layout.Column = [1, 2];

            tg = uigridlayout(timePanel, [2, 6]);
            tg.RowHeight   = {'fit', 'fit'};
            tg.ColumnWidth = {'fit', '1x', 'fit', '1x', 'fit', '1x'};
            tg.Padding     = [12 8 12 8];
            tg.ColumnSpacing = 14;

            % Row 1: Baseline — Start | End
            uilabel(tg, 'Text', 'Baseline Start:', 'FontSize', 13);
            app.BaselineStartField = uieditfield(tg, 'numeric');
            app.BaselineStartField.Value = -1000; app.BaselineStartField.FontSize = 13;

            uilabel(tg, 'Text', 'Baseline End:', 'FontSize', 13);
            app.BaselineEndField = uieditfield(tg, 'numeric');
            app.BaselineEndField.Value = -10; app.BaselineEndField.FontSize = 13;

            uilabel(tg, 'Text', ''); % placeholder
            uilabel(tg, 'Text', ''); % placeholder

            % Row 2: Post-Stim — Start | End (ROF) | End (RTF)
            uilabel(tg, 'Text', 'Post-Stim Start:', 'FontSize', 13);
            app.PostStimStartField = uieditfield(tg, 'numeric');
            app.PostStimStartField.Value = 20; app.PostStimStartField.FontSize = 13;

            uilabel(tg, 'Text', 'Post-Stim End (ROF):', 'FontSize', 13);
            app.PostStimEndROFField = uieditfield(tg, 'numeric');
            app.PostStimEndROFField.Value = 1000; app.PostStimEndROFField.FontSize = 13;

            uilabel(tg, 'Text', 'Post-Stim End (RTF):', 'FontSize', 13);
            app.PostStimEndRTFField = uieditfield(tg, 'numeric');
            app.PostStimEndRTFField.Value = 500; app.PostStimEndRTFField.FontSize = 13;

            % Note under time panel: ignore window is derived
            noteLabel = uilabel(gl);
            noteLabel.Text = '   Note: The ignore window around TMS is automatically defined as [Baseline End .. Post-Stim Start].';
            noteLabel.FontSize   = 12;
            noteLabel.FontAngle  = 'italic';
            noteLabel.FontColor  = [0.4 0.4 0.4];
            noteLabel.Layout.Row    = 4;
            noteLabel.Layout.Column = [1, 2];

            % ============================================================
            % Row 5: flexible spacer (pushes button to bottom)
            % ============================================================

            % ============================================================
            % Row 6: Centered Run button (spans 2 cols)
            % ============================================================
            runGrid = uigridlayout(gl, [1, 3]);
            runGrid.ColumnWidth = {'1x', 'fit', '1x'};
            runGrid.RowHeight   = {'fit'};
            runGrid.Padding     = [0 4 0 4];
            runGrid.Layout.Row    = 6;
            runGrid.Layout.Column = [1, 2];

            uilabel(runGrid, 'Text', ''); % left spacer

            app.RunAnalysisButton = uibutton(runGrid, 'push');
            app.RunAnalysisButton.Text            = '  Run Analysis  ';
            app.RunAnalysisButton.FontSize        = 16;
            app.RunAnalysisButton.FontWeight      = 'bold';
            app.RunAnalysisButton.FontColor       = [1 1 1];
            app.RunAnalysisButton.BackgroundColor = [0.18, 0.55, 0.34];
            app.RunAnalysisButton.ButtonPushedFcn = ...
                createCallbackFcn(app, @RunAnalysisButtonPushed, true);

            uilabel(runGrid, 'Text', ''); % right spacer

            % Hidden status label (kept for analysis callbacks)
            app.AnalysisStatusLabel = uilabel(gl);
            app.AnalysisStatusLabel.Text    = '';
            app.AnalysisStatusLabel.Visible = 'off';
            app.AnalysisStatusLabel.Layout.Row    = 6;
            app.AnalysisStatusLabel.Layout.Column = 1;
        end

        % -----------------------------------------------------------------
        function createResultsTab(app)
            app.ResultsTab       = uitab(app.TabGroup);
            app.ResultsTab.Title = '  Results  ';

            gl = uigridlayout(app.ResultsTab, [2, 4]);
            gl.RowHeight   = {'fit', '1x'};
            gl.ColumnWidth = {'fit', '1x', 'fit', 'fit'};
            gl.Padding     = [15 12 15 12];
            gl.RowSpacing  = 10;

            uilabel(gl, 'Text', 'View Result:', 'FontSize', 12);
            app.ResultDropDown = uidropdown(gl);
            app.ResultDropDown.Items    = {'--'};
            app.ResultDropDown.FontSize = 12;
            app.ResultDropDown.ValueChangedFcn = ...
                createCallbackFcn(app, @ResultDropDownChanged, true);

            app.ExportFigureButton = uibutton(gl, 'push');
            app.ExportFigureButton.Text     = 'Export Figure';
            app.ExportFigureButton.FontSize = 12;
            app.ExportFigureButton.ButtonPushedFcn = ...
                createCallbackFcn(app, @ExportFigureButtonPushed, true);

            app.ExportResultsButton = uibutton(gl, 'push');
            app.ExportResultsButton.Text     = 'Export Results (.mat)';
            app.ExportResultsButton.FontSize = 12;
            app.ExportResultsButton.ButtonPushedFcn = ...
                createCallbackFcn(app, @ExportResultsButtonPushed, true);

            app.ResultAxes = uiaxes(gl);
            app.ResultAxes.Layout.Row    = 2;
            app.ResultAxes.Layout.Column = [1, 4];
            title(app.ResultAxes, 'Results will appear here after analysis');
            app.ResultAxes.FontSize = 11;
        end
    end

    %% ===== Public Methods (Constructor / Destructor) ====================
    methods (Access = public)

        function app = TMS_EEG_Microstates_App
            % Construct the app
            app.DataLoaded      = false;
            app.VizEnabled      = false;
            app.AnalysisEnabled = false;
            app.ResultsEnabled  = false;
            app.AnalysisResults = struct( ...
                'name', {}, 'type', {}, 'data', {}, ...
                'site1', {}, 'site2', {}, ...
                'p_values_matrix', {}, 'data_struct_used', {});

            app.createComponents();
            registerApp(app, app.UIFigure);
            app.UIFigure.Visible = 'on';
        end

        function delete(app)
            delete(app.UIFigure);
        end
    end
end
