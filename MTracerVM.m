classdef MTracerVM < handle
    % 
    
    properties
        % UI
        app;
        mapFig;
        mapAxes;
        mapLayers struct = struct;
        
        % Data
        recId char = 'NP0_B0';
        
        rezFile;                    % path of rez.mat file
        rezOps struct;              % kilosort rez
        lfBinFile;                  % path of lf.bin file
        lfMeta struct;              % LFP metadata loaded from lf.meta
        channelTable table;         % channel info
        
        se MSessionExplorer;        % current data container, referenced to either seo or sek
        seo MSessionExplorer;       % container of original data
        sek MSessionExplorer;       % container of motion corrected data
        
        traces MTracerTrace;        % trace objects
        
        % Runtime variables
        focus = [0 0];
        currentTrace = NaN;
        isTracing = false;
    end
    
    properties(Dependent)
        hasLFP;
        hasSpikes;
        hasTrace;
        hasApp;
        hasMapAxes;
        layerNames;
        tLims;
        yLims;
    end
    
    methods
        function val = get.hasLFP(this)
            val = ismember('LFP', this.se.tableNames);
        end
        function val = get.hasSpikes(this)
            val = ismember('spikes', this.se.tableNames);
        end
        function val = get.hasTrace(this)
            val = ~isempty(this.traces);
        end
        function val = get.hasApp(this)
            val = ~isempty(this.app);
        end
        function val = get.hasMapAxes(this)
            h = this.mapAxes;
            val = ~isempty(h) && ishandle(h) && isvalid(h);
        end
        function val = get.layerNames(this)
            val = fieldnames(this.mapLayers);
        end
        function val = get.tLims(this)
            if this.hasSpikes
                val = [this.rezOps.tstart this.rezOps.tend] / this.rezOps.fs;
            elseif this.hasLFP
                val = [0 str2double(this.lfpMeta.fileTimeSecs)];
            else
                val = [0 1e4];
            end
        end
        function val = get.yLims(this)
            if ~isempty(this.channelTable)
                val = [min(this.channelTable.ycoords) max(this.channelTable.ycoords)];
            else
                val = [0 7660];
            end
        end
    end
    
    methods
        % Construction
        function this = MTracerVM(app)
            % Constructor
            
            if nargin > 0
                this.app = app;
            end
            
            this.mapLayers.focus = [];
            this.mapLayers.spikes = [];
            this.mapLayers.LFP = [];
            this.mapLayers.ap = [];
            this.mapLayers.anchors = [];
            this.mapLayers.interp = [];
            
            this.se = MSessionExplorer();
        end
        
        function LoadKilosortRez(this, filePath)
            % Load data from Kilosort output rez.mat
            
            if nargin < 2 || isempty(filePath)
                filePath = MBrowse.File([], 'Please select a rez.mat file', '*.mat');
            end
            if ~exist(filePath, 'file')
                return;
            end
            
            try
                % Load rez.mat
                load(filePath);
                this.rezFile = filePath;
                this.rezOps = rez.ops;
                
                % Make a table of spikes
                tb = table();
                tb.time = rez.st0(:,1) / rez.ops.fs;
                tb.y = rez.st0(:,2);
                tb.amp = rez.st0(:,3);
                tb = sortrows(tb, {'time', 'y'});
                
                % Save table to se
                C = mat2cell(tb{:,:}, height(tb), ones(1,width(tb)));
                tb = cell2table(C, 'VariableNames', tb.Properties.VariableNames);
                this.se.SetTable('spikes', tb, 'timeSeries', 0);
                
            catch e
                assignin('base', 'e', e);
                disp(e);
            end
        end
        
        function isLoaded = LoadLFP(this, filePath)
            % Load image
            
            isLoaded = false;
            
            if nargin < 2 || isempty(filePath)
                filePath = MBrowse.File([], 'Please select a lf.bin file', '*.bin');
            end
            if ~exist(filePath, 'file')
                return;
            end
            
            try
                % Load LFP data
                this.lfBinFile = filePath;
                lfMetaFile = strrep(filePath, '.bin', '.meta');
                [this.lfMeta, lfp, lfpTime] = MSpikeGLX.ReadLFP(lfMetaFile);
                
                % Sort channels by depth
                chanTb = this.LoadChannelConfig();
                [chanTb, I] = sortrows(chanTb, 'ycoords');
                this.channelTable = chanTb;
                lfp = lfp(:,I);
                
                % Make downsampled LFP timeSeries table
                lfpTb = table();
                r = 50; % downsample 50x from 2500Hz to 50Hz
                lfpTb.time{1} = downsample(lfpTime, r, r/2);
                lfpTb.v{1} = MMath.Decimate(lfp, r, r/2);
                this.se.SetTable('LFP', lfpTb, 'timeSeries', 0);
                
                isLoaded = true;
                
            catch e
                assignin('base', 'e', e);
                disp(e);
            end
        end
        
        function chanTb = LoadChannelConfig(this)
            % Load the 354-channel configuration
            s = load('NP1_NHP_HalfCol_kilosortChanMap.mat');
            chanTb = table;
            chanTb.id = s.chanMap;
            chanTb.xcoords = s.xcoords;
            chanTb.ycoords = s.ycoords;
        end
        
        function InitializeTrace(this)
            % 
            this.traces(end+1,1) = MTracerTrace(this);
            this.currentTrace = numel(this.traces);
        end
        
        function SelectTrace(this, t ,y)
            % 
            
            % No change in selection during tracing
            if this.isTracing || isempty(this.traces)
                return
            end
            
            % Find y-distance from traces to reference point
            yd = arrayfun(@(x) x.Dist2Point(t,y), this.traces);
            
            % Find trace with shortest distance and within tolerance
            [ydMin, I] = min(yd);
            if ydMin > 50
                this.currentTrace = NaN;
            else
                this.currentTrace = I;
            end
        end
        
        function DeleteTrace(this, ind)
            % Delete the selected trace
            if nargin < 2
                ind = this.currentTrace;
            end
            if ~isnan(ind)
                this.traces(ind) = [];
                this.currentTrace = NaN;
                this.ClearInvalidHandles();
                this.PlotTraces();
            end
            if this.hasApp
                this.app.TraceListBox.Items = arrayfun(@(x) x.dispName, this.traces, 'Uni', false);
            end
        end
        
        function delete(this)
            delete(this.mapFig);
        end
        
        % Display
        function InitializeMapAxes(this)
            
            % Create main figure if absent
            if isempty(this.mapFig) || ~isvalid(this.mapFig)
                this.mapFig = MPlot.Figure( ...
                    'Name', 'MTracer: Maps', ...
                    'NumberTitle', 'off', ...
                    'IntegerHandle', 'off', ...
                    'WindowKeyPressFcn', @this.KeyPress, ...
                    'WindowKeyReleaseFcn', @this.KeyRelease, ...
                    'WindowScrollWheelFcn', @this.Scroll, ...
                    'Interruptible', 'off', ...
                    'BusyAction', 'cancel');
            end
            
            % Delete old map axes
            if ishandle(this.mapAxes)
                delete(this.mapAxes);
            end
            
            % Create and format new map axes
            ax = axes(this.mapFig);
            ax.XLabel.String = 'Time (sec)';
            ax.YLabel.String = 'Distance from tip (um)';
            ax.XLim = this.tLims;
            ax.YLim = this.yLims;
            ax.LooseInset = [0 0 0 0];
            ax.ButtonDownFcn = @this.SetPoint;
            ax.BusyAction = 'cancel';
            hold(ax, 'on');
            MPlot.Axes(ax);
            this.mapAxes = ax;
        end
        
        function PlotFocus(this)
            % Plot the focus indicator
            
            if ~this.hasMapAxes
                return
            end
            
            xx = [this.focus([1 1])' this.tLims'];
            yy = [this.yLims' this.focus([2 2])'];
            
            delete(this.mapLayers.focus);
            
            this.mapLayers.focus = plot(xx, yy, '-', 'Color', [0 1 0 .5], 'LineWidth', 1, 'HitTest', 'off');
        end
        
        function PlotLFP(this)
            % Plot LFP siganl as a heatmap
            
            if ~this.hasLFP || ~this.hasMapAxes
                return
            end
            
            tb = this.se.GetTable('LFP');
            C = tb.v{1}';
            x = tb.time{1};
            y = this.channelTable.ycoords;
            
            delete(this.mapLayers.LFP)
            this.mapLayers.LFP = imagesc(this.mapAxes, x, y, C, 'HitTest', 'off');
            colormap(this.mapAxes, MPlot.PolarMap);
            
            MPlot.Axes(this.mapAxes);
            this.mapAxes.XLim = this.tLims;
            this.mapAxes.YLim = this.yLims;
            this.mapAxes.YDir = 'normal';
        end
        
        function PlotSpikes(this)
            % Plot spike events as a function of time and depth
            
            if ~this.hasSpikes || ~this.hasMapAxes
                return
            end
            
            tb = this.se.GetTable('spikes');
            t = tb.time{1};
            y = tb.y{1};
            a = tb.amp{1};
            
            ampRange = 8 : 100;
            hh = cell(numel(ampRange), 1);
            for k = 1 : numel(ampRange)
                % for each amplitude bin, plot all the spikes of that size in the
                % same shade of gray
                ind = a == ampRange(k); % the amplitudes are rounded to integers
                if ~any(ind)
                    continue
                end
                hh{k} = plot(this.mapAxes, t(ind), y(ind), ...
                    'LineStyle', 'none', ...
                    'Marker', '.', ...
                    'MarkerSize', 4, ...
                    'Color', [1 1 1] * max(0, 1-ampRange(k)/40), ... % the marker color here has been carefully tuned
                    'HitTest', 'off');
            end
            delete(this.mapLayers.spikes);
            this.mapLayers.spikes = cat(1, hh{:});
            
            MPlot.Axes(this.mapAxes);
            this.mapAxes.XLim = this.tLims;
            this.mapAxes.YLim = this.yLims;
            this.mapAxes.YDir = 'normal';
        end
        
        function PlotTraces(this)
            % Plot traces
            if ~this.hasMapAxes
                return
            end
            for i = 1 : numel(this.traces)
                this.traces(i).PlotTrace(ismember(i, this.currentTrace));
            end
        end
        
        function ToggleLayer(this, num)
            % Toggle layer on and off
            
            if num < 1 || num > numel(this.layerNames)
                return
            end
            
            mapLayer = this.mapLayers.(this.layerNames{num});
            
            for i = 1 : numel(mapLayer)
                mapLayer(i).Visible = ~mapLayer(i).Visible;
            end
        end
        
        function ClearInvalidHandles(this)
            % Remove invalid object handles in every layers
            n = fieldnames(this.mapLayers);
            for i = 1 : numel(n)
                if ~isempty(this.mapLayers.(n{i}))
                    isValid = isvalid(this.mapLayers.(n{i}));
                    this.mapLayers.(n{i})(~isValid) = [];
                end
            end
        end
        
        % Interactions
        function KeyPress(this, src, eventdata)
            
            K = eventdata.Key;
            M = eventdata.Modifier;
            
            if ismember('control', M)
                switch K
                    case 's'
                        % Save settings and traces
                        
                    case 'z'
                        % Undo/Redo point
                        if this.isTracing && ~isnan(this.currentTrace)
                            if ismember('shift', M)
                                this.traces(this.currentTrace).Redo();
                            else
                                this.traces(this.currentTrace).Undo();
                            end
                        end
                        
                    otherwise
                end
            else
                % Toggle layer in map
                if isscalar(K) && K >= '1' && K <= '9'
                    this.ToggleLayer(str2double(K))
                    pause(0.2);
                end
                
                % Change time and depth windows
                t = this.focus(1);
                y = this.focus(2);
                tWin = this.mapAxes.XLim;
                yWin = this.mapAxes.YLim;
                dt = diff(tWin) / 50;
                dy = diff(yWin) / 50;
                if ismember('shift', M)
                    dt = dt * 5;
                    dy = dy * 2;
                end
                
                switch K
                    case 't'
                        this.EnterTracingMode();
                    case {'f', 'escape'}
                        this.isTracing = false;
                        disp('Navigation Mode');
                        
                    case {'a', 'leftarrow'}
                        % Forward in time
                        tWin = tWin - dt;
                        this.UpdateMapLims(tWin, yWin);
                    case {'d', 'rightarrow'}
                        % Backward in time
                        tWin = tWin + dt;
                        this.UpdateMapLims(tWin, yWin);
                    case {'s', 'downarrow'}
                        % Downward in depth
                        yWin = yWin - dy;
                        this.UpdateMapLims(tWin, yWin);
                    case {'w', 'uparrow'}
                        % Upward in depth
                        yWin = yWin + dy;
                        this.UpdateMapLims(tWin, yWin);
                        
                    case 'leftbracket'
                        % Zoom-in in time
                        tWin = (tWin-t)/0.5 + t;
                        this.UpdateMapLims(tWin, yWin, t, y);
                    case 'rightbracket'
                        % Zoom-out in time
                        tWin = (tWin-t)*0.5 + t;
                        this.UpdateMapLims(tWin, yWin, t, y);
                    case {'equal', 'add'}
                        % Zoom-in in depth
                        yWin = (yWin-y)*0.8 + y;
                        this.UpdateMapLims(tWin, yWin, t, y);
                    case {'hyphen', 'subtract'}
                        % Zoom-out in depth
                        yWin = (yWin-y)/0.8 + y;
                        this.UpdateMapLims(tWin, yWin, t, y);
                        
                    case 'c'
                        % Center the reference point
                        this.UpdateMapLims(tWin, yWin, t, y);
                    case 'o'
                        % View full map
                        this.UpdateMapLims(this.tLims, this.yLims);
                        
                    case 'period'
                        % Increase spike marker size
                        this.ChangeSpikeMarkerSize(2, 'relative');
                    case 'comma'
                        % Decrease spike marker size
                        this.ChangeSpikeMarkerSize(0.5, 'relative');
                        
                    otherwise
                end
            end
            
%             disp(K);
        end
        
        function KeyRelease(this, src, eventdata)
            
            K = eventdata.Key;
            M = eventdata.Modifier;
            
%             if ismember('control', M)
%                 switch K
%                     case 's'
%                         this.SaveTraces();
%                     otherwise
%                 end
%             end
        end
        
        function Scroll(this, src, eventdata)
            tWin = this.mapAxes.XLim;
            dt = eventdata.VerticalScrollCount * eventdata.VerticalScrollAmount/3 * diff(tWin) / 50;
            this.UpdateMapLims(tWin+dt, this.mapAxes.YLim);
        end
        
        function UpdateMapLims(this, tWin, yWin, t, y)
            
            if nargin > 3
                tWin = tWin - mean(tWin) + t;
                yWin = yWin - mean(yWin) + y;
            end
            
            if tWin(1) < this.tLims(1)
                tWin = [this.tLims(1), this.tLims(1)+diff(tWin)];
            end
            if tWin(2) > this.tLims(2)
                tWin = [this.tLims(2)-diff(tWin), this.tLims(2)];
            end
            tWin = MMath.Bound(tWin, this.tLims);
            
            if yWin(1) < this.yLims(1)
                yWin = [this.yLims(1) this.yLims(1)+diff(yWin)];
            end
            if yWin(2) > this.yLims(2)
                yWin = [this.yLims(2)-diff(yWin), this.yLims(2)];
            end
            yWin = MMath.Bound(yWin, this.yLims);
            
            set(this.mapAxes, 'XLim', tWin, 'YLim', yWin);
        end
        
        function ChangeSpikeMarkerSize(this, r, valType)
            if ~this.hasSpikes || ~isfield(this.mapLayers, 'spikes')
                return
            end
            hh = this.mapLayers.spikes;
            sz = hh(1).MarkerSize;
            if strcmp(valType, 'relative')
                sz = sz * r;
            else
                sz = r;
            end
            sz = MMath.Bound(sz, [2 16]);
            for i = 1 : numel(hh)
                hh(i).MarkerSize = sz;
            end
        end
        
        function EnterTracingMode(this)
            if isscalar(this.currentTrace)
                this.isTracing = true;
                disp('Tracing Mode');
            else
                error('Cannot enter tracing mode with more than one selected trace\n');
            end
        end
        
        function SetPoint(this, src, eventdata)
            % Place a point on image
            
            % Read mouse position
            mousePos = get(this.mapAxes, 'CurrentPoint');
            
            % Update reference point
            t = MMath.Bound(mousePos(1), this.tLims);
            y = MMath.Bound(mousePos(3), this.yLims);
            this.focus = [t y];
            this.PlotFocus();
            
            if ~this.isTracing
                % Trace selection
                yd = arrayfun(@(x) x.Dist2Selection(t, y), this.traces);
                if all(isnan(yd))
                    this.currentTrace = NaN;
                else
                    [~, this.currentTrace] = min(yd);
                end
                this.PlotTraces();
                if this.hasApp
                    this.app.SelectTrace(this.currentTrace);
                end
            else
                % Adding or remove point
                if isnan(this.currentTrace)
                    set(this.mapLayers.anchors, 'Visible', 'on');
                    set(this.mapLayers.interp, 'Visible', 'on');
                    this.InitializeTrace();
                end
                tr = this.traces(this.currentTrace);
                if eventdata.Button == 1
                    tr.SetPoint(t, y);
                elseif eventdata.Button == 3
                    tr.RemovePoint(t, y);
                end
            end
        end
        
        function AddHandle2Layer(this, h, layerName)
            % Add a plot object handle to a layer
            
            if isempty(this.mapLayers.(layerName))
                % Add object to empty layer
                this.mapLayers.(layerName) = h;
            else
                % Append object handle to the end
                this.mapLayers.(layerName)(end+1,1) = h;
                this.ClearInvalidHandles();
            end
        end
        
        % Save and load traces
        function SaveTraces(this)
            % Save the data of each trace
            
            if ~numel(this.traces)
                uiwait(msgbox('No trace to save.', 'Save', 'modal'));
                return;
            end
            
            folderName = ['MTracer_' this.recId '_' datestr(now, 'yyyy-mm-dd_HH-MM-SS')];
            mkdir(folderName);
            
            for i = 1 : numel(this.traces)
                dataTb = this.traces(i).dataTb;
                save(fullfile(folderName, [this.traces(i).fileName '.mat']), 'dataTb');
            end
            
            uiwait(msgbox('Tracing results are saved successfully.', 'Save', 'modal'));
        end
        
        function LoadTraces(this)
            % Load saved traces
            
            filePaths = MBrowse.Files([], 'Please select one or more saved traces', '.mat');
            if isempty(filePaths)
                return
            end
            
            for i = 1 : numel(filePaths)
                try
                    load(filePaths{i}, '-mat');
                    obj = MTracerTrace(this, dataTb);
                    this.traces(end+1,1) = obj;
                    this.currentTrace = numel(this.traces);
                    disp(['Loaded trace: ' obj.dispName]);
                catch e
                    warning('Error occured when loading: \n%s', filePaths{i});
                    disp(e);
                end
            end
            
            this.PlotTraces();
            if this.hasApp && this.hasTrace
                this.app.TraceListBox.Items = arrayfun(@(x) x.dispName, this.traces, 'Uni', false);
            end
        end
        
        % Data operations
        function [CData, x, y] = GetSlice(this)
            % 
            
            % Slice the window of interest
            tWin = this.mapAxes.XLim;
            yWin = this.mapAxes.YLim;
            tb = this.se.SliceTimeSeries('LFP', tWin);
            
            % Filter across channels
            CData = tb.(2){1}';
            CData = medfilt1(CData, 5, [], 1);
            
            x = tb.time{1};
            y = this.channelTable.ycoords;
        end
    end
    
    methods(Static)
        function Video(ax, se)
            
            if isfield(se.userData, 'cameraInfo')
                atb = se.userData.cameraInfo.alignTb;
            else
                return;
            end
            
            tr = ax.UserData.trialNum;
            t = ax.UserData.time;
            if tr > height(atb) || tr < 1
                return;
            end
            
            % Initialize image object (must be done before initializing overlaying objects)
            if ~isfield(ax.UserData, 'img') || isempty(ax.UserData.img) || ~ishandle(ax.UserData.img)
                ax.UserData.img = imagesc(ax, zeros(1024, 1280));
                ax.CLim = [0 255];
                colormap(ax, 'gray');
                axis(ax, 'image');
                ax.XAxis.Visible = 'off';
                ax.YAxis.Visible = 'off';
                ax.TickLength(1) = 0;
                hold(ax, 'on');
            end
            
            % Find video file
            vidPath = [atb.vidNames{tr} '.avi'];
            if ~exist(vidPath, 'file')
                return;
            end
            
            if ~isfield(ax.UserData, 'vidObj')
                % Load the first video file
                ax.UserData.vidObj = VideoReader(vidPath);
                title(ax, ax.UserData.vidObj.Name, 'Interpreter', 'none');
                ax.TitleFontSizeMultiplier = 1.2;
            end
            vidObj = ax.UserData.vidObj;
            if ~strcmp(vidObj.Name, vidPath)
                % Load a new video file
                delete(vidObj);
                vidObj = VideoReader(vidPath);
                ax.UserData.vidObj = vidObj;
                title(ax, vidObj.Name, 'Interpreter', 'none');
                ax.TitleFontSizeMultiplier = 1.2;
            end
            
            % Find video time
            t = interp1(atb.tIntan{tr}, atb.tLed{tr}, t, 'linear', 'extrap');
            
            % Read video frame
            if t < 0 || t > vidObj.Duration
                return;
            end
            vidObj.CurrentTime = t;
            fr = rgb2gray(readFrame(vidObj));
            
            % Update plot
            ax.UserData.img.CData = fr;
        end
        
        function State(ax, se)
            
            if ~ismember('state', se.tableNames)
                return;
            end
            t0 = ax.UserData.time;
            z = 2^(-ax.UserData.zoom);
            tWin = t0 + 4 * [-z z];
            
            % Slice signals
            tb = se.SliceTimeSeries('state', tWin, 1, 'Fill', 'bleed');
            t = tb.time{1};
            s = tb.state{1};
            w = [-1 1];
            cc = lines(6);
            
            % Update plot
            cla(ax);
            hold(ax, 'on');
            plotRibbon(ax, s, w, cc, t, 1:6);
            ax.XLim = tWin;
            ax.YLim = w;
            ax.XLabel.String = 'Time (s)';
            ax.TickDir = 'out';
            ax.YAxis.Visible = 'off';
            
            plot(ax, [t0 t0]', ax.YLim', 'LineWidth', 2, 'Color', [0 0 0 .2])
            
            function plotRibbon(ax, xRange, yRange, colors, indVal, groups)
                
                % Standardize ranges
                xRange = findRange(xRange, groups, indVal);
                yRange = repmat({yRange(:)'}, size(xRange));
                
                function r = findRange(r, g, iv)
                    % Find ranges indices from a vector of discrete values
                    r = MMath.ValueBounds(r, g, 'Uni', false);
                    for k = 1 : numel(r)
                        r{k} = r{k} + [-.5 .5];
                    end
                    if ~isempty(iv)
                        % Convert indices to axis values
                        iv = iv(:);
                        ind = (1:numel(iv))';
                        for k = 1 : numel(r)
                            r{k} = interp1(ind, iv, r{k}, 'linear', 'extrap');
                        end
                    end
                end
                
                % Plot the ribbon
                for i = 1 : numel(xRange)
                    MPlot.Blocks(xRange{i}, yRange{i}, colors(i,:), 'Parent', ax);
                end
            end
        end
        
        function LFP(ax, se)
            if ~ismember('LFP', se.tableNames)
                return;
            end
            t0 = ax.UserData.time;
            z = 2^(-ax.UserData.zoom);
            tWin = t0 + 4 * [-z z];
            
            % Slice signals
            tb = se.SliceTimeSeries('LFP', tWin, 1, 'Fill', 'bleed');
            t = tb.time{1};
            v = tb.series1{1};
            
            % Update plot
            cla(ax);
            plot(ax, t, v, 'k');
            ax.XLim = tWin;
            ax.YLim = [-1 1] * 2e3;
            ax.XLabel.String = 'Time (s)';
            ax.TickLength = [0 0];
            ax.Box = 'off';
            grid(ax, 'on');
            
            hold(ax, 'on');
            plot(ax, [t0 t0]', ax.YLim', 'LineWidth', 2, 'Color', [0 0 0 .2])
        end
        
        function Spectrum1D(ax, se)
            
            if ~ismember('LFP', se.tableNames)
                return;
            end
            t0 = ax.UserData.time;
            z = 2^(-ax.UserData.zoom);
            tWin = t0 + 4 * [-z z];
            
            % Slice signals
            tb = se.SliceTimeSeries('LFP', tWin, 1, 'Fill', 'bleed');
            t = tb.time{1};
            v = tb.series1{1};
            
            % Transform signal
            fLims = [0 30];
            [P, F] = pspectrum(v, t, 'FrequencyLimits', fLims);
            
            % Plot
            plot(ax, F, P, 'k');
            ax.XLabel.String = 'Frequency (Hz)';
            ax.YLabel.String = 'Power';
            ax.TickDir = 'out';
        end
        
        function Spectrum2D(ax, se)
            
            if ~ismember('LFP', se.tableNames)
                return;
            end
            t0 = ax.UserData.time;
            z = 2^(-ax.UserData.zoom);
            tWin = t0 + 4 * [-z z];
            
            % Slice signals
            tb = se.SliceTimeSeries('LFP', tWin, 1, 'Fill', 'bleed');
            t = tb.time{1};
            v = tb.series1{1};
            
            % Transform signal
            fLims = [0 30];
            [P, F, T] = pspectrum(v, t, ...
                'spectrogram', ...
                'FrequencyLimits', fLims, ...
                'TimeResolution', 1.5, ...
                'OverlapPercent', 50);
            P = P.^(1/2);
            
            % Plot
            imagesc(ax, T, F, P);
            colormap(ax, 'bone');
            ax.CLim = [0 100];
            ax.XLim = tWin;
            ax.YLim = fLims;
            ax.XLabel.String = 'Time (s)';
            ax.YLabel.String = 'Frequency (Hz)';
            ax.YDir = 'normal';
            ax.TickDir = 'out';
        end
        
        function SpikeRaster(ax, se)
            
            k = ax.UserData.trialNum;
            if k > se.numEpochs || k < 1
                return;
            end
            if ~ismember('spikeTime', se.tableNames)
                return;
            end
            
            % Process spike times
            tb = se.SliceEventTimes('spikeTime', ax.UserData.timeLimits, k, 'Fill', 'bleed');
            [spkTimes, spkY, yTick] = SL.BehavFig.ConvertEventTimesForRasters(tb{1,:});
            
            % Update plot
            cla(ax);
            MPlot.PlotPointAsLine(spkTimes, spkY, .6, 'Color', [0 0 0], 'Parent', ax);
            hold(ax, 'on');
            
            ax.YTick = yTick;
            ax.XLabel.String = 'Time (s)';
            ax.YLabel.String = 'Unit';
            ax.TickLength = [0 0];
            ax.YGrid = 'on';
            ax.YLim = [0 width(tb)+1];
            
            % Callback
            ax.ButtonDownFcn = {@SL.MP.UnitPlots, se};
        end
        
        function isHit = ChangeState(d, se)
            
            % Process callback data
            modName = d.eventdata.Modifier;
            keyName = d.eventdata.Key;
            t0 = d.time;
            switch keyName
                case {'1', '2', '3', '4', '5', '6'}
                    val = str2double(keyName);
                case 'backquote'
                    val = 0;
                otherwise
                    isHit = false;
                    return
            end
            
            % Locate episode
            rt = se.GetReferenceTime();
            ep = find(t0 > rt, 1, 'last');
            stb = se.GetTable('state');
            t = stb.time{ep} + rt(ep);
            s = stb.state{ep};
            
            % Find range
            ids = find(diff([0; s]));
            [~, i0] = min(abs(t0 - t));
            ids1 = find(ids <= i0, 1, 'last');
            ids2 = find(ids >= i0, 1);
            
            if isempty(ids)
                % No transition labeled
                i1 = i0;
                i2 = i0;
            elseif ~isempty(ids1)
                % Has transition before
                i1 = ids(ids1);
                i2 = i0;
            else
                % Only has transition after
                i1 = i0;
                i2 = ids(ids2);
            end
            
            % Set state
            s(i1:i2) = val;
            stb.state{ep} = s;
            se.SetTable('state', stb);
            
            isHit = true;
        end
    end
    
end

