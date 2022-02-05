classdef MTracerVM < handle
    % 
    
    properties
        % UI
        app;                        % handle to app object
        appData;
        mapFig;                     % handle to Maps window
        mapAxes;                    % handle to the axes in Maps window
        mapLayers struct = struct;  % stores plot elements in mapAxes, each field is a layer
        
        % Data
        chanMapFile = 'NP1_NHP_HalfCol_kilosortChanMap.mat'; % channel map .mat file
        channelTable table;         % channel info
        apBinFile = '';             % path of ap.bin file
        imec Neuropixel.ImecDataset; % object used to access AP data
        lfBinFile = '';             % path of lf.bin file
        lfMeta struct;              % LFP metadata loaded from lf.meta
        ksFolder = '';              % path of Kilosort/Phy output folder
        rezOps struct;              % rez.ops struct
        ks Neuropixel.KilosortDataset; % object used to access sorting results
        ksm Neuropixel.KilosortMetrics; % object used to access sorting results
        kr NP.KilosortResult;
        
        recId char = 'NP0_B0';
        se MSessionExplorer;        % current data container, referenced to either seo or sek
        traces MTracerTrace;        % trace objects
        tracers struct;             % a list of auto tracers
                                    % field names are tracer names
                                    % each field saves the tracer data
        F1File = '';                % path of the saved motion scaling object
        F1;                         % movement scaling object for motion extrapolation
        F2File = '';                % path of the spatial-temporal interpolant
        F2 NP.MotionInterpolant;    % spatial-temporal interpolant for motion correction
        clustTb table;              % stores cluster data
        
        % Runtime variables
        focus = [0 0];
        isTracing = false;
        currentTrace = NaN;
        spikeMarkerSize = 4;
        currentCluster = [];
    end
    
    properties(Dependent)
        hasRez;
        hasChanMap;
        hasLFP;
        hasAP;
        hasTrace;
        hasInterp2;
        hasPhy;
        hasApp;
        hasMapAxes;
        layerNames;
        tLims;
        yLims;
    end
    
    methods
        function val = get.hasRez(this)
            val = ismember('spikes', this.se.tableNames);
        end
        function val = get.hasChanMap(this)
            val = ~isempty(this.channelTable);
        end
        function val = get.hasLFP(this)
            val = ismember('LFP', this.se.tableNames);
        end
        function val = get.hasAP(this)
            val = ismember('AP', this.se.tableNames);
        end
        function val = get.hasPhy(this)
            val = ~isempty(this.clustTb);
        end
        function val = get.hasTrace(this)
            val = ~isempty(this.traces);
        end
        function val = get.hasInterp2(this)
            val = ~isempty(this.F2);
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
            if this.hasRez
                val = [this.rezOps.tstart this.rezOps.tend] / this.rezOps.fs;
            elseif this.hasAP
                val = [0 this.imec.nSamplesAP] / 30e3;
            elseif this.hasLFP
                val = [0 str2double(this.lfMeta.fileTimeSecs)];
            elseif ~isempty(this.kr) && ~isempty(this.kr.mdat)
                val = [0 size(this.kr.mdat.Data.V, 2)-1] / 30e3;
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
            this.mapLayers.AP = [];
            this.mapLayers.anchors = [];
            this.mapLayers.interp = [];
            this.mapLayers.clusters = [];
            
            this.se = MSessionExplorer();
        end
        
        function obj = Duplicate(this)
            % Make a hard copy of the current MTracerVM
            
            pn = { ...
                'chanMapFile', 'channelTable', ...
                'apBinFile', 'imec', ...
                'lfBinFile', 'lfMeta', ...
                'ksFolder', 'rezOps', 'ks', 'ksm', 'kr', ...
                'recId', 'se', 'tracers', 'F1File', 'F1', 'F2File', 'F2', 'clustTb', ...
                'focus', 'currentTrace', 'appData'};
            
            obj = MTracerVM();
            
            for i = 1 : numel(pn)
                obj.(pn{i}) = this.(pn{i});
            end
            
            if this.hasTrace
                obj.traces = arrayfun(@(x) x.Duplicate(obj), this.traces);
            end
        end
        
        function delete(this)
            delete(this.mapFig);
        end
        
        % Data
        function LoadChannelMap(this, filePath)
            % Load the 354-channel configuration
            
            if nargin < 2 || isempty(filePath)
                filePath = MBrowse.File([], 'Please select the channel map file', '*.mat');
            end
            if ~exist(filePath, 'file')
                return
            end
            
            try
                s = load(filePath);
                
                tb = table;
                tb.ind = (1 : numel(s.chanMap))';
                tb.chanMap = s.chanMap;
                tb.xcoords = s.xcoords;
                tb.ycoords = s.ycoords;
                
                % Sort channels in ascending depth
                tb = sortrows(tb, 'ycoords');
                
                this.chanMapFile = filePath;
                this.channelTable = tb;
                
            catch e
                assignin('base', 'e', e);
                disp(e);
            end
        end
        
        function LoadAP(this, filePath)
            % Load and AP data to Neuropixel.util.imec
            
            if nargin < 2 || isempty(filePath)
                filePath = MBrowse.File([], 'Please select an ap.bin file', '*.bin');
            end
            if isempty(filePath)
                return
            end
            assert(this.hasChanMap, 'Channel Map must be loaded before loading ap.bin')
            
            try
                this.imec = Neuropixel.ImecDataset(filePath, 'channelMap', this.chanMapFile);
%                 this.mmap = this.imec.memmapAP_full();
                this.apBinFile = filePath;
                this.ExtractRecId(filePath);
                this.ReadApSlice();
            catch e
                assignin('base', 'e', e);
                disp(e);
            end
        end
        
        function ReadApSlice(this)
            % 
            
            if isempty(this.imec)
                return
            end
            
            % Limit time window size
            if this.hasMapAxes
                tWin = this.mapAxes.XLim;
            else
                tWin = [0 Inf];
            end
            maxSpan = 4;
            if diff(tWin) > maxSpan
                tWin = this.focus(1) + [-1 1]*maxSpan/2;
            end
            tWin = MMath.Bound(tWin, this.tLims);
            
            % Read data
            [v, ind] = this.imec.readAP_timeWindow(tWin);
            v = v(this.channelTable.ind, :)'; % reorder by depth
            t = ind' / 30e3; % consistent with rez.mat
            
            % Set to table
            tb = table;
            tb.time = {t};
            tb.v = {v};
            this.se.SetTable('AP', tb, 'timeSeries');
        end
        
        function LoadLFP(this, filePath)
            % Load and preprocess LFP data
            
            if nargin < 2 || isempty(filePath)
                filePath = MBrowse.File([], 'Please select a lf.bin file', '*.bin');
            end
            if ~exist(filePath, 'file')
                return
            end
            assert(this.hasChanMap, 'Channel Map must be loaded before loading lf.bin')
            
            try
                % Load LFP data
                disp('Loading lf.bin ...');
                lfMetaFile = strrep(filePath, '.bin', '.meta');
                [meta, v, t] = MSpikeGLX.ReadLFP(lfMetaFile);
                v = v(:, this.channelTable.ind); % sort channels by depth
                
                % Make downsampled LFP timeSeries table
                disp('Downsampling LFP to 50Hz ...');
                lfpTb = table();
                r = 50; % downsample 50x from 2500Hz to 50Hz
                lfpTb.time{1} = downsample(t, r, r/2);
                lfpTb.v{1} = MMath.Decimate(v, r, r/2);
                
                this.lfBinFile = filePath;
                this.lfMeta = meta;
                this.se.SetTable('LFP', lfpTb, 'timeSeries', 0);
                disp('LFP data loaded');
                
                % Extract recording ID
                this.ExtractRecId(filePath);
                
            catch e
                assignin('base', 'e', e);
                disp(e);
            end
        end
        
        function LoadKilosortRez(this, filePath)
            % Load spiking data from Kilosort output rez.mat
            
            if nargin < 2 || isempty(filePath)
                filePath = MBrowse.File([], 'Please select a rez.mat file', '*.mat');
            end
            if ~exist(filePath, 'file')
                return
            end
            
            try
                % Load rez.mat
                disp('Loading rez.mat ...');
                load(filePath, 'rez');
                
                % Make a table of spikes
                tb = table();
                tb.time = (rez.st0(:,1) + rez.ops.tstart) / rez.ops.fs;
                tb.y = rez.st0(:,2);
                tb.amp = rez.st0(:,3);
                tb = sortrows(tb, {'time', 'y'});
                C = mat2cell(tb{:,:}, height(tb), ones(1,width(tb)));
                tb = cell2table(C, 'VariableNames', tb.Properties.VariableNames);
                
                this.ksFolder = fileparts(filePath);
                this.rezOps = rez.ops;
                this.se.SetTable('spikes', tb, 'timeSeries', 0);
                disp('rez.mat data loaded');
                
                % Extract recording ID
                this.ExtractRecId(filePath);
                
            catch e
                assignin('base', 'e', e);
                disp(e);
            end
        end
        
        function LoadPhy(this, ksDir)
            % 
            
            if nargin < 2 || isempty(ksDir)
                ksDir = MBrowse.Folder([], 'Please select the Kilosort output folder');
            end
            if ~exist(ksDir, 'dir')
                return
            end
            assert(this.hasChanMap, 'Channel Map must be loaded before loading ap.bin')
            
            try
                % Construct NP.KilosortResult
                if this.hasRez
                    this.kr = NP.KilosortResult(ksDir, this.chanMapFile, this.rezOps.tstart);
                else
                    warning('If sorting was performed on temporally truncated data, please load rez.mat first for the sample offset value.');
                    this.kr = NP.KilosortResult(ksDir, this.chanMapFile);
                end
                this.kr.ComputeAll();
                
                % Gather spike data
                spkTb = table;
                spkTb.tSpk = double(this.kr.spkTb.timeInd) / 30e3;
                spkTb.ySpk = this.kr.spkTb.covCentCoords(:,2);
                
                % Group spikes by clusters
                c = this.kr.spkTb.clusIds;
                [G, ID] = findgroups(c);
                cTb = table;
                for i = 1 : width(spkTb)
                    cTb.(i) = splitapply(@(x) {x}, spkTb.(i), G);
                end
                cTb.Properties.VariableNames = spkTb.Properties.VariableNames;
                cTb.id = ID;
                cTb = [cTb this.kr.clusTb(:, {'group', 'SNR', 'RPV', 'contam'})];
                cTb.n_spikes = cellfun(@numel, cTb.tSpk);
                cTb.depth = round(cellfun(@(x) median(x,'omitnan'), cTb.ySpk));
                cTb(strcmp(cTb.group,'noise'), :) = [];
                
                % Sort rows by depth
                cTb = sortrows(cTb, 'depth', 'descend');
                
                % Columns about plotting
                cTb.color = lines(height(cTb));
                cTb.handle{1} = [];
                
                this.ksFolder = ksDir;
                this.clustTb = cTb;
                this.currentCluster = [];
                this.PlotClusters();
                
                % Extract recording ID
                this.ExtractRecId(ksDir);
                
            catch e
                assignin('base', 'e', e);
                disp(e);
            end
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
                    
                    % Extract recording ID
                    this.ExtractRecId(filePaths{i});
                catch e
                    warning('Error occured when loading: \n%s', filePaths{i});
                    disp(e);
                end
            end
            
            this.PlotTraces();
        end
        
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
            
            disp('Tracing results are saved to the current directory.');
            if this.hasApp
                uialert(this.app.UIFigure, 'Tracing results are saved to the current directory.', 'Success', 'Icon', 'Info');
            end
        end
        
        function s = SliceMap(this, varargin)
            % Get map data in the current view
            % 
            %   s = obj.SliceMap()
            %   s = obj.SliceMap(tWin)
            %   s = obj.SliceMap(tWin, yWin)
            %   s = obj.SliceMap(..., 'DataNames', {'traces', 'spikes', 'LFP', 'AP'})
            %   s = obj.SliceMap('TraceType', 'interp')
            % 
            % Inputs
            %   tWin            tbw
            %   yWin            tbw
            %   'DataNames'     tbw
            %   'TraceType'     tbw
            % Output
            %   s               tbw
            
            p = inputParser;
            p.addOptional('tWin', this.mapAxes.XLim, @(x) numel(x)==2 && isnumeric(x));
            p.addOptional('yWin', this.mapAxes.YLim, @(x) numel(x)==2 && isnumeric(x));
            p.addParameter('DataNames', {'traces', 'spikes', 'LFP', 'AP'});
            p.addParameter('TraceType', 'interp', @(x) ismember(x, {'anchor', 'interp'}));
            p.parse(varargin{:});
            tWin = p.Results.tWin;
            yWin = p.Results.yWin;
            dn = cellstr(p.Results.DataNames);
            traceType = p.Results.TraceType;
            
            s.tWin = tWin;
            s.yWin = yWin;
            s.focus = this.focus;
            
            if ~isempty(this.channelTable)
                tb = this.channelTable;
                chanInd = tb.ycoords >= yWin(1) & tb.ycoords <= yWin(2);
                tb = tb(chanInd,:);
                s.channelTable = tb;
            end
            
            if ismember('spikes', dn) && this.hasRez
                tb = this.se.SliceTimeSeries('spikes', tWin);
                tb = DenestTable(tb);
                spInd = tb.y >= yWin(1) & tb.y <= yWin(2);
                s.spikes = tb(spInd,:);
            end
            
            if ismember('LFP', dn) && this.hasLFP
                tb = this.se.SliceTimeSeries('LFP', tWin);
                tb = DenestTable(tb);
                tb.v = tb.v(:,chanInd);
                s.LFP = tb;
            end
            
            if ismember('AP', dn) && this.hasAP
                tb = this.se.SliceTimeSeries('AP', tWin);
                tb = DenestTable(tb);
                tb.v = tb.v(:,chanInd);
                s.AP = tb;
            end
            
            if ismember('traces', dn) && this.hasTrace
                tb = table;
                
                for k = numel(this.traces) : -1 : 1
                    [t, y, ti, yi] = this.traces(k).GetCoords();
                    if strcmp(traceType, 'interp')
                        t = ti;
                        y = yi;
                    end
                    ptInd = t >= tWin(1) & t <= tWin(2) & y >= yWin(1) & y <= yWin(2);
                    if sum(ptInd) > 1
                        tb.time{k} = t(ptInd);
                        tb.y{k} = y(ptInd);
                    end
                end
                s.traces = tb(~cellfun(@isempty, tb.time), :);
            end
            
            function tbOut = DenestTable(tbIn)
                tbOut = table;
                for i = 1 : width(tbIn)
                    n = tbIn.Properties.VariableNames{i};
                    tbOut.(n) = tbIn.(n){1};
                end
            end
        end
        
        function ApplyCorrection(this)
            % 
            
            if ~this.hasInterp2
                if this.hasApp
                    uialert(this.app.UIFigure, 'The interpolant is not loaded', 'Oops...', 'Icon', 'warning');
                end
                error('A spatial-temporal interpolant must be set to the property F2 before applying correction.');
            end
            
            if this.hasRez
                tb = this.se.GetTable('spikes');
                tb.y = cellfun(@(t,y) this.F2.CorrectDepths(t,y), tb.time, tb.y, 'Uni', false);
                this.se.SetTable('spikes', tb);
            end
            
            if this.hasLFP
                tb = this.se.GetTable('LFP');
                y = this.channelTable.ycoords';
                tb.v = cellfun(@(t,v) this.F2.CorrectVoltArray(t,y,v), tb.time, tb.v, 'Uni', false);
                this.se.SetTable('LFP', tb);
            end
            
            if this.hasAP
                tb = this.se.GetTable('AP');
                y = this.channelTable.ycoords';
                tb.v = cellfun(@(t,v) this.F2.CorrectVoltArray(t,y,v,'Highpass',true), tb.time, tb.v, 'Uni', false);
                this.se.SetTable('AP', tb);
            end
            
            if this.hasTrace
                this.traces = arrayfun(@(x) this.F2.CorrectTraces(x), this.traces);
            end
            
            if this.hasMapAxes
                delete(this.mapAxes);
                this.PlotAll();
            end
        end
        
        function ExtractRecId(this, s)
            % Extract recording ID from a string, typically a path
            
            % Don't overwrite existing ID
            if ~strcmp(this.recId, 'NP0_B0') || isempty(s)
                return
            end
            
            % Extract ID
            matches = regexp(s, 'NP\d{2,}_B\d+', 'match');
            if ~isempty(matches)
                this.recId = matches{1};
                if this.hasApp
                    this.app.RecEditField.Value = this.recId;
                end
            end
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
        
        function PlotAll(this)
            this.InitializeMapAxes();
            this.PlotLFP();
            this.PlotSpikes();
            this.PlotClusters();
            this.PlotTraces();
            this.PlotFocus();
        end
        
        function PlotFocus(this)
            % Plot the focus indicator
            
            if this.hasApp
                this.app.FocusTimeEditField.Value = this.focus(1);
                this.app.FocusDepthEditField.Value = this.focus(2);
            end
            
            if ~this.hasMapAxes
                return
            end
            
            xx = [this.focus([1 1])' this.tLims'];
            yy = [this.yLims' this.focus([2 2])'];
            
            if this.isTracing
                cc = [0 1 0 .5];
            else
                cc = [0 1 1 .7];
            end
            
            h = this.mapLayers.focus;
            if isempty(h) || any(~isvalid(h))
                h = plot(xx, yy, 'Color', cc, 'LineWidth', 1, 'HitTest', 'off');
                delete(this.mapLayers.focus);
                this.AddHandle2Layer(h, 'focus');
            else
                set(h(1), 'XData', this.focus([1 1])', 'YData', this.yLims', 'Color', cc);
                set(h(2), 'XData', this.tLims', 'YData', this.focus([2 2])', 'Color', cc);
            end
            uistack(h, 'top', Inf);
        end
        
        function PlotSpikes(this)
            % Plot spike events as a function of time and depth
            
            if ~this.hasRez || ~this.hasMapAxes
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
                hh{k} = plot(this.mapAxes, t(ind), y(ind), '.', ...
                    'MarkerSize', this.spikeMarkerSize, ...
                    'Color', [1 1 1] * max(0, 1-ampRange(k)/40), ... % the marker color here has been carefully tuned
                    'HitTest', 'off');
            end
            
            MPlot.Axes(this.mapAxes);
            this.UpdateMapLims(this.tLims, this.yLims);
            
            delete(this.mapLayers.spikes);
            this.AddHandle2Layer(cat(1, hh{:}), 'spikes');
        end
        
        function PlotClusters(this)
            % Plot all or selected clusters
            
            if ~this.hasPhy || ~this.hasMapAxes
                return
            end
            
            tb = this.clustTb;
            
            if isempty(this.currentCluster)
                ids = tb.id;
            else
                ids = this.currentCluster;
            end
            
            % Plot the selected cluster later, thus on top of unselected ones
            ind = ismember(tb.id, ids);
            tb = [tb(~ind,:); tb(ind,:)];
            
            % Clear existing clusters in map
            delete(this.mapLayers.clusters);
            
            % Plot clusters
            for i = 1 : height(tb)
                t = tb.tSpk{i};
                y = tb.ySpk{i};
                c = tb.color(i,:);
                k = tb.id(i);
                g = tb.group(i);
                
                % Add alpha to unselected clusters, with "mua" lighter than "good"
                if g == "good"
                    r = 0.4;
                else % "mua"
                    r = 0.2;
                end
                if ~ismember(k, ids)
                    w = ones(size(c));
                    c = r*c + (1-r)*w;
                end
                
                % Plot spikes
                h = plot(this.mapAxes, t, y, '.', 'Color', c, 'MarkerSize', this.spikeMarkerSize, 'HitTest', 'off');
                
                % Label selected plots
                if ismember(k, ids)
                    h.UserData.id = k;
                    h.UserData.group = g;
                end
                
                tb.handle{i} = h;
            end
            this.AddHandle2Layer(cat(1, tb.handle{:}), 'clusters');
            this.UpdateSpikeMarkerSize(1, 'relative');
            
            % Initialize brush callback
            brushObj = brush(this.mapFig);
            brushObj.ActionPostCallback = @this.PlotBrushedSpikes;
        end
        
        function PlotBrushedSpikes(this, fig, axesStruct)
            % 
            
            h = axesStruct.Axes.Children;
            cid = NaN(numel(h), 1);
            spkMask = cell(size(cid));
            for i = 1 : numel(h)
                % Check if the handle is cluster spike plot
                s = h(i).UserData;
                if ~all(isfield(s, {'id', 'group'}))
                    continue
                end
                
                % Check if selecting any data
                b = logical(h(i).BrushData);
                if ~any(b)
                    continue
                end
                
                cid(i) = s.id;
                spkMask{i} = b;
            end
            isBrushed = ~isnan(cid);
            if ~any(isBrushed)
                return
            end
            cid = cid(isBrushed);
            spkMask = spkMask(isBrushed);
            
            % Get cluster colors
            for i = numel(cid) : -1 : 1
                m = cid(i)==this.clustTb.id;
                cc(i,:) = this.clustTb.color(m,:);
            end
%             if isscalar(cid)
%                 cc = [0 0 0];
%             end
            
            % Make plot
            f = MPlot.Figure(301); clf
            f.Name = 'Spike Waveform';
            this.kr.PlotClusterWaveform(cid, 50, spkMask, 'NumChannels', 32, 'Color', cc, 'ShowCentroids', true);
            ax = MPlot.Axes(gca);
            ax.Title.String = ['Cluster ' num2str(cid(:)')];
        end
        
        function PlotCCG(this, cid)
            % 
            
            if ~this.hasPhy
                return
            end
            
            if nargin < 2
                cid = this.currentCluster;
            end
            
            % No plot when selecting to many clusters
            if numel(cid) > 5
                return
            end
            
            tb = this.clustTb;
            tb = tb(ismember(tb.id, cid), :);
            
            f = MPlot.Figure(200); clf
            f.Name = 'CCG';
            this.kr.PlotCCG(tb.id, 'Color', tb.color, 'ShowMetrics', true);
            
%             tb = this.clustTb;
%             tb = tb(ismember(tb.id, cid), :);
%             
%             f = MPlot.Figure(200); clf
%             f.Name = 'CCG';
%             tEdge = -0.02 : 0.5e-3 : 0.02;
%             ccg = MNeuro.CCG(tEdge, tb.tSpk{:});
%             N = height(tb);
%             tBin = tEdge(2:end) - diff(tEdge)/2;
%             for i = 1 : N
%                 for j = i : N
%                     ax = subplot(N, N, (i-1)*N+j);
%                     h = bar(tBin*1e3, squeeze(ccg(i,j,:)), 'histc');
%                     h.EdgeColor = 'none';
%                     h.FaceColor = (tb.color(i,:)+tb.color(j,:)) / 2;
%                     MPlot.Axes(ax);
%                     ax.XLim = tEdge([1 end])*1e3;
%                 end
%             end
        end
        
        function PlotWaveform(this, cid)
            % Plot 
            
            if ~this.hasPhy
                return
            end
            
            if nargin < 2
                cid = this.currentCluster;
            end
            
            % No plot when selecting to many clusters
            if numel(cid) > 5
                return
            end
            
            % Find cluster colors
            for i = numel(cid) : -1 : 1
                m = cid(i)==this.clustTb.id;
                cc(i,:) = this.clustTb.color(m,:);
            end
%             if isscalar(ids)
%                 cc = [0 0 0];
%             end
            
            % Plot waveforms
            f = MPlot.Figure(300); clf
            f.Name = 'Spike Waveform';
            this.kr.PlotClusterWaveform(cid, 50, 'NumChannels', 10, 'Color', cc);
            ax = MPlot.Axes(gca);
            ax.Title.String = ['Cluster ' num2str(cid(:)')];
%             ax.LooseInset = [0 0 0 0];
        end
        
        function PlotLFP(this)
            % Plot LFP as a heatmap
            
            if ~this.hasLFP || ~this.hasMapAxes
                return
            end
            
            tb = this.se.GetTable('LFP');
            C = tb.v{1}';
            x = tb.time{1};
            y = this.channelTable.ycoords;
            
            h = imagesc(this.mapAxes, x, y, C, 'HitTest', 'off');
            colormap(this.mapAxes, MPlot.PolarMap);
            MPlot.Axes(this.mapAxes);
            this.mapAxes.YDir = 'normal';
            this.UpdateMapLims(this.tLims, this.yLims);
            
            delete(this.mapLayers.LFP)
            this.AddHandle2Layer(h, 'LFP');
        end
        
        function PlotAP(this)
            % Plot AP signal
            
            if ~this.hasMapAxes
                return
            end
            
            % Read a temporal slice of AP data from file
            this.ReadApSlice();
            if ~this.hasAP
                return
            end
            tb = this.se.GetTable('AP');
            V = tb.v{1};
            t = tb.time{1};
            y = this.channelTable.ycoords;
            
            % Limit y window to 2000um maximal
            yWin = this.mapAxes.YLim;
            if diff(yWin) > 2000
                yWin = this.focus(2) + [-1 1]*1000;
            end
            
            % Choose channel interval based on y window size
            dy = y(2) - y(1);
            itvl = MMath.Bound(diff(yWin)/dy/50, [1 2 3 4]);
            isChan = y > yWin(1) & y <= yWin(2);
            chanInd = find(isChan, 1, 'first') : itvl : find(isChan, 1, 'last');
            
            % Use every other sample in time
            tInd = 1 : 2 : numel(t);
            t = t(tInd);
            V = V(tInd, chanInd);
            y = y(chanInd);
            
            % Highpass filtering
            isHp = true;
            if isHp
                % Prepare parameters for high-pass filtering
                %   Phy uses bandpass between 500Hz and 14.25kHz(.475*sample_rate) to visulize waveform;
                %   here uses 300Hz highpass to be consistent with Kilosort configuration
                D = designfilt('highpassiir', ...
                    'PassbandFrequency', 300, ...
                    'StopbandFrequency', 250, ...
                    'StopbandAttenuation', 60, ...
                    'PassbandRipple', 0.1, ...
                    'SampleRate', 1/diff(t([1 2])), ...
                    'DesignMethod', 'ellip'); % order of this filter is 8
                
                V = double(V);
                for i = 1 : size(V,2)
                    V(:,i) = filtfilt(D, V(:,i));
                end
            end
            
            % Scale voltage timeseries
            V = zscore(V) * 3;
            
            % Plot
            hh = line(t, V + y', 'Color', [0 0 0 .3], 'Parent', this.mapAxes, 'HitTest', 'off');
            delete(this.mapLayers.AP)
            this.AddHandle2Layer(hh, 'AP');
            
%             tWin = this.mapAxes.XLim;
        end
        
        function PlotTraces(this)
            % Plot traces
            if this.hasApp
                this.app.TraceListBox.Items = arrayfun(@(x) x.dispName, this.traces, 'Uni', false);
            end
            if this.hasMapAxes
                for i = 1 : numel(this.traces)
                    this.traces(i).PlotTrace(ismember(i, this.currentTrace));
                end
            end
        end
        
        function ToggleLayer(this, num)
            % Toggle layer on and off
            if num < 1 || num > numel(this.layerNames)
                return
            end
            hh = this.mapLayers.(this.layerNames{num});
            for i = 1 : numel(hh)
                hh(i).Visible = ~hh(i).Visible;
            end
        end
        
        function AddHandle2Layer(this, h, layerName)
            % Add plot object handle(s) to a layer
            
            if ~isfield(this.mapLayers, layerName) || isempty(this.mapLayers.(layerName))
                % Add object to the new or empty layer
                this.mapLayers.(layerName) = h;
            else
                % Append object handle to the end
                this.mapLayers.(layerName)(end+1:end+numel(h), 1) = h;
                this.ClearInvalidHandles();
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
        
        function DeleteLayers(this, varargin)
            for i = 1 : numel(varargin)
                layerName = varargin{i};
                if isfield(this.mapLayers, layerName)
                    delete(this.mapLayers.(layerName));
                    this.mapLayers = rmfield(this.mapLayers, layerName);
                end
            end
            this.UpdateLayerTree();
        end
        
        function UpdateMapLims(this, tWin, yWin, t, y)
            
            if ~this.hasMapAxes
                return
            end
            
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
            
            if this.hasApp
                tWin = round(tWin, 3);
                yWin = round(yWin, 3);
                this.app.TimeWin1EditField.Value = tWin(1);
                this.app.TimeWin2EditField.Value = tWin(2);
                this.app.DepthWin1EditField.Value = yWin(1);
                this.app.DepthWin2EditField.Value = yWin(2);
            end
        end
        
        function UpdateSpikeMarkerSize(this, r, valType)
            
            if ~this.hasMapAxes
                return
            end
            
            % Change the current spike marker size
            sz = this.spikeMarkerSize;
            if strcmp(valType, 'relative')
                sz = sz * r;
            else
                sz = r;
            end
            sz = MMath.Bound(sz, [2 16]);
            this.spikeMarkerSize = sz;
            
            % Apply marker size to spikes
            if this.hasRez
                hh = this.mapLayers.spikes;
                for i = 1 : numel(hh)
                    hh(i).MarkerSize = sz;
                end
            end
            
            % Apply marker size to clusters
            if this.hasPhy
                hh = this.mapLayers.clusters;
                for i = 1 : numel(hh)
                    hh(i).MarkerSize = sz;
                end
            end
        end
        
        % Tracing
        function SelectTrace(this, t ,y)
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
        end
        
        function EnterTracingMode(this)
            if isscalar(this.currentTrace)
                this.isTracing = true;
                this.PlotFocus();
                disp('Tracing Mode');
            else
                error('Cannot enter tracing mode with more than one selected trace\n');
            end
        end
        
        function EnterSelectionMode(this)
            this.isTracing = false;
            this.PlotFocus();
            disp('Selection Mode');
        end
        
        function InitializeTrace(this)
            % Create a trace and set it as current
            this.traces(end+1,1) = MTracerTrace(this);
            this.currentTrace = numel(this.traces);
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
        
        % User Input Callbacks
        function KeyPress(this, src, eventdata)
            
            K = eventdata.Key;
            M = eventdata.Modifier;
            
            if ismember('control', M)
                switch K
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
                    case 'v'
                        this.PlotAP();
                        
                    case 't'
                        this.EnterTracingMode();
                    case 'escape'
                        this.EnterSelectionMode();
                        
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
                        this.UpdateSpikeMarkerSize(2, 'relative');
                    case 'comma'
                        % Decrease spike marker size
                        this.UpdateSpikeMarkerSize(0.5, 'relative');
                        
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
        
        function SetPoint(this, src, eventdata)
            % Place a point on image
            
            % Read mouse position
            mousePos = get(this.mapAxes, 'CurrentPoint');
            
            % Update reference point
            t = MMath.Bound(mousePos(1), this.tLims);
            y = MMath.Bound(mousePos(3), this.yLims);
            this.focus = [t y];
            this.PlotFocus();
            
            if this.isTracing
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
            else
                % Trace selection
                this.SelectTrace(t, y);
            end
        end
        
%         function PlotWaveformOld(this, ids)
%             
%             if ~this.hasPhy
%                 return
%             end
%             if numel(ids) > 10
%                 return
%             end
%             
%             f = MPlot.Figure(100); clf
%             f.Name = 'Spike Waveform';
%             snippetSet = this.ks.getWaveformsFromRawData('cluster_ids', ids, ...
%                 'num_waveforms', 50, 'best_n_channels', 10, 'car', false);
%             snippetSet.plotAtProbeLocations('alpha', 0.8);
%             ax = MPlot.Axes(gca);
%             ax.LooseInset = [0 0 0 0];
%         end
        
    end
    
end

