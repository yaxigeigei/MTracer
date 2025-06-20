classdef MTracerVM < MTracer.LayeredFigure
    % 
    
    properties(Constant)
        cacheFolderName = 'mtracer_cache';
    end
    
    properties
        % UI
        app                         % handle to app object
        appData                     % 
        
        % Data
        recId char = 'NP0_B0'
        chanMapFile = ''            % channel map .mat file
        channelTable table          % channel info
        apBinFile = ''              % path of ap.bin file
        imec Neuropixel.ImecDataset % object used to access AP data
        lfBinFile = ''              % path of lf.bin file
        lfMeta struct               % LFP metadata loaded from lf.meta
        se MSessionExplorer         % container of voltage and raw spike data, all stored in one epoch
        
        ksFolder = ''               % path of Kilosort/Phy output folder
        rezOps struct               % rez.ops struct
        clustering MTracer.ClusteringVM % view-model that manages clusters and clustering
        
        traces MTracerTrace         % trace objects
        tracers struct              % a list of auto tracers, field names are tracer names, each field saves the tracer data
        
        F1File = ''                 % path of the saved motion scaling object
        F1                          % movement scaling object for motion extrapolation
        F2File = ''                 % path of the spatial-temporal interpolant
        F2                          % spatial-temporal interpolant for motion correction
        
        % Runtime variables
        mapMode = 'selection'       % 'selection', 'tracing', or 'clustering'
        apSource = 'imec.ap.bin'    % must be 'imec.ap.bin' or 'temp_wh.dat'
        focus = [0 0]
        currentTrace = NaN          % currently selected trace(s)
        currentTracer = 'none'      % currently selected auto tracer
        spikeMarkerSize = 4
    end
    
    properties(Dependent)
        hasApp
        hasChanMap
        hasAP
        hasLFP
        hasRez
        hasClus
        hasTrace
        hasInterp2
        tLims
        yLims
    end
    methods
        function val = get.hasApp(this)
            val = ~isempty(this.app);
        end
        function val = get.hasChanMap(this)
            val = ~isempty(this.channelTable);
        end
        function val = get.hasAP(this)
            val = ~isempty(this.imec);
        end
        function val = get.hasLFP(this)
            val = ismember('LFP', this.se.tableNames);
        end
        function val = get.hasRez(this)
            val = ismember('spikes', this.se.tableNames);
        end
        function val = get.hasClus(this)
            val = this.clustering.hasClus;
        end
        function val = get.hasTrace(this)
            val = ~isempty(this.traces);
        end
        function val = get.hasInterp2(this)
            val = ~isempty(this.F2);
        end
        function val = get.tLims(this)
            if this.hasRez
                val = [this.rezOps.tstart this.rezOps.tend] / this.rezOps.fs;
            elseif this.hasAP
                val = [0 this.imec.nSamplesAP] / 30e3;
            elseif this.hasLFP
                val = [0 str2double(this.lfMeta.fileTimeSecs)];
            elseif this.clustering.hasClus
                ops = this.clustering.sr.rezLite.ops;
                val = ops.trange;
                val(2) = min(val(2), ops.tend / ops.fs);
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
            
            this.se = MSessionExplorer();
            this.clustering = MTracer.ClusteringVM(this);
            
            this.layers.focus = [];      % 1
            this.layers.spikes = [];     % 2
            this.layers.clusters = [];   % 3
            this.layers.AP = [];         % 4
            this.layers.LFP = [];        % 5
            this.layers.anchors = [];    % 6
            this.layers.interp = [];     % 7
        end
        
        function vm = Duplicate(this)
            % Make a hard copy of the current MTracerVM
            
            % Construct an empty object
            vm = MTracerVM();
            
            % Fields to copy directly
            pn = { ...
                'name', 'figPos', 'axPos', ... % superclass properties
                'appData', ...
                'recId', 'chanMapFile', 'channelTable', 'apBinFile', 'imec', 'lfBinFile', 'lfMeta', 'se', ...
                'ksFolder', 'rezOps', ...
                'tracers', ...
                'F1File', 'F1', 'F2File', 'F2', ...
                'apSource', 'focus', 'currentTrace', 'spikeMarkerSize'};
            for i = 1 : numel(pn)
                vm.(pn{i}) = this.(pn{i});
            end
            
            % Make a hard copy of the sorting result object
            vm.clustering = this.clustering.Duplicate(vm);
            
            % Make hard copies of the traces associated with the new main VM
            if this.hasTrace
                vm.traces = arrayfun(@(x) x.Duplicate(vm), this.traces);
            end
        end
        
        function delete(this)
            delete(this.hFig);
            delete(this.clustering);
        end
        
        % Data
        function LoadChannelMap(this, chanMapPath)
            % Load the 354-channel configuration
            
            if nargin < 2 || isempty(chanMapPath)
                chanMapPath = MBrowse.File([], 'Please select the channel map file', '*.mat');
            end
            if ~exist(chanMapPath, 'file')
                return
            end
            
            try
                chanTb = MKilosort2.LoadChanMap2Table(chanMapPath);
                [chanTb, I] = sortrows(chanTb, {'shankInd', 'ycoords', 'xcoords'}, {'ascend', 'descend', 'ascend'});
                chanTb.sortInd = I;
                
                this.chanMapFile = chanMapPath;
                this.channelTable = chanTb;
                
            catch e
                assignin('base', 'e', e);
                disp(e);
            end
        end
        
        function LinkAP(this, filePath)
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
            
            % Limit time window size
            if this.isOpen
                tWin = this.hAxes.XLim;
            else
                tWin = [0 Inf];
            end
            maxSpan = 2;
            if diff(tWin) > maxSpan
                tWin = this.focus(1) + [-1 1]*maxSpan/2;
            end
            tWin = MMath.Bound(tWin, this.tLims);
            
            % Read data
            if strcmp(this.apSource, 'imec.ap.bin') && this.hasAP
                [V, ind] = this.imec.readAP_timeWindow(tWin);
                V = V(this.channelTable.sortInd, :)'; % reorder by depth
                y = this.channelTable.ycoords;
            elseif strcmp(this.apSource, 'temp_wh.dat') && this.hasClus
                spWin = round(tWin * 30e3);
                V = this.clustering.sr.ReadSnippets(0, spWin)';
                ind = spWin(1) : spWin(2);
                y = this.clustering.sr.chanTb.ycoords;
            else
                return
            end
            V = double(V);
            t = ind' / 30e3; % consistent with rez.mat
            
            % Set to table
            tb = table;
            tb.time = {t};
            tb.v = {V};
            tb.Properties.UserData.y = y;
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
                v = v(:, this.channelTable.sortInd); % sort channels by depth
                
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
        
        function LoadRez(this, ksDir)
            % Load spiking data from Kilosort output rez.mat
            
            if nargin < 2 || isempty(ksDir)
                ksDir = MBrowse.Folder([], 'Please select the Kilosort output folder');
            end
            if ~exist(ksDir, 'dir')
                return
            end
            
            % Try loading from cache
            [s, info] = this.LoadCache(fullfile(ksDir, this.cacheFolderName, 'rez_lite_*.mat'));
            if strcmp(info, 'cancelled')
                return
            end
            
            try
                if ~isempty(s)
                    rezLite = s.rezLite;
                else
                    % Check source file
                    rezPath = fullfile(ksDir, 'rez.mat');
                    if ~exist(rezPath, 'file')
                        if this.hasApp
                            uialert(this.app.UIFigure, 'Cannot find rez.mat file in the KS output folder', 'Hmmm...', 'Icon', 'warning')
                        end
                        return
                    end
                    
                    % Load rez.mat
                    disp('Loading rez.mat ...');
                    s = load(rezPath, 'rez');
                    
                    % Cache a lightweight version of rez
                    if isfield(s.rez, 'st0')
                        rezLite.st0 = s.rez.st0;
                    else
                        fprintf("Raw drift map data (rez.st0) is not available. Construct one from sorted spikes (rez.st3).\n");
                        st = s.rez.st3;
                        t = st(:,1) - s.rez.ops.tstart;
                        iTemp = st(:,2);
                        iChan = s.rez.iNeighPC(16, iTemp);
                        y = s.rez.ycoords(iChan);
                        a = round(st(:,3));
                        rezLite.st0 = [t, y, a];
                    end
                    rezLite.ops = s.rez.ops;
                    cacheFolder = fullfile(ksDir, this.cacheFolderName);
                    if ~exist(cacheFolder, 'dir')
                        mkdir(cacheFolder);
                    end
                    save(fullfile(cacheFolder, ['rez_lite_' datestr(now, 'yyyy-mm-dd_HH-MM-SS') '.mat']), 'rezLite');
                end
                
                % Make a table of spikes
                tb = table();
                tb.time = (rezLite.st0(:,1) + rezLite.ops.tstart) / rezLite.ops.fs;
                tb.y = rezLite.st0(:,2);
                tb.amp = rezLite.st0(:,3);
                tb = sortrows(tb, {'time', 'y'});
                C = mat2cell(tb{:,:}, height(tb), ones(1,width(tb)));
                tb = cell2table(C, 'VariableNames', tb.Properties.VariableNames);
                
                this.ksFolder = ksDir;
                this.rezOps = rezLite.ops;
                this.se.SetTable('spikes', tb, 'timeSeries', 0);
                disp('rez.mat data loaded');
                
                % Update UI
                this.ExtractRecId(ksDir);
                this.PlotRezSpikes();
                
            catch e
                assignin('base', 'e', e);
                disp(e);
                if this.hasApp
                    uialert(this.app.UIFigure, 'Unable to load rez.mat', 'Hmmm...', 'Icon', 'warning')
                end
            end
        end
        
        function LoadClusResults(this, ksDir)
            % Load clustering results using MTracer.KilosortResult class
            
            if nargin < 2 || isempty(ksDir)
                ksDir = MBrowse.Folder([], 'Please select the Kilosort output folder');
            end
            if ~exist(ksDir, 'dir')
                return
            end
            this.ksFolder = ksDir;
            this.ExtractRecId(ksDir);
            
            % Try loading from cache
            [s, info] = this.LoadCache(fullfile(ksDir, this.cacheFolderName, 'sr_*.mat'));
            if strcmp(info, 'cancelled')
                return
            end
            
            try
                if ~isempty(s)
                    % Use cached sr
                    sr = s.sr;
                    sr.mdat.Filename = fullfile(this.ksFolder, 'temp_wh.dat'); % update the mdat mapping
                else
                    % Construct MTracer.KilosortResult
                    sr = MTracer.KilosortResult();
                    assert(this.hasChanMap, 'Channel Map must be loaded before importing sorting results.');
                    sr.ImportData(ksDir, this.chanMapFile);
                    sr.ComputeAll();
                    
                    % Cache the new sr object
                    cacheFolder = fullfile(ksDir, this.cacheFolderName);
                    if ~exist(cacheFolder, 'dir')
                        mkdir(cacheFolder);
                    end
                    cacheFile = fullfile(cacheFolder, ['sr_' this.recId '_' datestr(now, 'yyyy-mm-dd_HH-MM-SS') '.mat']);
                    save(cacheFile, 'sr', '-v7.3');
                end
                
                % Replace the old
                delete(this.clustering);
                delete(this.layers.clusters);
                this.clustering = MTracer.ClusteringVM(this, sr);
                this.clustering.UpdateAll();
                
            catch e
                assignin('base', 'e', e);
                disp(e);
                if this.hasApp
                    uialert(this.app.UIFigure, 'Unable to find or load cluster data', 'Hmmm...', 'Icon', 'warning')
                end
            end
        end
        
        function [s, info] = LoadCache(~, pattern)
            % Check caches and load the selected file
            
            % Default return values
            s = [];
            info = '';
            
            % Check for caches
            cacheSearch = MBrowse.Dir2Table(pattern);
            if isempty(cacheSearch)
                return
            end
            
            options = [string(cacheSearch.name); "From original files"];
            [idx, isSelected] = listdlg( ...
                'ListString', options, ...
                'Name', 'Select a cache or recompute', ...
                'SelectionMode', 'single', ...
                'ListSize', [300 160]);
                
            if ~isSelected
                info = 'cancelled';
                return
            end
            
            if idx == numel(options)
                info = 'recompute';
            else
                info = fullfile(cacheSearch.folder{idx}, cacheSearch.name{idx});
                s = load(info);
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
                    this.traces = [this.traces; obj];
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
            
            % Prompt user to select a parent folder
            parentFolder = uigetdir('', 'Select Parent Folder');
            if parentFolder == 0
                return; % User canceled the dialog
            end
            
            folderName = fullfile(parentFolder, ['MTracer_' this.recId '_' datestr(now, 'yyyy-mm-dd_HH-MM-SS')]);
            mkdir(folderName);
            
            for i = 1 : numel(this.traces)
                dataTb = this.traces(i).dataTb;
                save(fullfile(folderName, [this.traces(i).fileName '.mat']), 'dataTb');
            end
            
            disp('Tracing results are saved to the selected directory.');
            if this.hasApp
                uialert(this.app.UIFigure, 'Tracing results are saved to the selected directory.', 'Success', 'Icon', 'Info');
            end
        end
        
        function tracerName = LoadAutoTracer(this)
            % Load an auto-tracer from file
            filePath = MBrowse.File([], 'Please select a file with an auto-tracer object', '*.mat');
            if ~exist(filePath, 'file')
                tracerName = [];
                return
            end
            s = load(filePath);
            [~, tracerName] = fileparts(filePath);
            fn = fieldnames(s);
            assert(numel(fn)==1, "Tracer file should only contain a single tracer object variable.");
            if ismember(tracerName, fieldnames(this.tracers))
                warning("Loading this auto-tracer overwrites the existing one with the same file name.")
            end
            obj = s.(fn{1});
            obj.vm = this;
            this.tracers(1).(tracerName) = obj;
            this.currentTracer = tracerName;
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
            p.addOptional('tWin', [], @(x) numel(x)==2 && isnumeric(x));
            p.addOptional('yWin', [], @(x) numel(x)==2 && isnumeric(x));
            p.addParameter('DataNames', {'traces', 'spikes', 'LFP', 'AP'});
            p.addParameter('TraceType', 'interp', @(x) ismember(x, {'anchor', 'interp'}));
            p.parse(varargin{:});
            tWin = p.Results.tWin;
            yWin = p.Results.yWin;
            dn = cellstr(p.Results.DataNames);
            traceType = p.Results.TraceType;
            
            if isempty(tWin)
                tWin = this.hAxes.XLim;
            end
            if isempty(yWin)
                tWin = this.hAxes.YLim;
            end
            
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
            
            if this.isOpen
                cla(this.hAxes);
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
        
        function cfg = SaveLayout(this)
            % Save window and axes positions
            
            if this.hasApp
                cfg.layout.app.figPos = this.app.UIFigure.Position;
            else
                cfg.layout.app.figPos = [];
            end
            cfg.layout.main = this.positions;
            cfg.layout.ftFig = this.clustering.ftFig.positions;
            cfg.layout.ccgFig = this.clustering.ccgFig.positions;
            cfg.layout.waveFig = this.clustering.waveFig.positions;
            
            [fileName, fileDir] = uiputfile('*.mat', "Save configuration file", 'MTracer_config.mat');
            if fileName
                save(fullfile(fileDir, fileName), 'cfg');
            end
        end
        
        function LoadLayout(this)
            % Load and apply window and axes positions
            
            filePath = MBrowse.File([], "Select a saved configuration file", {'*.mat'});
            if ~exist(filePath, 'file')
                return
            end
            load(filePath, 'cfg');
            
            if this.hasApp && ~isempty(cfg.layout.app.figPos)
                this.app.UIFigure.Position(1:2) = cfg.layout.app.figPos(1:2);
            end
            this.positions = cfg.layout.main;
            this.clustering.ftFig.positions = cfg.layout.ftFig;
            this.clustering.ccgFig.positions = cfg.layout.ccgFig;
            this.clustering.waveFig.positions = cfg.layout.waveFig;
        end
        
        % Display
        function InitializeMapAxes(this)
            
            % Clear stored Axes position since (maybe?) it interfereces with LooseInset
            this.axPos = [];
            
            % Create main figure if absent
            [f, ax] = this.Open( ...
                'Name', 'MTracer: Maps', ...
                'NumberTitle', 'off', ...
                'IntegerHandle', 'off', ...
                'Menubar', 'none', ...
                'Toolbar', 'figure', ...
                'WindowKeyPressFcn', @this.KeyPress, ...
                'WindowKeyReleaseFcn', @this.KeyRelease, ...
                'WindowScrollWheelFcn', @this.Scroll, ...
                'Interruptible', 'off', ...
                'BusyAction', 'cancel');
            
            % Create and format new map axes
            ax.XLabel.String = 'Time (sec)';
            ax.YLabel.String = 'Distance from tip (um)';
            ax.XLim = this.tLims;
            ax.YLim = this.yLims;
            ax.TickLength(1) = 0.002;
            ax.LooseInset = [0 0 0 0];
            ax.ButtonDownFcn = @this.SetPoint;
            ax.BusyAction = 'cancel';
            ax.Tag = 'Y-Position'; % tag axes by the feature name
            hold(ax, 'on');
            MPlot.Axes(ax);
            
            if isempty(this.layers)
                disp("This may be an old MTracerVM? Initialize layers struct");
                this.layers.focus = [];      % 1
                this.layers.spikes = [];     % 2
                this.layers.clusters = [];   % 3
                this.layers.AP = [];         % 4
                this.layers.LFP = [];        % 5
                this.layers.anchors = [];    % 6
                this.layers.interp = [];     % 7
            end
        end
        
        function PlotAll(this)
            this.InitializeMapAxes();
            this.PlotLFP();
            this.PlotRezSpikes();
            this.clustering.PlotSpikes();
            this.PlotTraces();
            this.PlotFocus();
        end
        
        function PlotFocus(this)
            % Plot the focus indicator
            
            if this.hasApp
                this.app.FocusTimeEditField.Value = this.focus(1);
                this.app.FocusDepthEditField.Value = this.focus(2);
            end
            
            if ~this.isOpen
                return
            end
            
            xx = [this.focus([1 1])' this.tLims'];
            yy = [this.yLims' this.focus([2 2])'];
            
            switch this.mapMode
                case 'selection'
                    cc = [0 1 1 .7];
                case 'tracing'
                    cc = [0 1 0 .5];
                case 'clustering'
                    cc = [0 0 1 .5];
            end
            
            h = this.layers.focus;
            if isempty(h) || any(~isvalid(h))
                h = plot(xx, yy, 'Color', cc, 'LineWidth', 1, 'HitTest', 'off');
                delete(this.layers.focus);
                this.AddHandle2Layer(h, 'focus');
            else
                set(h(1), 'XData', this.focus([1 1])', 'YData', this.yLims', 'Color', cc);
                set(h(2), 'XData', this.tLims', 'YData', this.focus([2 2])', 'Color', cc);
            end
            uistack(h, 'top', Inf);
        end
        
        function PlotRezSpikes(this)
            % Plot spike events as a function of time and depth
            
            if ~this.hasRez || ~this.isOpen
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
                hh{k} = plot(this.hAxes, t(ind), y(ind), '.', ...
                    'MarkerSize', this.spikeMarkerSize, ...
                    'Color', [1 1 1] * max(0, 1-ampRange(k)/40), ... % the marker color here has been carefully tuned
                    'HitTest', 'off');
            end
            
            MPlot.Axes(this.hAxes);
            this.SetMapROI(this.tLims, this.yLims);
            
            delete(this.layers.spikes);
            this.AddHandle2Layer(cat(1, hh{:}), 'spikes');
        end
        
        function PlotLFP(this)
            % Plot LFP as a heatmap
            
            if ~this.hasLFP || ~this.isOpen
                return
            end
            
            tb = this.se.GetTable('LFP');
            C = tb.v{1}';
            x = tb.time{1};
            y = this.channelTable.ycoords;
            
            h = imagesc(this.hAxes, x, y, C, 'HitTest', 'off');
            colormap(this.hAxes, MPlot.PolarMap);
            MPlot.Axes(this.hAxes);
            this.hAxes.YDir = 'normal';
            this.SetMapROI(this.tLims, this.yLims);
            
            delete(this.layers.LFP);
            this.AddHandle2Layer(h, 'LFP');
        end
        
        function PlotAP(this)
            % Plot AP signal
            
            if ~this.isOpen
                return
            end
            
            % Read a temporal slice of AP data from file
            this.ReadApSlice();
            if ~ismember('AP', this.se.tableNames)
                return
            end
            tb = this.se.GetTable('AP');
            v = tb.v{1};
            t = tb.time{1};
            y = tb.Properties.UserData.y;
            
            % Limit y window to 2000um maximal
            yWin = this.hAxes.YLim;
            if diff(yWin) > 2000
                yWin = this.focus(2) + [-1 1]*1000;
            end
            
            % Choose channel interval based on y window size
            dy = abs(y(2)-y(1));
            iy = MMath.Bound(diff(yWin)/dy/50, [1 2 3 4]);
            isChan = y > yWin(1) & y <= yWin(2);
            chanInd = find(isChan, 1, 'first') : iy : find(isChan, 1, 'last');
            
            % Use every other sample in time
            tInd = 1 : 2 : numel(t);
            t = t(tInd);
            v = v(tInd, chanInd);
            y = y(chanInd);
            
            % Highpass filtering
            isHp = strcmp(this.apSource, 'imec.ap.bin');
            if isHp
                % Prepare parameters for high-pass filtering
                %   Phy uses bandpass between 500Hz and 14.25kHz(.475*sample_rate) to visulize waveform;
                %   here uses 300Hz highpass to be consistent with Kilosort configuration
                D = designfilt('highpassiir', ...
                    'PassbandFrequency', 300, ...
                    'StopbandFrequency', 250, ...
                    'StopbandAttenuation', 60, ...
                    'PassbandRipple', 0.1, ...
                    'SampleRate', 30e3, ...
                    'DesignMethod', 'ellip'); % order of this filter is 8
                
%                 v = double(v);
                for i = 1 : size(v,2)
                    v(:,i) = filtfilt(D, v(:,i));
                end
            end
            
            % Scale voltage timeseries
            v = zscore(v) * 3;
            
            % Plot
            hh = line(t, v + y', 'Color', [0 0 0 .3], 'Parent', this.hAxes, 'HitTest', 'off');
            delete(this.layers.AP)
            this.AddHandle2Layer(hh, 'AP');
        end
        
        function SetMapROI(this, tWin, yWin, t, y)
            % Go to the ROI on the Map and update related app UI components
            % 
            %   this.SetMapROI(tWin, yWin)
            %   this.SetMapROI(tWin, yWin, t, y)
            % 
            % When reference coordinates t and y are not provided, tWin and yWin are used as absolute positions.
            % When t and y are provided, tWin and yWin are used as relative positions (i.e. diff(tWin) and diff(yWin))
            %   and are re-centered around (t,y).
            % If the given or derived tWin or yWin ends up out of limit (exceeding this.tLims or this.yLims), 
            %   the ROI size will be maintained by adjusting the center.
            % 
            
            if nargin > 3
                % Recenter ROI around (t,y)
                tWin = tWin - mean(tWin) + t;
                yWin = yWin - mean(yWin) + y;
            end
            
            % Check and adjust out-of-limit in time-axis
            if tWin(1) < this.tLims(1)
                tWin = [this.tLims(1), this.tLims(1)+diff(tWin)];
            end
            if tWin(2) > this.tLims(2)
                tWin = [this.tLims(2)-diff(tWin), this.tLims(2)];
            end
            tWin = MMath.Bound(tWin, this.tLims);
            
            % Check and adjust out-of-limit in depth-axis
            if yWin(1) < this.yLims(1)
                yWin = [this.yLims(1) this.yLims(1)+diff(yWin)];
            end
            if yWin(2) > this.yLims(2)
                yWin = [this.yLims(2)-diff(yWin), this.yLims(2)];
            end
            yWin = MMath.Bound(yWin, this.yLims);
            
            % Apply axes limits
            if this.isOpen
                set(this.hAxes, 'XLim', tWin, 'YLim', yWin);
            end
            if this.hasClus && this.clustering.ftFig.isOpen
                set(this.clustering.ftFig.hAxes, 'XLim', tWin);
            end
            
            if this.hasApp
                % Update app UI components
                tWin = round(tWin, 3);
                yWin = round(yWin, 3);
                this.app.TimeWin1EditField.Value = tWin(1);
                this.app.TimeWin2EditField.Value = tWin(2);
                this.app.DepthWin1EditField.Value = yWin(1);
                this.app.DepthWin2EditField.Value = yWin(2);
            end
        end
        
        function UpdateSpikeMarkerSize(this, r, valType)
            
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
            if this.isOpen && this.hasRez
                hh = this.layers.spikes;
                for i = 1 : numel(hh)
                    hh(i).MarkerSize = sz;
                end
            end
            
            % Apply marker size to clusters
            if this.isOpen && this.hasClus
                hh = this.layers.clusters;
                for i = 1 : numel(hh)
                    hh(i).MarkerSize = sz;
                end
            end
            
            % 
            if this.hasClus
                this.clustering.UpdateSpikeMarkerSize(sz);
            end
        end
        
        % Tracing
        function PlotTraces(this)
            % Plot traces
            if this.hasApp
                this.app.TraceListBox.Items = arrayfun(@(x) x.dispName, this.traces, 'Uni', false);
            end
            if this.isOpen
                for i = 1 : numel(this.traces)
                    this.traces(i).PlotTrace(ismember(i, this.currentTrace));
                end
            end
        end
        
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
                this.mapMode = 'tracing';
                this.PlotFocus();
                disp('Tracing Mode');
            else
                warning('Cannot enter tracing mode with more than one selected trace\n');
            end
        end
        
        function EnterClusteringMode(this)
            if this.hasClus
                this.mapMode = 'clustering';
                this.PlotFocus();
                disp('Clustering Mode');
            else
                fprintf("Cannot enter clustering mode since there is no cluster data\n");
            end
        end
        
        function EnterSelectionMode(this)
            if this.hasClus
                this.clustering.plg.ClearPolygon();
            end
            this.mapMode = 'selection';
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
        
        function MergeTraces(this)
            % Merge selected traces into one
            if numel(this.currentTrace) < 2
                return
            end
            tb = cat(1, this.traces(this.currentTrace).dataTb);
            tb = sortrows(tb, 't');
            tr = MTracerTrace(this, tb);
            this.traces(this.currentTrace) = [];
            this.traces = [this.traces; tr];
            this.currentTrace = numel(this.traces);
            this.PlotTraces();
        end
        
        % User Input Callbacks
        function KeyPress(this, ~, eventdata)
            
            K = eventdata.Key;
            M = eventdata.Modifier;
            
            if ismember('control', M)
                switch K
                    case 'z'
                        % Undo/Redo
                        if strcmp(this.mapMode, 'tracing') && ~isnan(this.currentTrace)
                            if ismember('shift', M)
                                this.traces(this.currentTrace).Redo();
                            else
                                this.traces(this.currentTrace).Undo();
                            end
                        end
                end
            elseif any(ismember({'alt', 'option'}, M))
                switch K
                    % Set cluster label
                    case 'g', this.clustering.LabelClusters('good');
                    case 'm', this.clustering.LabelClusters('mua');
                    case 'n', this.clustering.LabelClusters('noise');
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
                tWin = this.hAxes.XLim;
                yWin = this.hAxes.YLim;
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
                    case 'q'
                        this.AutoTrace();
                        
                    case 'k'
                        if strcmp(this.mapMode, 'selection')
                            this.EnterClusteringMode();
                        elseif strcmp(this.mapMode, 'clustering')
                            this.clustering.CutClusters();
                        end
                    case 'escape'
                        this.EnterSelectionMode();
                    case 'f'
                        this.clustering.NextFeature();
                        
                    case {'a', 'leftarrow'}
                        % Forward in time
                        tWin = tWin - dt;
                        this.SetMapROI(tWin, yWin);
                    case {'d', 'rightarrow'}
                        % Backward in time
                        tWin = tWin + dt;
                        this.SetMapROI(tWin, yWin);
                    case {'s', 'downarrow'}
                        % Downward in depth
                        yWin = yWin - dy;
                        this.SetMapROI(tWin, yWin);
                    case {'w', 'uparrow'}
                        % Upward in depth
                        yWin = yWin + dy;
                        this.SetMapROI(tWin, yWin);
                        
                    case 'leftbracket'
                        % Zoom-in in time
                        tWin = (tWin-t)/0.5 + t;
                        this.SetMapROI(tWin, yWin, t, y);
                    case 'rightbracket'
                        % Zoom-out in time
                        tWin = (tWin-t)*0.5 + t;
                        this.SetMapROI(tWin, yWin, t, y);
                    case {'equal', 'add'}
                        % Zoom-in in depth
                        yWin = (yWin-y)*0.8 + y;
                        this.SetMapROI(tWin, yWin, t, y);
                    case {'hyphen', 'subtract'}
                        % Zoom-out in depth
                        yWin = (yWin-y)/0.8 + y;
                        this.SetMapROI(tWin, yWin, t, y);
                        
                    case 'c'
                        % Center the reference point
                        this.SetMapROI(tWin, yWin, t, y);
                    case 'o'
                        % View full map
                        this.SetMapROI(this.tLims, this.yLims);
                        
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
        
        function KeyRelease(this, ~, eventdata)
            
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
        
        function Scroll(this, ~, eventdata)
            tWin = this.hAxes.XLim;
            dt = eventdata.VerticalScrollCount * eventdata.VerticalScrollAmount/3 * diff(tWin) / 50;
            this.SetMapROI(tWin+dt, this.hAxes.YLim);
        end
        
        function SetPoint(this, src, eventdata)
            % Place a point on image
            
            % Read mouse position
            mousePos = get(src, 'CurrentPoint');
            
            % Update reference point
            t = MMath.Bound(mousePos(1), this.tLims);
            y = MMath.Bound(mousePos(3), this.yLims);
            this.focus = [t y];
            this.PlotFocus();
            
            if strcmp(this.mapMode, 'tracing')
                % Add or remove trace point
                if isnan(this.currentTrace)
                    set(this.layers.anchors, 'Visible', 'on');
                    set(this.layers.interp, 'Visible', 'on');
                    this.InitializeTrace();
                end
                tr = this.traces(this.currentTrace);
                if eventdata.Button == 1
                    tr.SetPoint(t, y);
                elseif eventdata.Button == 3
                    tr.RemovePoint(t, y);
                end
            elseif strcmp(this.mapMode, 'clustering')
                % Add or remove polygon vertex
                if eventdata.Button == 1
                    this.clustering.plg.AddPoint(src, t, y);
                elseif eventdata.Button == 3
                    this.clustering.plg.RemoveLastPoint();
                end
            else
                % Trace selection
                this.SelectTrace(t, y);
            end
        end
        
        function AutoTrace(this, src, eventdata)
            % 
            if isnan(this.currentTrace) || strcmp(this.currentTracer, 'none')
                return
            end
            [t, y] = this.traces(this.currentTrace).AutoSetPoint();
            this.focus = [t y];
            this.PlotFocus();
        end
        
    end
    
end