classdef ClusteringVM < handle
    % 
    
    properties
        vm MTracerVM;                   % handle of the main app view-model
        sr;                             % object managing sorting results
        om MTracer.OperationManager;    % object managing clustering operation history
        coi = [];       % clusters of interest, usually aka the selected clusters
        cutCoords;      % a n-by-2 numeric array that stores the (t,y) coordinates of the cutting polygon
        hCut;           % handle of the cutting polygon
    end
    properties(Constant)
        appTbColNames = {'clusId', 'group', 'depth', 'numSpikes', 'RPV', 'contam', 'SNR'};
        noiseColor = [0 0 0] + .5;
        maxClusOp = 5;
    end
    properties(Dependent)
        hasClus;
        hasBin;
        isShowNoise;
    end
    methods
        function val = get.hasClus(this)
            val = ~(isempty(this.sr) || isempty(this.sr.spkTb));
        end
        function val = get.hasBin(this)
            val = ~isempty(this.sr) && ~isempty(this.sr.mdat);
        end
        function val = get.isShowNoise(this)
            if this.vm.hasApp
                val = this.vm.app.ShowNoiseCheckBox.Value;
            else
                val = false;
            end
        end
    end
    
    methods
        % Object construction
        function this = ClusteringVM(vm, sr)
            % Constructor
            this.vm = vm;
            if nargin > 1
                this.sr = sr;
                this.om = MTracer.OperationManager(this.GetStateData([], []));
                this.SetClusterColor();
            end
        end
        
        function vmClus = Duplicate(this, vmApp)
            % Make a hard copy of the object
            
            vmClus = MTracer.ClusteringVM(vmApp);
            
            if ~isempty(this.sr)
                return
            end
            
            % Make a hard copy
            srCopy = this.sr.Duplicate();

            % Remove added columns
            isAdded = ismember(srCopy.clusTb.Properties.VariableNames, {'handle', 'memmap'});
            srCopy.clusTb(:,isAdded) = [];

            vmClus.sr = srCopy;
            vmClus.om = this.om;
            vmClus.coi = this.coi;
        end
        
        function SaveObject(this, filePath)
            % Save sorting result object
            
            if ~this.hasClus
                return
            end
            
            % Make a hard copy
            srCopy = this.sr.Duplicate();
            
            % Remove added columns
            isAdded = ismember(srCopy.clusTb.Properties.VariableNames, {'handle', 'memmap'});
            srCopy.clusTb(:,isAdded) = [];
            
            % Remove clusters without spikes
            srCopy.clusTb(~srCopy.clusTb.numSpikes,:) = [];
            
            % Save object
            if nargin < 2 || isempty(filePath)
                filePath = fullfile(this.vm.ksFolder, ['sr_' this.vm.recId '_' datestr(now, 'yyyy-mm-dd_HH-MM-SS')]);
            end
            sr = srCopy;
            uisave('sr', filePath);
        end
        
        function SaveData(this)
            % Save clustering result to files
            if ~this.hasClus
                return
            end
            this.sr.ExportData();
        end
        
        function delete(this)
            % Destructor
            delete(this.hCut);
        end
        
        % Cluster operations
        function s = GetStateData(this, preCOI, postCOI)
            % Assemble a struct of state data for this.om to manage
            s.spkClusId = this.sr.spkTb.clusId;
            s.preClusId = preCOI(:)';
            s.postClusId = postCOI(:)';
        end
        
        function SetClusterColor(this)
            % Assign color to clusters
            %   1) Colors are generated by lines()
            %   2) To initialize, colors are assigned to clusters in the order of descending depth 
            %      such that clusters close in space can maximally avoid the same color.
            %   3) Black [0 0 0] is reserved for unassigned clusters, and will be assigned a color 
            %      upon calling this method.
            %   4) Once assigned, the color of a cluster will never be changed (unless by users).
            %   5) When assigning color to an unassigned (e.g. new) cluster, a color that is different 
            %      from its 6 closest clusters will be used.
            %   6) 'noise' clusters, though having colors, are always plotted in gray (this.noiseColor).
            
            % Find the forward and backward sorting indices
            tb = this.sr.clusTb;
            [~, fwInd] = sort(tb.depth, 'descend');
            [~, bwInd] = sort(fwInd);
            
            % Initialize colors if empty
            if ~ismember('color', tb.Properties.VariableNames)
                cc = lines(height(tb));
                tb.color = cc(bwInd,:);
            end
            
            % Find unassigned clusters (i.e. [0 0 0])
            indNoColor = find(~any(tb.color, 2));
            
            % Assign colors
            cMap = lines(7);
            for i = indNoColor(:)'
                % Find the colors of the 6 closest clusters
                d = abs(tb.depth - tb.depth(i));
                [~, I] = sort(d);
                cUsed = tb.color(I(2:min(7,end)),:);
                
                % Find and assign the first unused color
                isUsed = ismember(cMap, cUsed, 'rows');
                k = find(~isUsed, 1);
                tb.color(i,:) = cMap(k,:);
            end
            
            this.sr.clusTb = tb;
        end
        
        function cc = GetClusterColor(this, cid)
            % Return cluster color in a n-by-3 array. Noise clusters are always in gray (this.noiseColor).
            tb = this.sr.clusTb;
            cc = zeros(numel(cid), 3);
            for i = 1 : numel(cid)
                m = tb.clusId==cid(i);
                if strcmp(tb.group{m}, 'noise')
                    c = this.noiseColor;
                else
                    c = tb.color(m, :);
                end
                cc(i,:) = c;
            end
        end
        
        function LabelClusters(this, newGroup)
            % Change the group label of specified clusters
            cid = this.coi;
            if ~this.hasClus || numel(cid) < 1 || numel(cid) > this.maxClusOp
                return
            end
            this.sr.LabelClusters(cid, newGroup);
            this.UpdateAll();
        end
        
        function MergeClusters(this, cid)
            % Merge selected clusters and update UI
            if nargin < 2
                cid = this.coi;
            end
            if ~this.hasClus || numel(cid) < 2 || numel(cid) > this.maxClusOp
                return
            end
            
            newId = this.sr.MergeClusters(cid);
            this.SetClusterColor();
            
            s = this.GetStateData(cid, newId);
            this.om.AddOperation(s);
            
            this.coi = newId;
            this.UpdateAll();
        end
        
        function CutClusters(this)
            % Cut spikes encircled by the polygon and update UI
            if ~this.hasClus
                return
            end
            
            spkInd = this.FindSpikesInPolygon();
            this.ClearPolygon();
            if isempty(spkInd)
                return
            end
            
            [newId, reId, oldId] = this.sr.CutClusters(spkInd);
            this.SetClusterColor();
            
            s = this.GetStateData(oldId, [newId reId]);
            this.om.AddOperation(s);
            
            this.coi = [newId reId];
            this.UpdateAll();
            this.vm.EnterSelectionMode();
        end
        
        function Undo(this)
            % Undo the last operation
            
            [fromState, toState] = this.om.Undo();
            if isempty(toState)
                return
            end
            
            % Reassign the previous set of spike cluster IDs
            this.sr.spkTb.clusId = toState.spkClusId;
            
            % Recompute the number of spikes
            this.sr.ComputeNumSpikes([fromState.preClusId, fromState.postClusId]);
            
            % Select clusters
            this.coi = fromState.preClusId;
            this.UpdateAll();
        end
        
        function Redo(this)
            % Redo the next operation
            
            [~, toState] = this.om.Redo();
            if isempty(toState)
                return
            end
            
            % Reassign the previous set of spike cluster IDs
            this.sr.spkTb.clusId = toState.spkClusId;
            
            % Recompute the number of spikes
            this.sr.ComputeNumSpikes([toState.preClusId, toState.postClusId]);
            
            % Select clusters
            this.coi = toState.postClusId;
            this.UpdateAll();
        end
        
        function spkInd = FindSpikesInPolygon(this)
            % Return row indices in spkTb for spikes encircled by the polygon
            if isempty(this.coi) || size(this.cutCoords, 1) < 3
                spkInd = [];
                return
            end
            m = ismember(this.sr.spkTb.clusId, this.coi);
            spkInd = find(m);
            xq = this.sr.spkTb.timeSec(spkInd);
            yq = this.sr.spkTb.centCoords(spkInd,2);
            xv = this.cutCoords(:,1);
            yv = this.cutCoords(:,2);
            in = inpolygon(xq, yq, xv, yv);
            spkInd = spkInd(in);
        end
        
        function SetPoint(this, tPt, yPt)
            % Add a new point
            this.cutCoords = [this.cutCoords; tPt yPt];
            this.PlotPolygon();
        end
        
        function RemoveLastPoint(this)
            % Remove the last point
            if size(this.cutCoords, 1) > 0
                this.cutCoords(end,:) = [];
                this.PlotPolygon();
            end
        end
        
        % Interaction
        function InitializeAppClusTable(this)
            % Initialize ClusTable in the Clustering tab
            if ~this.vm.hasApp
                return
            end
            app = this.vm.app;
            app.ClusTable.ColumnName = this.appTbColNames;
            app.ClusTable.ColumnWidth = repmat({'1x'}, size(this.appTbColNames));
            app.ClusTable.SelectionType = 'row';
        end
        
        function UpdateAppClusTable(this)
            % Update ClusTable in the Clustering tab
            
            if ~this.vm.hasApp
                return
            end
            app = this.vm.app;
            if ~this.hasClus
                app.ClusTable.Data = [];
                return
            end
            
            % Prepare the table
            tb = this.sr.clusTb;
            tb = tb(tb.numSpikes > 0, [this.appTbColNames {'color'}]);
            isNoise = strcmp(tb.group, 'noise');
            if this.isShowNoise
                tb.color(isNoise,:) = repmat(this.noiseColor, [sum(isNoise) 1]);
            else
                tb(isNoise,:) = [];
            end
            tb = sortrows(tb, 'depth', 'descend');
            
            % Set the table
            app.ClusTable.Data = tb(:, this.appTbColNames);
            
            % Apply text colors
            removeStyle(app.ClusTable);
            [cc, ~, ic] = unique(tb.color, 'rows');
            for i = 1 : size(cc, 1)
                s = uistyle("FontColor", cc(i,:));
                addStyle(app.ClusTable, s, 'row', find(ic==i));
            end
            
            % Set selected row(s) and scroll to the first item
            ind = find(ismember(tb.clusId, this.coi));
            app.ClusTable.Selection = ind(:)'; % must be a row vector
            if ~isempty(ind)
                scroll(app.ClusTable, 'row', max(1, ind(1)-2));
            end
        end
        
        function UpdateCOI(this)
            % 
            this.PlotSpikes();
            if ~this.vm.hasApp
                return
            end
            app = this.vm.app;
            if app.ROIfollowCheckBox.Value
                this.FindOnMap();
            end
            if app.AutoCCGCheckBox.Value
                this.PlotCCG();
            end
            if app.AutoWaveformCheckBox.Value
                this.PlotWaveform();
            end
            if app.AutoAmpCheckBox.Value
                this.PlotAmplitude();
            end
        end
        
        function UpdateAll(this)
            % 
            this.UpdateCOI();
            this.UpdateAppClusTable();
        end
        
        function FindOnMap(this, cid)
            % Zoom to the selected cluster(s) in the Maps Window
            
            if ~this.hasClus || ~this.vm.hasMapAxes
                return
            end
            if nargin < 2 || isempty(cid)
                cid = this.coi;
            end
            if isempty(cid)
                cid = this.sr.clusTb.clusId;
            end
            
            % Find the range of mean cluster depths
            m = ismember(this.sr.clusTb.clusId, cid);
            yClus = this.sr.clusTb.depth(m);
            yWin = [min(yClus) max(yClus)] + [-1 1]*200;
            
            % Set ROI
            this.vm.SetMapROI(this.vm.mapAxes.XLim, yWin);
        end
        
        function FindInTable(this)
            % Select the cluster that's closest in depth to the focus
            if ~this.hasClus || ~this.vm.hasApp || ~this.vm.hasMapAxes
                return
            end
            y = this.vm.app.ClusTable.Data.depth;
            y0 = this.vm.focus(2);
            [~, I] = sort(abs(y - y0));
            cid = this.vm.app.ClusTable.Data.clusId(I(1));
            this.coi = cid;
            this.UpdateAll();
        end
        
        % Plotting
        function PlotSpikes(this)
            % Plot all or selected clusters
            
            if ~this.hasClus || ~this.vm.hasMapAxes
                return
            end
            
            % Plot handles should not be accessed from a table for it is extremly slow
            tb = this.sr.clusTb;
            if ~ismember('handle', tb.Properties.VariableNames)
                % Initialize plot handles with placeholders
                h = matlab.graphics.GraphicsPlaceholder();
                delete(h);
                hh = repmat(h, size(tb.clusId));
            else
                hh = tb.handle;
            end
            
            % Include all clusters when none is currently selected
            if isempty(this.coi)
                cid = tb.clusId;
            else
                cid = this.coi;
            end
            
            % Highlight included non-empty clusters
            isNoise = strcmp(tb.group, "noise");
            isHighlight = ismember(tb.clusId, cid) & tb.numSpikes > 0;
            if ~this.isShowNoise
                isHighlight = isHighlight & ~isNoise;
            end
            
            for i = 1 : height(tb)
                % Get cluster info
                k = tb.clusId(i);
                g = tb.group(i);
                c = tb.color(i,:);
                
                % Remove plot if the cluster is empty
                h = hh(i);
                if ~tb.numSpikes(i)
                    delete(h);
                    continue
                end
                
                % Overwrite the color of 'noise' cluster
                if isNoise(i)
                    c = this.noiseColor;
                end
                
%                 % Add alpha to unselected clusters, with "mua" lighter than "good"
%                 if strcmp(g, "good")
%                     r = 0.4;
%                 else % "mua"
%                     r = 0.2;
%                 end
%                 if ~isHighlight(i)
%                     w = ones(size(c));
%                     c = r*c + (1-r)*w;
%                 end
                
                % Make or update plot
                if isempty(properties(h)) || ~isvalid(h) % placeholder does not have properties
                    % Get spike data
                    m = this.sr.spkTb.clusId == k;
                    t = this.sr.spkTb.timeSec(m);
                    y = this.sr.spkTb.centCoords(m,2);
                    
                    % Plot spikes
                    h = plot(this.vm.mapAxes, t, y, ...
                        'Color', c, ...
                        'LineStyle', 'none', ...
                        'Marker', '.', ...
                        'MarkerSize', this.vm.spikeMarkerSize, ...
                        'HitTest', 'off', ...
                        'Visible', isHighlight(i));
                else
                    if any(c ~= h.Color)
                        h.Color = c;
                    end
                    h.Visible = isHighlight(i);
                end
                
                % Label line object
                if isHighlight(i)
                    h.UserData.clusId = k;
                    h.UserData.group = g;
                    h.UserData.color = c;
                else
                    h.UserData = [];
                end
                
                hh(i) = h;
            end
            
            % Move 'noise' clusters to the bottom
            uistack(hh(isHighlight & isNoise), 'bottom');
            
            % Save handles
            this.sr.clusTb.handle = hh;
            this.vm.mapLayers.clusters = hh(isvalid(hh));
            
            % Initialize brush callback
            brushObj = brush(this.vm.mapFig);
            brushObj.ActionPostCallback = @this.PlotBrushedSpikes;
        end
        
        function PlotBrushedSpikes(this, fig, axesStruct)
            % 
            
            if ~this.hasBin
                warning('Cannot plot waveform as the binary file is not available.');
                return
            end
            
            % Find clusters and spikes to plot
            hh = axesStruct.Axes.Children;
            cid = NaN(numel(hh), 1);
            cc = NaN(numel(hh), 3);
            spkMask = cell(size(cid));
            for i = 1 : numel(hh)
                % Check if the handle is cluster spike plot
                s = hh(i).UserData;
                if ~isfield(s, 'clusId')
                    continue
                end
                
                % Check if selecting any data
                b = logical(hh(i).BrushData);
                if ~any(b)
                    continue
                end
                
                cid(i) = s.clusId;
                cc(i,:) = s.color;
                spkMask{i} = b;
            end
            isBrushed = ~isnan(cid);
            if ~any(isBrushed)
                return
            end
            cid = cid(isBrushed);
            cc = cc(isBrushed,:);
            spkMask = spkMask(isBrushed);
            
            % Create figure if absent
            f = this.vm.waveFig;
            if isempty(f) || ~isvalid(f)
                f = MPlot.Figure( ...
                    'Name', 'Waveform', ...
                    'NumberTitle', 'off', ...
                    'Menubar', 'none', ...
                    'Toolbar', 'figure', ...
                    'IntegerHandle', 'off', ...
                    'Interruptible', 'off', ...
                    'BusyAction', 'cancel');
                this.vm.waveFig = f;
            end
            f.Position(3) = 200 * numel(cid);
            
            ax = f.Children;
            if isempty(ax) || ~isvalid(ax)
                ax = MPlot.Axes(f);
            end
            hold(ax, 'off');
            
            this.sr.PlotClusterWaveform(cid, 50, spkMask, 'NumChannels', 10, 'Color', cc, 'Axes', ax);
            ax = MPlot.Axes(gca);
            ax.Title.String = ['Cluster ' num2str(cid(:)')];
        end
        
        function PlotWaveform(this, cid)
            
            if ~this.hasBin
                warning('Cannot plot waveform as the binary file is not available.');
                return
            end
            if nargin < 2
                cid = this.coi;
            end
            if isempty(cid) || numel(cid) > this.maxClusOp
                return
            end
            
            f = this.vm.waveFig;
            if isempty(f) || ~isvalid(f)
                % Create figure if absent
                f = MPlot.Figure( ...
                    'Name', 'Waveform', ...
                    'NumberTitle', 'off', ...
                    'Menubar', 'none', ...
                    'Toolbar', 'none', ...
                    'IntegerHandle', 'off', ...
                    'Interruptible', 'off', ...
                    'BusyAction', 'cancel');
                this.vm.waveFig = f;
            end
            f.Position(3) = 200 * numel(cid);
            
            ax = f.Children;
            if isempty(ax) || ~isvalid(ax)
                ax = MPlot.Axes(f);
            end
            hold(ax, 'off');
            
            this.sr.PlotClusterWaveform(cid, 50, 'NumChannels', 10, 'Color', this.GetClusterColor(cid), 'Axes', ax);
            this.sr.PlotClusterTemplate(cid, 'Axes', ax);
            ax.Title.String = ['Cluster ' num2str(cid(:)')];
            ax.LooseInset([2 4]) = [0 0];
        end
        
        function PlotCCG(this, cid)
            
            if nargin < 2
                cid = this.coi;
            end
            if ~this.hasClus || isempty(cid) || numel(cid) > this.maxClusOp
                return
            end
            
            f = this.vm.ccgFig;
            if isempty(f) || ~isvalid(f)
                % Create figure if absent
                f = MPlot.Figure( ...
                    'Name', 'Cross-correlograms', ...
                    'NumberTitle', 'off', ...
                    'Menubar', 'none', ...
                    'Toolbar', 'none', ...
                    'IntegerHandle', 'off', ...
                    'Interruptible', 'off', ...
                    'BusyAction', 'cancel');
                this.vm.ccgFig = f;
            end
            
            this.sr.PlotCCG(cid, 'Color', this.GetClusterColor(cid), 'ShowMetrics', true, 'Figure', f);
        end
        
        function PlotFeatureTime(this, cid, featName)
            
            if nargin < 2
                cid = this.coi;
            end
            if ~this.hasClus || isempty(cid) || numel(cid) > this.maxClusOp
                return
            end
            
            f = this.vm.ftFig;
            if isempty(f) || ~isvalid(f)
                % Create figure if absent
                f = MPlot.Figure( ...
                    'Name', featName, ...
                    'NumberTitle', 'off', ...
                    'Menubar', 'none', ...
                    'Toolbar', 'figure', ...
                    'IntegerHandle', 'off', ...
                    'Interruptible', 'off', ...
                    'BusyAction', 'cancel');
                this.vm.ftFig = f;
                if this.vm.hasMapAxes
                    f.Position([1 3]) = this.vm.mapFig.Position([1 3]);
                    f.Position(4) = 150;
                end
            end
            
            ax = f.Children;
            if isempty(ax) || ~isvalid(ax)
                ax = axes(f);
            end
            hold(ax, 'off');
            
            cTb = this.sr.clusTb;
            for i = 1 : numel(cid)
                % Get cluster info
                k = find(cTb.clusId == cid(i));
                g = cTb.group(k);
                if strcmp(g, 'noise')
                    c = this.noiseColor; % overwrite the color of 'noise' cluster
                else
                    c = cTb.color(k,:);
                end
                
                if ~cTb.numSpikes(k)
                    continue
                end
                
                % Get spike data
                m = this.sr.spkTb.clusId == cid(i);
                t = this.sr.spkTb.timeSec(m);
                switch featName
                    case 'Amplitude'
                        a = this.sr.spkTb.tempAmp(m);
                    case 'PC1'
                        
                    case 'PC2'
                        
                end
                
                % Plot spikes
                h = plot(ax, t, a, ...
                    'Color', c, ...
                    'LineStyle', 'none', ...
                    'Marker', '.', ...
                    'MarkerSize', this.vm.spikeMarkerSize, ...
                    'HitTest', 'off');
                hold(ax, 'on');
                
                % Label line object
                h.UserData.clusId = cid(i);
                h.UserData.group = g;
                h.UserData.color = c;
            end
            
            MPlot.Axes(ax);
            if this.vm.hasMapAxes
                ax.XLim = this.vm.mapAxes.XLim;
            end
            ax.TickLength(1) = 0.002;
            ax.YLabel.String = 'Amplitude (au)';
            ax.LooseInset = [0 0 0 0];
            
            % Initialize brush callback
            brushObj = brush(this.vm.ampFig);
            brushObj.ActionPostCallback = @this.PlotBrushedSpikes;
        end
        
        function PlotAmplitude(this, cid)
            
            if nargin < 2
                cid = this.coi;
            end
            if ~this.hasClus || isempty(cid) || numel(cid) > this.maxClusOp
                return
            end
            
            f = this.vm.ampFig;
            if isempty(f) || ~isvalid(f)
                % Create figure if absent
                f = MPlot.Figure( ...
                    'Name', 'Amplitude', ...
                    'NumberTitle', 'off', ...
                    'Menubar', 'none', ...
                    'Toolbar', 'figure', ...
                    'IntegerHandle', 'off', ...
                    'Interruptible', 'off', ...
                    'BusyAction', 'cancel');
                this.vm.ampFig = f;
                if this.vm.hasMapAxes
                    f.Position([1 3]) = this.vm.mapFig.Position([1 3]);
                    f.Position(4) = 150;
                end
            end
            
            ax = f.Children;
            if isempty(ax) || ~isvalid(ax)
                ax = axes(f);
            end
            hold(ax, 'off');
            
            cTb = this.sr.clusTb;
            for i = 1 : numel(cid)
                % Get cluster info
                k = find(cTb.clusId == cid(i));
                g = cTb.group(k);
                if strcmp(g, 'noise')
                    c = this.noiseColor; % overwrite the color of 'noise' cluster
                else
                    c = cTb.color(k,:);
                end
                
                if ~cTb.numSpikes(k)
                    continue
                end
                
                % Get spike data
                m = this.sr.spkTb.clusId == cid(i);
                t = this.sr.spkTb.timeSec(m);
                a = this.sr.spkTb.tempAmp(m);
                
                % Plot spikes
                h = plot(ax, t, a, ...
                    'Color', c, ...
                    'LineStyle', 'none', ...
                    'Marker', '.', ...
                    'MarkerSize', this.vm.spikeMarkerSize, ...
                    'HitTest', 'off');
                hold(ax, 'on');
                
                % Label line object
                h.UserData.clusId = cid(i);
                h.UserData.group = g;
                h.UserData.color = c;
            end
            
            MPlot.Axes(ax);
            if this.vm.hasMapAxes
                ax.XLim = this.vm.mapAxes.XLim;
            end
            ax.TickLength(1) = 0.002;
            ax.YLabel.String = 'Amplitude (au)';
            ax.LooseInset = [0 0 0 0];
            
            % Initialize brush callback
            brushObj = brush(this.vm.ampFig);
            brushObj.ActionPostCallback = @this.PlotBrushedSpikes;
        end
        
        function PlotPolygon(this)
            % Plot cluster cutting polygon
            
            if size(this.cutCoords, 1) < 1
                return
            end
            
            t = this.cutCoords([1:end 1],1);
            y = this.cutCoords([1:end 1],2);
            
            ax = this.vm.mapAxes;
            if isempty(this.hCut) || ~isvalid(this.hCut)
                this.hCut = plot(ax, t, y, '-', 'LineWidth', 1, 'Color', [0 0 1]);
            else
                set(this.hCut, 'XData', t, 'YData', y);
            end
        end
        
        function ClearPolygon(this)
            this.cutCoords = [];
            delete(this.hCut);
        end
        
    end
    
end
