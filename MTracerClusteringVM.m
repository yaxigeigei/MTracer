classdef MTracerClusteringVM < handle
    % 
    
    properties
        vm MTracerVM;           % Handle of the main app model
        sr;                     % Object for accessing sorting results
        
        currentCluster = [];
        cutCoords;
        logTb table;
        
        hClusTb table;          % a table of handles for cluster spike plots
        hCut;                   % handle of the cutting polygon
        waveFig;                % handle of the waveform window
        ccgFig;                 % handle of the CCG window
        ampFig;                 % handle of the amplitude window
    end
    properties(Constant)
        appTbColNames = {'clusId', 'group', 'depth', 'numSpikes', 'SNR', 'RPV', 'contam'};
        noiseColor = [0 0 0] + .5;
        maxClusPlot = 5;
    end
    properties(Dependent)
        hasClus;
        hasBin;
        isShowNoise;
        clusTb;
    end
    methods
        function val = get.hasClus(this)
            val = ~isempty(this.sr);
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
        function tb = get.clusTb(this)
            if this.hasClus
                tb = this.sr.clusTb;
                
                isNoise = strcmp(tb.group, 'noise');
                hasSpk = ismember(tb.clusId, this.sr.spkTb.clusId);
                isShow = hasSpk & (~isNoise | this.isShowNoise);
                
                tb = tb(isShow, this.appTbColNames);
                tb = sortrows(tb, 'depth', 'descend');
                
                isNoise = strcmp(tb.group, 'noise');
                tb.color(~isNoise,:) = lines(sum(~isNoise));
                tb.color(isNoise,:) = repmat(this.noiseColor, [sum(isNoise) 1]);
            else
                tb = [];
            end
        end
    end
    
    methods
        % Object construction
        function this = MTracerClusteringVM(vm, sr)
            % Constructor
            this.vm = vm;
            if nargin > 1
                this.sr = sr;
            end
            this.hClusTb = table([], {}, 'VariableNames', {'clusId', 'spikes'});
        end
        
        function delete(this)
            % Destructor
            delete(this.hCut);
            delete(this.waveFig);
            delete(this.ccgFig);
            delete(this.ampFig);
        end
        
        % Data operations
        function UpdateClusHandleTb(this)
            % 
            if ~this.hasClus
                return
            end
            if isempty(this.hClusTb)
                this.hClusTb = table([], [], 'VariableNames', {'clusId', 'handle'});
            end
            
            oldId = this.hClusTb.clusId;
            oldHandle = this.hClusTb.handle;
            
            newId = this.clusTb.clusId;
            h = matlab.graphics.GraphicsPlaceholder();
            delete(h);
            newHandle = repmat(h, size(newId));
            
            for i = 1 : numel(oldId)
                m = newId == oldId(i);
                h = oldHandle(i);
                if ~any(m)
                    % Delete plot if the cluster is no longer displayed
                    delete(h);
                elseif ~isempty(h) && isvalid(h)
                    % Copy handle to new table
                    newHandle(m) = h;
                end
            end
            
            this.hClusTb = table(newId, newHandle, 'VariableNames', this.hClusTb.Properties.VariableNames);
        end
        
        function LabelClusters(this, newGroup)
            % Change the group label of specified clusters
            cid = this.currentCluster;
            if ~this.hasClus || numel(cid) < 1 || numel(cid) > this.maxClusPlot
                return
            end
            this.sr.LabelClusters(cid, newGroup);
            this.UpdateAppUI();
            this.PlotSpikes();
        end
        
        function MergeClusters(this, cid)
            % Merge selected clusters and update UI
            if nargin < 2
                cid = this.currentCluster;
            end
            if ~this.hasClus || numel(cid) < 2 || numel(cid) > this.maxClusPlot
                return
            end
            
            this.LogBeforeOp();
            
            newId = this.sr.MergeClusters(cid);
            this.currentCluster = newId;
            
            this.LogAfterOp();
            
            this.UpdateAppUI();
            this.PlotSpikes();
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
            
            this.LogBeforeOp();
            
            [newId, reId] = this.sr.CutClusters(spkInd);
            this.currentCluster = [newId reId];
            
            this.LogAfterOp();
            
            this.UpdateAppUI();
            this.PlotSpikes();
            this.vm.EnterSelectionMode();
        end
        
        function spkInd = FindSpikesInPolygon(this)
            % Return row indices in spkTb for spikes encircled by the polygon
            if isempty(this.currentCluster) || size(this.cutCoords, 1) < 3
                spkInd = [];
                return
            end
            m = ismember(this.sr.spkTb.clusId, this.currentCluster);
            spkInd = find(m);
            xq = this.sr.spkTb.timeSec(spkInd);
            yq = this.sr.spkTb.covCentCoords(spkInd,2);
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
        
        function LogBeforeOp(this)
            % Capture the current state of clusters, usually right before operation
            if ~height(this.logTb)
                this.LogAfterOp();
            else
                this.logTb.spkClusId{end} = this.sr.spkTb.clusId;
                this.logTb.currentCluster{end} = this.currentCluster;
                this.logTb.isDone(end) = true;
            end
        end
        
        function LogAfterOp(this)
            % Capture the new state of clusters, usually right after operation
            warning('off', 'MATLAB:table:RowsAddedExistingVars');
            k = height(this.logTb) + 1;
            this.logTb.spkClusId{k} = this.sr.spkTb.clusId;
            this.logTb.currentCluster{k} = this.currentCluster;
            this.logTb.isDone(k) = true;
        end
        
        function Undo(this)
            % Undo the last operation
            ind = find(this.logTb.isDone, 2, 'last');
            if numel(ind) >= 2
                this.logTb.isDone(ind(end)) = false;
                this.sr.spkTb.clusId = this.logTb.spkClusId{ind(end-1)};
                this.currentCluster = this.logTb.currentCluster{ind(end-1)};
                this.UpdateAppUI();
                this.PlotSpikes();
            end
        end
        
        function Redo(this)
            % Redo the last operation
            idx = find(~this.logTb.isDone, 1, 'first');
            if ~isempty(idx)
                this.logTb.isDone(idx) = true;
                this.sr.spkTb.clusId = this.logTb.spkClusId{idx};
                this.currentCluster = this.logTb.currentCluster{idx};
                this.UpdateAppUI();
                this.PlotSpikes();
            end
        end
        
        % Display
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
        
        function UpdateAppUI(this)
            % Update ClusTable in the Clustering tab
            
            if ~this.vm.hasApp
                return
            end
            app = this.vm.app;
            if ~this.hasClus
                app.ClusTable.Data = [];
                return
            end
            
            % Set the table
            tb = this.clusTb;
            app.ClusTable.Data = tb(:,this.appTbColNames);
            
            % Apply text colors
            removeStyle(app.ClusTable);
            for i = 1 : height(tb)
                s = uistyle("FontColor", tb.color(i,:));
                addStyle(app.ClusTable, s, 'row', i);
            end
            
            % Set selected row(s) and scroll to the first item
            ind = find(ismember(tb.clusId, this.currentCluster));
            app.ClusTable.Selection = ind(:)'; % must be a row vector
            if ~isempty(ind)
                scroll(app.ClusTable, 'row', ind(1));
            end
        end
        
        function FindOnMap(this, cid)
            % Zoom to the selected cluster(s) in the Maps Window
            
            if ~this.hasClus || ~this.vm.hasApp || isempty(cid)
                return
            end
            
            % Find the range of spike times
            m = ismember(this.sr.spkTb.clusId, cid);
            tSpk = this.sr.spkTb.timeSec(m);
            tWin = [min(tSpk), max(tSpk)];
            
            % Find the range of mean cluster depths
            m = ismember(this.sr.clusTb.clusId, cid);
            yClus = this.sr.clusTb.depth(m);
            yWin = [min(yClus), max(yClus)];
            
            % Set ROI
            this.vm.SetMapROI(tWin, yWin+[-500 500]);
        end
        
        function PlotSpikesBk(this)
            % Plot all or selected clusters
            
            if ~this.hasClus || ~this.vm.hasMapAxes
                return
            end
            
            tb = this.clusTb;
            
            if isempty(this.currentCluster)
                cid = tb.clusId;
            else
                cid = this.currentCluster;
            end
            
            % Plot the selected cluster later such that they are on top of the unselected ones
            ind = ismember(tb.clusId, cid);
            tb = [tb(~ind,:); tb(ind,:)];
            
            % Clear existing clusters in map
            delete(this.vm.mapLayers.clusters);
            
            % Plot clusters
            hh = cell(height(tb), 1);
            for i = 1 : height(tb)
                % Get cluster info
                k = tb.clusId(i);
                g = tb.group(i);
                c = tb.color(i,:);
                
                % Add alpha to unselected clusters, with "mua" lighter than "good"
                if strcmp(g, "good")
                    r = 0.4;
                else % "mua"
                    r = 0.2;
                end
                if ~ismember(k, cid)
                    w = ones(size(c));
                    c = r*c + (1-r)*w;
                end
                
                % Get spike data
                m = this.sr.spkTb.clusId == k;
                t = this.sr.spkTb.timeSec(m);
                y = this.sr.spkTb.covCentCoords(m,2);
                
                % Plot spikes
                h = line(this.vm.mapAxes, t, y, ...
                    'Color', c, ...
                    'LineStyle', 'none', ...
                    'Marker', '.', ...
                    'MarkerSize', this.vm.spikeMarkerSize, ...
                    'HitTest', 'off');
                
                % Label selected plots
                if ismember(k, cid)
                    h.UserData.id = k;
                    h.UserData.group = g;
                    h.UserData.color = c;
                end
                
                hh{i} = h;
            end
            this.vm.AddHandle2Layer(cat(1, hh{:}), 'clusters');
            this.vm.UpdateSpikeMarkerSize(1, 'relative');
            
            % Initialize brush callback
            brushObj = brush(this.vm.mapFig);
            brushObj.ActionPostCallback = @this.PlotBrushedSpikes;
        end
        
        function PlotSpikes(this)
            % Plot all or selected clusters
            
            if ~this.hasClus || ~this.vm.hasMapAxes
                return
            end
            
            this.UpdateClusHandleTb();
            tb = this.clusTb;
            hh = this.hClusTb.handle;
            
            if isempty(this.currentCluster)
                cid = tb.clusId;
            else
                cid = this.currentCluster;
            end
            
            % Plot clusters
            isSelected = ismember(tb.clusId, cid);
            for i = 1 : height(tb)
                % Get cluster info
                k = tb.clusId(i);
                g = tb.group(i);
                c = tb.color(i,:);
                
                % Add alpha to unselected clusters, with "mua" lighter than "good"
                if strcmp(g, "good")
                    r = 0.4;
                else % "mua"
                    r = 0.2;
                end
                if ~isSelected(i)
                    w = ones(size(c));
                    c = r*c + (1-r)*w;
                end
                
                % Make or update plot
                h = hh(i);
                if isempty(h) || ~isvalid(h)
                    % Get spike data
                    m = this.sr.spkTb.clusId == k;
                    t = this.sr.spkTb.timeSec(m);
                    y = this.sr.spkTb.covCentCoords(m,2);
                    
                    % Plot spikes
                    h = plot(this.vm.mapAxes, t, y, ...
                        'Color', c, ...
                        'LineStyle', 'none', ...
                        'Marker', '.', ...
                        'MarkerSize', this.vm.spikeMarkerSize, ...
                        'HitTest', 'off');
                else
                    if any(h.Color ~= c)
                        h.Color = c;
                    end
%                     if h.MarkerSize ~= this.vm.spikeMarkerSize
%                         h.MarkerSize = this.vm.spikeMarkerSize;
%                     end
                end
                
                % Label line object
                if isSelected(i)
                    h.UserData.id = k;
                    h.UserData.group = g;
                    h.UserData.color = c;
                else
                    h.UserData = [];
                end
                
                hh(i) = h;
            end
            
            uistack(hh(isSelected), 'top');
            this.hClusTb.handle = hh;
            this.vm.mapLayers.clusters = hh;
%             this.vm.AddHandle2Layer(hh, 'clusters');
%             this.vm.UpdateSpikeMarkerSize(1, 'relative');
            
            % Initialize brush callback
            brushObj = brush(this.vm.mapFig);
            brushObj.ActionPostCallback = @this.PlotBrushedSpikes;
        end
        
        function PlotBrushedSpikes(this, fig, axesStruct)
            % 
            
            % Find clusters and spikes to plot
            hh = axesStruct.Axes.Children;
            cid = NaN(numel(hh), 1);
            cc = NaN(numel(hh), 3);
            spkMask = cell(size(cid));
            for i = 1 : numel(hh)
                % Check if the handle is cluster spike plot
                s = hh(i).UserData;
                if ~isfield(s, 'id')
                    continue
                end
                
                % Check if selecting any data
                b = logical(hh(i).BrushData);
                if ~any(b)
                    continue
                end
                
                cid(i) = s.id;
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
            if isempty(this.waveFig) || ~isvalid(this.waveFig)
                this.waveFig = MPlot.Figure( ...
                    'Name', 'Waveform', ...
                    'NumberTitle', 'off', ...
                    'IntegerHandle', 'off', ...
                    'Interruptible', 'off', ...
                    'BusyAction', 'cancel');
            end
            
            figure(this.waveFig);
            hold off
            this.waveFig.Position(3) = 250*numel(cid);
            
            this.sr.PlotClusterWaveform(cid, 50, spkMask, 'NumChannels', 10, 'Color', cc, 'ShowCentroids', false);
            ax = MPlot.Axes(gca);
            ax.Title.String = ['Cluster ' num2str(cid(:)')];
        end
        
        function PlotWaveform(this, cid)
            
            if ~this.hasBin
                warning('Cannot plot waveform as the binary file is not available.');
                return
            end
            if nargin < 2
                cid = this.currentCluster;
            end
            if isempty(cid) || numel(cid) > this.maxClusPlot
                return
            end
            
            % Create figure if absent
            if isempty(this.waveFig) || ~isvalid(this.waveFig)
                this.waveFig = MPlot.Figure( ...
                    'Name', 'Waveform', ...
                    'NumberTitle', 'off', ...
                    'IntegerHandle', 'off', ...
                    'Interruptible', 'off', ...
                    'BusyAction', 'cancel');
            end
            
            figure(this.waveFig);
            hold off
            this.waveFig.Position(3) = 200*numel(cid);
            
            tb = this.clusTb;
            tb = tb(ismember(tb.clusId, cid),:);
            this.sr.PlotClusterWaveform(tb.clusId, 50, 'NumChannels', 10, 'Color', tb.color);
            ax = MPlot.Axes(gca);
            ax.Title.String = ['Cluster ' num2str(cid(:)')];
%             ax.LooseInset = [0 0 0 0];
        end
        
        function PlotCCG(this, cid)
            
            if nargin < 2
                cid = this.currentCluster;
            end
            if ~this.hasClus || isempty(cid) || numel(cid) > this.maxClusPlot
                return
            end
            
            % Create figure if absent
            if isempty(this.ccgFig) || ~isvalid(this.ccgFig)
                this.ccgFig = MPlot.Figure( ...
                    'Name', 'Cross-correlograms', ...
                    'NumberTitle', 'off', ...
                    'IntegerHandle', 'off', ...
                    'Interruptible', 'off', ...
                    'BusyAction', 'cancel');
            end
            
            figure(this.ccgFig);
            hold off
%             clf(this.ccgFig);
            
            tb = this.clusTb;
            tb = tb(ismember(tb.clusId, cid),:);
            this.sr.PlotCCG(tb.clusId, 'Color', tb.color, 'ShowMetrics', true);
        end
        
        function PlotAmplitude(this, cid)
            
            if nargin < 2
                cid = this.currentCluster;
            end
            if ~this.hasClus || isempty(cid) || numel(cid) > this.maxClusPlot
                return
            end
            
            % Create figure if absent
            if isempty(this.ampFig) || ~isvalid(this.ampFig)
                this.ampFig = MPlot.Figure( ...
                    'Name', 'Amplitude', ...
                    'NumberTitle', 'off', ...
                    'IntegerHandle', 'off', ...
                    'Interruptible', 'off', ...
                    'BusyAction', 'cancel');
            end
            
            figure(this.ampFig);
            hold off
            
            tb = this.clusTb;
            tb = tb(ismember(tb.clusId, cid),:);
            this.sr.PlotAmplitude(tb.clusId, 'Color', tb.color);
%             ax.LooseInset = [0 0 0 0];
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
