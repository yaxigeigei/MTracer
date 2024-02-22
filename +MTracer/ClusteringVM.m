classdef ClusteringVM < handle
    % 
    
    properties(Constant)
        appTbColNames = {'clusId', 'group', 'depth', 'numSpikes', 'RPV', 'contam', 'SNR'}
        noiseColor = [0 0 0] + .5
        maxClusOp = 8
%         featList = {'Amplitude', 'X-Position', 'PC1'};
        featList = {'Amplitude', 'X-Position'};
    end
    
    properties
        vm MTracerVM                % handle of the main app view-model
        sr                          % object managing sorting results
        om MTracer.OperationManager % object managing clustering operation history
        plg MTracer.Polygon         % object managing spikes selection
        coi = []                    % clusters of interest, usually aka the selected clusters
        
        waveFig                     % object of the waveform window
        ccgFig                      % object of the CCG window
        ftFig                       % a vector of objects for feature-time windows
    end
    
    properties(Dependent)
        hasClus
        hasBin
        isShowNoise
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
            this.plg = MTracer.Polygon(vm);
            
            this.waveFig = MTracer.LayeredFigure();
            
            this.ccgFig = MTracer.LayeredFigure();
            
            this.ftFig = MTracer.LayeredFigure();
            this.ftFig.name = this.featList{1};
            
            if nargin > 1
                this.sr = sr;
                this.om = MTracer.OperationManager(this.GetStateData([], []));
                this.SetClusterColor();
            end
        end
        
        function delete(this)
            delete(this.waveFig);
            delete(this.ccgFig);
            delete(this.ftFig);
        end
        
        function cvm = Duplicate(this, appvm)
            % Make a hard copy of the object
            
            cvm = MTracer.ClusteringVM(appvm);
            
            if isempty(this.sr)
                return
            end
            
            % Make a hard copy
            srCopy = this.sr.Duplicate();
            
            % Remove added columns
            isAdded = ismember(srCopy.clusTb.Properties.VariableNames, {'handle', 'memmap'});
            srCopy.clusTb(:,isAdded) = [];
            
            cvm.sr = srCopy;
            cvm.om = this.om.Duplicate;
            cvm.coi = this.coi;
            cvm.waveFig = this.waveFig.Duplicate;
            cvm.ccgFig = this.ccgFig.Duplicate;
            cvm.ftFig = this.ftFig.Duplicate;
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
            this.sr.ExportData(this.vm.ksFolder);
        end
        
        % Cluster operations
        function s = GetStateData(this, preCOI, postCOI)
            % Assemble a struct of state data for this.om to manage
            s.spkClusId = this.sr.spkTb.clusId;
            s.preClusId = preCOI(:)';
            s.postClusId = postCOI(:)';
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
        
        function SetClusterColor(this)
            % Assign color to clusters
            %   1) Colors are generated by lines()
            %   2) To initialize, colors are assigned to clusters in the order of descending depth 
            %      such that clusters close in space do not have the same color.
            %   3) Black [0 0 0] is reserved for unassigned clusters, and will be assigned a color 
            %      upon calling this method.
            %   4) When assigning a color to a cluster, the color will be different from other 6 
            %      closest non-noise neighboring clusters.
            %   5) 'noise' clusters, though having colors, are always plotted in gray (this.noiseColor).
            
            tb = this.sr.clusTb;
            
            % Initialize colors if empty
            if ~ismember('color', tb.Properties.VariableNames)
                % Find the forward and backward sorting indices
                [~, fwInd] = sort(tb.depth, 'descend');
                [~, bwInd] = sort(fwInd);
                cc = lines(height(tb));
                tb.color = cc(bwInd,:);
            end
            
            % Find unassigned clusters (i.e. [0 0 0])
            indNoColor = find(~any(tb.color, 2));
            
            % Assign colors
            cMap = lines(7);
            isExclude = tb.group == "noise" | tb.numSpikes == 0;
            for i = indNoColor(:)'
                % Find the colors of the 6 closest clusters
                d = abs(tb.depth - tb.depth(i));
                d(i) = Inf;
                d(isExclude) = Inf;
                [~, I] = sort(d);
                cUsed = tb.color(I(1:min(6,end)),:);
                
                % Find and assign the first unused color
                isUsed = ismember(cMap, cUsed, 'rows');
                k = find(~isUsed, 1);
                tb.color(i,:) = cMap(k,:);
            end
            
            this.sr.clusTb = tb;
        end
        
        function RecolorClusters(this)
            % Reassign color to clusters
            
            cid = this.coi;
            if ~this.hasClus || numel(cid) < 1
                return
            end
            
            % Unassign colors
            tb = this.sr.clusTb;
            m = ismember(tb.clusId, cid);
            tb.color(m,:) = 0;
            this.sr.clusTb = tb;
            
            % Assign colors
            this.SetClusterColor();
            
            % Refresh
            this.UpdateAll();
        end
        
        function LabelClusters(this, newGroup)
            % Change the group label of specified clusters
            cid = this.coi;
            if ~this.hasClus || numel(cid) < 1
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
            this.plg.ClearPolygon();
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
            
            if isempty(this.coi) || this.plg.numPoints < 3
                spkInd = [];
                return
            end
            m = ismember(this.sr.spkTb.clusId, this.coi);
            spkInd = find(m);
            
            % Get feature points
            xq = this.sr.spkTb.timeSec(spkInd);
            switch this.plg.target
                case 'X-Position'
                    yq = this.sr.spkTb.centCoords(spkInd,1);
                case 'Y-Position'
                    yq = this.sr.spkTb.centCoords(spkInd,2);
                case 'Amplitude'
                    yq = this.sr.spkTb.tempAmp(spkInd);
                otherwise
                    error("'%s' is not a valid object name", this.plg.target);
            end
            
            % Find enclosed points
            in = this.plg.IsInPolygon(xq, yq);
            spkInd = spkInd(in);
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
            if app.AutoFtCheckBox.Value
                this.PlotFeatureTime();
            end
        end
        
        function UpdateAll(this)
            this.UpdateCOI();
            this.UpdateAppClusTable();
        end
        
        function FindOnMap(this, cid)
            % Zoom to the selected cluster(s) in the Maps Window
            
            if ~this.hasClus || ~this.vm.isOpen
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
            yWin = [min(yClus) max(yClus)] + [-1 1]*250;
            
            % Set ROI
            this.vm.SetMapROI(this.vm.hAxes.XLim, yWin);
        end
        
        function FindInTable(this)
            % Select the cluster that's closest in depth to the focus
            if ~this.hasClus || ~this.vm.hasApp || ~this.vm.isOpen
                return
            end
            y = this.vm.app.ClusTable.Data.depth;
            y0 = this.vm.focus(2);
            [~, I] = sort(abs(y - y0));
            cid = this.vm.app.ClusTable.Data.clusId(I(1));
            this.coi = cid;
            this.UpdateAll();
        end
        
        function NextFeature(this)
            % Cycle to the next feature
            featName = this.ftFig.name;
            I = find(strcmp(featName, this.featList));
            if I == numel(this.featList)
                I = 0;
            end
            this.ftFig.name = this.featList{I+1};
            this.PlotFeatureTime();
        end
        
        function SetPoint(this, src, eventdata)
            % Add or remove polygon vertex from the caller axes
            
            if ~strcmp(this.vm.mapMode, 'clustering')
                return
            end
            
            % Read mouse position
            mousePos = get(src, 'CurrentPoint');
            t = mousePos(1);
            y = mousePos(3);
            
            if eventdata.Button == 1
                this.plg.AddPoint(src, t, y);
            elseif eventdata.Button == 3
                this.plg.RemoveLastPoint();
            end
        end
        
        % Plotting
        function PlotSpikes(this)
            % Plot all or selected clusters
            
            if ~this.hasClus || ~this.vm.isOpen
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
                
                % Make or update plot
                if isempty(properties(h)) || ~isvalid(h) % placeholder does not have properties
                    % Get spike data
                    m = this.sr.spkTb.clusId == k;
                    t = this.sr.spkTb.timeSec(m);
                    y = this.sr.spkTb.centCoords(m,2);
                    
                    % Plot spikes
                    h = plot(this.vm.hAxes, t, y, ...
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
            this.vm.layers.clusters = hh(isvalid(hh));
            
            % Initialize brush callback
            brushObj = brush(this.vm.hFig);
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
            
            % Make plot
            this.PlotWaveform(cid, spkMask);
        end
        
        function PlotFeatureTime(this, cid, featName)
            
            if nargin < 3
                featName = this.ftFig.name;
            end
            if nargin < 2
                cid = this.coi;
            end
            if ~this.hasClus || isempty(cid) || numel(cid) > this.maxClusOp
                return
            end
            
            % Get figure and axes
            if ~this.ftFig.isOpen
                % Create figure and axes if absent
                [f, ax] = this.ftFig.Open( ...
                    'Name', 'Feature-Time', ...
                    'NumberTitle', 'off', ...
                    'IntegerHandle', 'off', ...
                    'Menubar', 'none', ...
                    'Toolbar', 'figure', ...
                    'WindowKeyPressFcn', @(varargin) this.vm.KeyPress(varargin{:}), ...
                    'WindowKeyReleaseFcn', @(varargin) this.vm.KeyRelease(varargin{:}), ...
                    'WindowScrollWheelFcn', @(varargin) this.vm.Scroll(varargin{:}), ...
                    'Interruptible', 'off', ...
                    'BusyAction', 'cancel');
                
                % Align new figure to the main map
                f.Position(4) = 150; % height
                if this.vm.isOpen
                    inPos = this.vm.hFig.Position;
                    outPos = this.vm.hFig.OuterPosition;
                    f.Position([1 3]) = inPos([1 3]); % origin x, width
                    f.Position(2) = sum(outPos([2 4])); % origin y
                    
%                     mainPos = this.vm.hAxes.Position;
%                     ax.Position([1 3]) = mainPos([1 3]);
                end
            else
                ax = this.ftFig.hAxes;
            end
            hold(ax, 'off');
            
            % Get feature data
            switch featName
                case 'Amplitude'
                    feat = this.sr.spkTb.tempAmp;
                    yLabel = featName + " (au)";
                case 'X-Position'
                    feat = this.sr.spkTb.centCoords(:,1);
                    yLabel = 'X-Position (um)';
%                 case 'PC1'
%                     feat = cellfun(@(x) x(16,1), this.sr.spkTb.pcWeights);
%                     yLabel = featName + " (au)";
                otherwise
                    error("'%s' is not a valid feature name", featName);
            end
            
            % Plot through clusters
            cTb = this.sr.clusTb;
            hh = cell(size(cid));
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
                y = feat(m);
                
                % Plot spikes
                h = plot(ax, t, y, ...
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
                hh{i} = h;
            end
            
            % Save handles
            this.ftFig.layers.clusters = [hh{:}];
            
            % Tag axes by the feature name
            ax.Tag = featName;
            
            MPlot.Axes(ax);
            if this.vm.isOpen
                ax.XLim = this.vm.hAxes.XLim;
            end
            ax.TickLength(1) = 0.002;
            ax.XTickLabels = [];
            ax.YLabel.String = yLabel;
            ax.LooseInset = [0 0 0 0];
            
            % Mouse button callback
            ax.ButtonDownFcn = @this.SetPoint;
            ax.BusyAction = 'cancel';
            
            % Initialize brush callback
            brushObj = brush(this.ftFig.hFig);
            brushObj.ActionPostCallback = @this.PlotBrushedSpikes;
        end
        
        function UpdateSpikeMarkerSize(this, sz)
            % Update spike marker size in feature-time figures
            if this.ftFig.isOpen
                hh = this.ftFig.layers.clusters;
                for i = 1 : numel(hh)
                    hh(i).MarkerSize = sz;
                end
            end
        end
        
        function PlotWaveform(this, cid, spkMask)
            
            if ~this.hasBin
                fprintf('Cannot plot waveform since the binary file is not available.\n');
                return
            end
            if nargin < 3
                spkMask = [];
            end
            if nargin < 2
                cid = this.coi;
            end
            if isempty(cid) || numel(cid) > this.maxClusOp
                return
            end
            
            [f, ax] = this.waveFig.Open( ...
                'Name', 'Waveform', ...
                'NumberTitle', 'off', ...
                'Menubar', 'none', ...
                'Toolbar', 'none', ...
                'IntegerHandle', 'off', ...
                'Interruptible', 'off', ...
                'BusyAction', 'cancel');
            
            hold(ax, 'off');
            this.sr.PlotClusterWaveform(cid, 30, spkMask, 'NumChannels', 16, 'Color', this.GetClusterColor(cid), 'Axes', ax);
%             this.sr.PlotClusterTemplate(cid, 'Axes', ax);
            ax.Title.String = ['Cluster ' num2str(cid(:)')];
            ax.LooseInset([2 4]) = [0 0];
            f.Position(3) = diff(ax.XLim) * 30 + 15;
        end
        
        function PlotCCG(this, cid)
            
            if nargin < 2
                cid = this.coi;
            end
            if ~this.hasClus || isempty(cid) || numel(cid) > this.maxClusOp
                return
            end
            
            f = this.ccgFig.Open( ...
                'Name', 'Cross-correlograms', ...
                'NumberTitle', 'off', ...
                'Menubar', 'none', ...
                'Toolbar', 'none', ...
                'IntegerHandle', 'off', ...
                'Interruptible', 'off', ...
                'BusyAction', 'cancel');
            
            this.sr.PlotCCG(cid, 'Color', this.GetClusterColor(cid), 'ShowMetrics', true, 'Figure', f);
        end
        
    end
    
end
