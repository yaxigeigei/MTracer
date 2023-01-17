classdef SortingResult < handle
    %SortingResult Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        chanTb table;
        spkTb table;
        clusTb table;
        tempTb table;
        mdat;
        samplingRate = 30e3;
        sampleOffset = 0;
        pcBases;
    end
    
    methods
        % Constructor
        function this = SortingResult()
            % do nothing
        end
        
        function sr = Duplicate(this)
            % Make a hard copy of the object
            sr = NP.SortingResult();
            p = properties(this);
            for i = 1 : numel(p)
                sr.(p{i}) = this.(p{i});
            end
        end
        
        % Data IO
        function [W, chWins, tmWins] = ExtractWaveform(this, spkInd, varargin)
            % Extract waveform for the given spike indices (not spike time indices)
            
            p = inputParser();
            p.addParameter('NumSamples', 82, @(x) isscalar(x) && isnumeric(x));
            p.addParameter('NumChannels', 32, @(x) isscalar(x) && isnumeric(x));
            p.addParameter('ChunkSize', 1e4, @(x) isscalar(x) && isnumeric(x));
            p.parse(varargin{:});
            nSp = p.Results.NumSamples;
            nCh = p.Results.NumChannels;
            tmWin = [-floor(nSp/2), ceil(nSp/2)-1];
            chWin = [-floor(nCh/2), ceil(nCh/2)-1];
            chunkSz = p.Results.ChunkSize;
            
            % Find spike info
            if isempty(spkInd)
                spkInd = (1 : height(this.spkTb))';
            elseif islogical(spkInd)
                spkInd = find(spkInd);
            end
            spkTime = this.spkTb.timeInd(spkInd) - this.sampleOffset;
            spkTpInd = this.spkTb.tempId(spkInd) + 1;
            spkTpCent = this.tempTb.centIdx(spkTpInd);
            
            % Chunking
            nSpk = numel(spkInd);
            nChunk = ceil(nSpk / chunkSz);
            W = cell(nChunk,1);
            
            % Extract waveform by chunks to reduce memory usage
            for i = 1 : nChunk
                a = (i-1)*chunkSz + 1;
                b = min(i*chunkSz, nSpk);
                W{i} = MKilosort2.ReadSnippets(this.mdat, spkTime(a:b), tmWin, spkTpCent(a:b), chWin, ...
                    'ChannelOrder', this.chanTb.sortInd);
                disp([num2str(b) ' waveform extracted']);
            end
            W = cat(3, W{:});
            
            % Return the windows used
            tmWins = spkTime + tmWin;
            chWins = spkTpCent + chWin;
        end
        
        function sn = ReadSnippets(this, spInd, spWin, varargin)
            % A wrapper of MKilosort2.ReadSnippets with preset parameters
            spInd = spInd - this.sampleOffset;
            sn = MKilosort2.ReadSnippets(this.mdat, spInd, spWin, varargin{:}, ...
                    'ChannelOrder', this.chanTb.sortInd);
        end
        
        function WriteSpikeAudio(this, fileName, cid, tWin)
            % Extract spike waveform and plot 
            
            % Audio specs
            Fs = 44100;
            tAud = tWin(1) : 1/Fs : tWin(2);
            
            % Find spikes that belong to the cluster and are within the time window
            isClus = this.spkTb.clusId == cid;
            tSpk = this.spkTb.timeSec;
            spkInd = isClus & tSpk > tWin(1) & tSpk < tWin(2);
            
            % Extract single-channel waveform
            [wf, ~, wfWin] = this.ExtractWaveform(spkInd, 'NumChannels', 1);
            
            if ~isempty(wf)
                % Find spike waveform timestamps for audio
                wf = permute(wf, [2 3 1]);
                tWf = zeros(size(wf));
                for i = 1 : size(wf,2)
                    tWf(:,i) = wfWin(i,1) : wfWin(i,2);
                end
                wf = wf(:);
                tWf = tWf(:) / this.samplingRate;
                
                % Remove overlapping samples
                [tWf, ind] = unique(tWf);
                wf = wf(ind);
                
                % Generate full audio by interpolation
                wAud = interp1(tWf, double(wf), tAud, 'linear', 0);
                wAud = wAud / max(abs(wAud)) * 0.1; % adjust amplitude, 0.1 is just an emperical factor
            else
                warning('Cluster %i has no spikes in the specified time window. The audio will be silent.', cid);
                wAud = zeros(size(tAud));
            end
            
            % Write audio
            audiowrite(fileName, wAud, Fs);
        end
        
        % Computing
        function ComputeAll(this)
            % Compute all metrics for all clusters
            
            % Compute cluster spiking metrics
            disp('Computing spiking metrics');
            this.ComputeContamStats([]);
            
%             % Compute wavefrom centeroids
%             disp('Localizing spikes');
%             this.ComputeWaveformCenter([]);
            
%             % Extract and cache center waveform in batches
%             nSpk = height(this.spkTb);
%             batchSize = 1e5; % number of spikes
%             nBatches = ceil(nSpk/batchSize);
%             centW = cell(nBatches, 1);
%             for i = 1 : nBatches
%                 fprintf('Batch %i/%i\n', i, nBatches);
%                 spkInd = (i-1)*batchSize+1 : min(i*batchSize, nSpk);
%                 centW{i} = this.ExtractWaveform(spkInd, 'NumChannels', 10);
%             end
%             centW = cat(3, centW{:});
            
            % Compute cluster waveform metrics
            disp('Computing waveform metrics');
%             [centW, chWin, tmWin] = this.ExtractWaveform([], 'NumChannels', 10);
%             this.ComputeWaveformStats([], centW, chWin);
            this.ComputeWaveformStats();
            this.ComputeClusterDepth();
        end
        
        function spkXY = ComputeWaveformCenter(this, spkInd)
            % Compute centroids of spike waveform
            % 
            %	ComputeWaveformCenter()
            %	ComputeWaveformCenter(spkInd)
            % 
            
            if nargin < 2 || isempty(spkInd)
                spkInd = (1:height(this.spkTb))';
            end
            
            sTb = this.spkTb(spkInd,:);
            tTb = this.tempTb;
            
            % Preallocate center coordinates
            if ~ismember('centCoords', sTb.Properties.VariableNames)
                sTb.centCoords = NaN(height(sTb), 2);
            end
            
            % Loop through spikes by templates
            tidList = unique(sTb.tempId);
            for tid = tidList(:)'
                % Check if template is present
                isTemp = tTb.tempId == tid;
                if isempty(tTb.temp{isTemp})
                    continue
                end
                
                % Find spikes that use this template
                isSpk = sTb.tempId == tid;
%                 if ~any(isSpk)
%                     continue
%                 end
                
                % Read waveform
                [W, chWin] = this.ExtractWaveform(isSpk, 'NumChannels', 32);
                
                % Get channel coordinates
                chXY = [tTb.chanX{isTemp} tTb.chanY{isTemp}];
                
                % Estimate centroids
                spkXY = MNeuro.ComputeWaveformCenter(W, chXY, 'power', 'centroid');
                
                % Save to spike table
                sTb.centCoords(isSpk,:) = spkXY;
            end
            spkXY = sTb.centCoords;
            this.spkTb.centCoords(spkInd,:) = spkXY;
        end
        
        function spkXY = ComputeWaveformCenterRecon(this, spkInd)
            % Compute centroids from denoised (reverse PCA) spike waveform
            % 
            %	ComputeWaveformCenter()
            %	ComputeWaveformCenter(spkInd)
            % 
            
            if nargin < 2 || isempty(spkInd)
                spkInd = (1:height(this.spkTb))';
            end
            
            sTb = this.spkTb(spkInd,:);
            tTb = this.tempTb;
            
            % Get principal component bases
            B = single(this.pcBases);
            
            % Get projection weights
            P = cat(3, sTb.pcWeights{:}); % concatenate along spikes
            P = permute(P, [2 1 3]); % to pc-by-channel-by-spike
            
            % Preallocate center coordinates
            if ~ismember('centCoords', sTb.Properties.VariableNames)
                sTb.centCoords = NaN(height(sTb), 2);
            end
            
            % Loop through spikes by templates
            tidList = unique(sTb.tempId);
            for tid = tidList(:)'
                % Check template
                isTemp = tTb.tempId == tid;
                if isempty(tTb.temp{isTemp})
                    continue
                end
                
                % Find spikes of the template
                isSpk = sTb.tempId == tid;
                if ~any(isSpk)
                    continue
                end
                
                % Reconstruct waveform from bases and projection weights
                clear W
                spkP = P(:,:,isSpk);
                for k = size(spkP,3) : -1 : 1
                    W(:,:,k) = B * spkP(:,:,k);
                end
%                 W = pagemtimes(B, P(:,:,isSpk));
                W = permute(W, [2 1 3]); % to channel-by-time-by-spike
                
                % Find channel coordinates
                chXY = [tTb.chanX{isTemp} tTb.chanY{isTemp}];
                
                % Estimate centroids
                spkXY = MNeuro.ComputeWaveformCenter(W, chXY, 'power', 'centroid');
                
                % Save to spike table
                sTb.centCoords(isSpk,:) = spkXY;
            end
            spkXY = sTb.centCoords;
            this.spkTb.centCoords(spkInd,:) = spkXY;
        end
        
        function ComputeNumSpikes(this, cid)
            % Compute the number of spikes for specified clusters. Results will be added/updated in clusTb.
            %
            %   ComputeNumSpikes()
            %	ComputeNumSpikes(cid)
            % 
            % Inputs
            %   cid         A vector of cluster IDs. If not provided or an empty array, the methods will 
            %               compute for all clusters.
            
            tb = this.clusTb;
            
            if nargin < 2 || isempty(cid)
                cid = tb.clusId;
            end
            rowInd = find(ismember(tb.clusId, cid));
            
            for i = rowInd(:)'
                m = this.spkTb.clusId == tb.clusId(i);
                tb.numSpikes(i) = sum(m);
            end
            
            this.clusTb = tb;
        end
        
        function ComputeContamStats(this, cid)
            % Compute contamination metrics such as refractory period violation rate (RPV) and contamination rate 
            % for specified clusters. Results will be added/updated in clusTb.
            %
            %   ComputeContamStats()
            %	ComputeContamStats(cid)
            % 
            % Inputs
            %   cid         A vector of cluster IDs. If not provided or an empty array, the methods will 
            %               compute for all clusters.
            % 
            % Output
            %   The following variables will be added/updated in clusTb for each cluster
            %   numSpikes           The number of spikes.
            %   isiEdges            Bin edges of inter-spike interval (ISI) histogram in seconds.
            %   isiCount            Spike count of ISI histogram.
            %   RPV                 Refractory period violation rate (%). The RP threshold is NP.Param.minISI.
            %   meanActiveRate      Mean spike rate when unit is active.
            %   activeThreshold     The threshold for active firing, computed from NP.Param.activePrct.
            %   contam              Contamination rate (%) (Hill, Mehta and Kleinfeld, J Neuro, 2011).
            %
            
            tb = this.clusTb;
            
            if nargin < 2 || isempty(cid)
                cid = tb.clusId;
            end
            rowInd = find(ismember(tb.clusId, cid));
            
            for i = rowInd(:)'
                % Find spikes of this cluster
                m = this.spkTb.clusId == tb.clusId(i);
                
                % Compute basic stats
                tb.numSpikes(i) = sum(m);
                
                % Compute comtamination
                tSpk = this.spkTb.timeInd(m) / 30e3;
                s = NP.Unit.ComputeContamStats(tSpk);
                
                % Add results to table
                fns = fieldnames(s);
                for k = 1 : numel(fns)
                    fn = fns{k};
                    if isscalar(s.(fn))
                        tb.(fn)(i) = s.(fn);
                    else
                        tb.(fn){i} = s.(fn);
                    end
                end
            end
            
            this.clusTb = tb;
        end
        
        function ComputeWaveformStats(this, cid, wf, chWins)
            % Compute waveform location and statistics such as median, SD, and SNR for specified clusters. 
            % Results will be added/updated in clusTb.
            % 
            %   ComputeWaveformStats()
            %	ComputeWaveformStats(cid)
            %	ComputeWaveformStats(cid, W, chWins)
            % 
            % Inputs
            %   cid         A vector of cluster IDs. If not provided or an empty array, the methods will 
            %               compute for all clusters.
            %   wf          A #channels-by-#timepoints-by-#spikes array for all spikes of interest. 
            %               If not provided, the method will extract them from the binary file.
            % 
            % Output
            %   The following variables will be added/updated in clusTb for each cluster
            %   numSpikes       The number of spikes.
            %   waveformMed     Median waveform in #channels-by-#timepoints matrix.
            %   waveformSD      Standard deviation of the waveform in #channels-by-#timepoints matrix. 
            %                   For robustness, this is computed using MAD/0.6745.
            %   SNR             Signal-to-noise ratio of the waveform.
            %   depth           If 'centCoords' is available.
            % 
            
            cTb = this.clusTb;
            if ~exist('cid', 'var') || isempty(cid)
                cid = cTb.clusId;
            end
            
            % Find spikes of the selected clusters
            isClus = ismember(this.spkTb.clusId, cid);
            sTb = this.spkTb(isClus,:);
            
            % Only keep waveform inputs for the selected clusters
            if exist('wf', 'var') && size(wf,3) == height(this.spkTb)
                wf = wf(:,:,isClus);
            end
            if exist('chWins', 'var') && size(chWins,1) == height(this.spkTb)
                chWins = chWins(isClus,:);
            end
            
            for i = 1 : numel(cid)
                % Find spikes of the current cluster
                m = sTb.clusId == cid(i);
                
                % Compute basic cluster stats
                cTb.numSpikes(i) = sum(m);
                if ~cTb.numSpikes(i)
                    warning('Cluster %i has no spikes', cid(i));
                    continue
                end
                
                % Get waveform
                if ~exist('wf', 'var') || isempty(wf)
                    [w, cw] = this.ExtractWaveform(this.spkTb.clusId == cid(i), 'NumChannels', 10);
                elseif isnumeric(W) && size(W,3) == numel(m)
                    w = wf(:,:,m);
                    cw = chWins(m,:);
                end
                w = double(w);
                
                % Average waveform
                k = cTb.clusId == cid(i);
                [wMed, ~, wMad] = MMath.MedianStats(w, 3);
                cTb.waveformMed{k} = wMed;
                cTb.waveformSD{k} = wMad / 0.6745;
                
                % Compute SNR
                ic = ceil(size(w,1)/2); % find the middle channel
                r = w(ic,:,:) - wMed(ic,:); % 1-by-nTm-by-nW residuals
                [~, ~, rMad] = MMath.MedianStats(r(:));
                cTb.SNR(k) = (max(wMed(ic,:)) - min(wMed(ic,:))) / (rMad/0.6745);
                
                % Compute waveform centroid
                chInd = arrayfun(@(x,y) x:y, cw(:,1), cw(:,2), 'Uni', false);
                chInd = cat(1, chInd{:})';
                chXY = cat(3, this.chanTb.xcoords(chInd), this.chanTb.ycoords(chInd));
                chXY = permute(chXY, [1 3 2]);
                spkXY = MNeuro.ComputeWaveformCenter(w, chXY, 'power', 'centroid');
                sTb.centCoords(m,:) = spkXY;
                cTb.depth(k) = round(median(sTb.centCoords(m,2), 'omitnan'));
            end
            
            % Update cluster table
            this.clusTb = cTb;
            
            % Update centroid coordinates to spike table
            if ~ismember('centCoords', this.spkTb.Properties.VariableNames)
                this.spkTb.centCoords = NaN(height(this.spkTb), 2);
            end
            this.spkTb.centCoords(isClus,:) = sTb.centCoords;
        end
        
        function ComputeClusterDepth(this, cid)
            % Compute waveform statistics such as median, SD, and SNR for specified clusters. 
            % Results will be added/updated in clusTb.
            %
            %   ComputeClusterDepth()
            %	ComputeClusterDepth(cid)
            % 
            % Inputs
            %   cid         A vector of cluster IDs. If not provided or an empty array, the methods will 
            %               compute for all clusters.
            % Outputs
            %   The following variables will be added/updated in clusTb for each cluster
            %   numSpikes       The number of spikes.
            %   depth           The median spike depth.
            
            if ~ismember('centCoords', this.spkTb.Properties.VariableNames)
                error('Spikes have not been localized. Run ComputeWaveformCenter method first.');
            end
            
            cTb = this.clusTb;
            if ~exist('cid', 'var') || isempty(cid)
                cid = cTb.clusId;
            end
            
            for i = 1 : numel(cid)
                % Find spikes of the current cluster
                m = this.spkTb.clusId == cid(i);
                k = cTb.clusId == cid(i);
                
                % Compute basic cluster stats
                cTb.numSpikes(k) = sum(m);
                
                % Compute median cluster depth
                if cTb.numSpikes(k)
                    cTb.depth(k) = round(median(this.spkTb.centCoords(m,2), 'omitnan'));
                else
                    warning('Cluster #%i has no spikes.', cid(i));
                end
            end
            
            this.clusTb = cTb;
        end
        
        function tps = IVectorizeTemplates(this, chWin)
            % Crop and vectorize templates
            if ~exist('chWin', 'var')
                chWin = -4 : 5;
            end
            tps = this.tempTb.temp;
            for i = 1 : numel(tps)
                tp = tps{i};
                
                if all(~tp, 'all')
                    warning('Template #%i (zero-based indexing) is empty. %i spikes are associated with this template.', ...
                        i-1, sum(this.spkTb.tempId==i-1));
                    tps{i} = [];
                    continue
                end
                
                % Crop template to a few channels in the middle
                chInd = this.tempTb.centIdx(i) + chWin;
                tp = tp(chInd,:);
                
                % Vectorize template by concatenating channels along time
                tp = tp';
                tp = tp(:)'; % turn template into a row vector by concatenating channels
                tps{i} = double(tp);
            end
        end
        
        % Cluster operations
        function LabelClusters(this, cid, newGroup)
            % Change the group label of specified clusters
            newGroup = cellstr(newGroup);
            if numel(newGroup) == 1 && numel(cid) > 1
                newGroup = repmat(newGroup, size(cid));
            end
            for i = 1 : numel(cid)
                this.clusTb.group{this.clusTb.clusId==cid(i)} = newGroup{i};
            end
        end
        
        function newId = MergeClusters(this, oldIds, newGroup)
            % Merge two or more clusters to one
            % 
            %   newId = MergeClusters(oldIds)
            %   newId = MergeClusters(oldIds, newGroup)
            % 
            % Inputs
            %   oldIds              A vector of cluster IDs, or a cell array of such vectors.
            %   newGroup            'auto' (default), 'good', 'mua', or 'noise'
            % Output
            %   newId               The ID of the merged cluster (max ID plus one).
            % 
            
            if ~exist('newGroup', 'var') || isempty(newGroup)
                newGroup = 'auto';
            end
            
            if iscell(oldIds)
                % Merge each set of clusters respectively with recursion
                if ischar(newGroup)
                    newGroup = repmat({newGroup}, size(oldIds));
                end
                for i = 1 : numel(oldIds)
                    this.MergeClusters(oldIds{i}, newGroup{i});
                end
                return
            end
            
            if isempty(oldIds)
                newId = [];
                return
            end
            
            % Determine new cluster ID
            newId = max(this.clusTb.clusId) + 1;
            
            % Determine new group label
            isOldClus = ismember(this.clusTb.clusId, oldIds);
            if strcmp(newGroup, 'auto')
                groupTypes = unique(this.clusTb.group(isOldClus));
                if numel(groupTypes) == 1
                    newGroup = groupTypes{1};
                else
                    newGroup = 'mua';
                end
            end
            
            % Update spike cluster ID
            isOldSpk = ismember(this.spkTb.clusId, oldIds);
            this.spkTb.clusId(isOldSpk) = newId;
            
            % Append new cluster to clusTb
            disp(['Merge #' num2str(oldIds(:)') ' to #' num2str(newId) ' (' newGroup ')']);
            warning('off', 'MATLAB:table:RowsAddedExistingVars');
            this.clusTb.clusId(end+1) = newId;
            this.clusTb.group{end} = newGroup;
            
            % Compute cluster metrics
            disp('Recomputing cluster metrics');
            this.clusTb.numSpikes(isOldClus) = 0; % clear spike number of the old clusters
            this.ComputeContamStats(newId);
            this.ComputeClusterDepth(newId);
        end
        
        function [newId, reId, oldId] = CutClusters(this, spkInd, newGroup)
            % Cut spikes out as a new cluster
            %   The cluster ID of the selected spikes will be the current max ID plus one.
            %   The ID of the clusters being cut will be updated following the new cluster ID.
            % 
            %   [newId, reId, oldId] = CutClusters(spkInd, newGroup)
            % 
            % Inputs
            %   spkInd              A vector of spike indices, i.e. row indices of spkTb
            %   newGroup            'good', 'mua', or 'noise'
            % Outputs
            %   newId               The cluster ID of the cut-out spikes.
            %   reId                Renamed ID(s) of the cluster(s) being cut.
            %   oldId               Original ID(s) of the cluster(s) being cut.
            
            if isempty(spkInd) || ~any(spkInd)
                newId = [];
                reId = [];
                return
            end
            
            % Find clusters involved in the cut
            oldId = unique(this.spkTb.clusId(spkInd)'); % needs row vector
            isOldClus = ismember(this.clusTb.clusId, oldId);
            
            % Determine new group label
            if ~exist('newGroup', 'var') || isempty(newGroup)
                newGroup = 'auto';
            end
            if strcmp(newGroup, 'auto')
                groupTypes = unique(this.clusTb.group(isOldClus));
                if numel(groupTypes) == 1
                    newGroup = groupTypes{1};
                else
                    newGroup = 'mua';
                end
            end
            
            % Create new cluster from the cut
            warning('off', 'MATLAB:table:RowsAddedExistingVars');
            newId = max(this.clusTb.clusId) + 1;
            disp(['Cut out a new cluster #' num2str(newId) ' (' newGroup ')']);
            this.spkTb.clusId(spkInd) = newId;
            this.clusTb.clusId(end+1) = newId;
            this.clusTb.group{end} = newGroup;
            
            % Rename remaining clusters
            reId = max(this.clusTb.clusId) + (1:numel(oldId)); % needs row vector
            disp(['Rename old clusters #' num2str(oldId) ' to #' num2str(reId) ', respectively']);
            for k = 1 : numel(oldId)
                this.spkTb.clusId(this.spkTb.clusId==oldId(k)) = reId(k);
                this.clusTb.clusId(end+1) = reId(k);
                this.clusTb.group(end) = this.clusTb.group(this.clusTb.clusId==oldId(k));
            end
            
            % Compute new cluster metrics
            disp('Computing new cluster metrics');
            this.clusTb.numSpikes(isOldClus) = 0; % clear spike number of the old clusters
            this.ComputeContamStats([newId reId]);
            this.ComputeClusterDepth([newId reId]);
        end
        
        % Plotting
        function PlotWaveform(this, spkInd, varargin)
            % 
            
            p = inputParser;
            p.KeepUnmatched = true;
            p.addOptional('N', min(50, numel(spkInd)), @(x) isnumeric(x) && isscalar(x));
            p.addParameter('Color', [0 0 0], @(x) ismember(numel(x), [3 4]));
            p.parse(varargin{:});
            N = p.Results.N;
            cc = p.Results.Color;
            
            if numel(cc) == 3
                cc(4) = 0.15;
            end
            
            spkInd = randsample(spkInd, N);
            [W, chWins, tmWins] = this.ExtractWaveform(spkInd, p.Unmatched);
            W = double(W);
            
            [nCh, nTm, nW] = size(W);
            
            t = ((0:nTm-1)'-nTm/2) / 30; % convert to millisecond
            for i = 1 : nW
                % Find Channel
                chInd = chWins(i,1) : chWins(i,2);
                y = this.chanTb.ycoords(chInd);
                
                % Plot a wavefrom
                MPlot.PlotTraceLadder(t, W(:,:,i)'*3/250, y', 'ColorArray', cc);
                hold on
            end
            
            ax = MPlot.Axes(gca);
            axis tight
        end
        
        function W = PlotClusterWaveform(this, cid, varargin)
            % Plot example waveform of given clusters
            % 
            %   PlotClusterWaveform(cid)
            %   PlotClusterWaveform(cid, N)
            %   PlotClusterWaveform(cid, N, spkInd)
            %   PlotClusterWaveform(..., 'Memmap', {})
            %   PlotClusterWaveform(..., 'Color', [0 0 0])
            %   PlotClusterWaveform(..., 'Style', 'each')
            %   PlotClusterWaveform(..., 'ShowCentroids', false)
            %   PlotClusterWaveform(..., 'Axes', gca)
            % 
            % Inputs
            %   cid                 A vector of cluster IDs (zero-based).
            %   N                   The number of waveform to randomly sample and plot for each cluster. 
            %                       The default is 50.
            %   spkInd              A vector (for a single cluster) of relative spike indices (logical 
            %                       or numeric) within a cluster, or a cell array of such vectors (for 
            %                       multiple clusters). This mask is applied before randomly sampling the 
            %                       N waveform. Empty vector(s) means no masking (default).
            %   'Memmap'            A cell array of memmapfile objects for cached waveform data. If any 
            %                       element is empty, waveform of that cluster will be extracted from 
            %                       this.mdat.
            %   'Color'             
            %   'Style'             
            %   'ShowCentroids'     
            %   'Axes'              
            % Output
            %   W                   Extracted or loaded waveform.
            
            p = inputParser;
            p.KeepUnmatched = true;
            p.addOptional('N', 50, @(x) isnumeric(x) && isscalar(x));
            p.addOptional('spkInd', [], @(x) ~ischar(x));
            p.addParameter('Color', [0 0 0], @(x) size(x,2)==3);
            p.addParameter('Style', 'each', @(x) ismember(lower(x), {'each', 'average'}));
            p.addParameter('ShowCentroids', false, @islogical);
            p.addParameter('Axes', gca, @(x) isa(x, 'matlab.graphics.axis.Axes'));
            p.addParameter('Memmap', {});
            p.parse(varargin{:});
            N = p.Results.N;
            I = p.Results.spkInd;
            isShowCentroids = p.Results.ShowCentroids;
            cc = p.Results.Color;
            ax = p.Results.Axes;
            
            nClus = numel(cid);
            if ~iscell(I)
                I = repmat({I}, [nClus 1]);
            end
            if size(cc,1) < numel(cid)
                cc = repmat(cc, [nClus 1]);
            end
            
            W = cell(size(cid));
            
            for k = 1 : nClus
                % Find spike indices
                spkMask = this.spkTb.clusId==cid(k);
                spkInd = find(spkMask);
                assert(~isempty(spkInd), 'No spike is found for cluster #%i', cid(k));
                if ~isempty(I{k})
                    spkInd = spkInd(I{k});
                end
                if numel(spkInd) > N
                    spkInd = randsample(spkInd, N);
                end
                
                % Extract waveform
                [w, chWins, tmWins] = this.ExtractWaveform(spkInd, p.Unmatched);
                W{k} = w;
                w = double(w);
                
                % Append NaN to the end of each trace as separator
                w = permute(w, [2 1 3]);
                w(end+1,:,:) = NaN;
                [nTm, nCh, nW] = size(w);
                
                % Scale waveform and add depth
                for i = 1 : nW
                    chInd = chWins(i,1) : chWins(i,2);
                    y = this.chanTb.ycoords(chInd);
                    w(:,:,i) = w(:,:,i)*3/250 + y(:)';
                end
                
                % Make time coordinates
                t = ((0:nTm-1)'-nTm/2) / 30; % convert to millisecond
                tShift = (k-1)*3;
                t = t + tShift;
                T = repmat(t, [1 nCh nW]);
                
                % Compute transparency
                a = min(1, 7.5/nW); % 7.5 is just a value that works fine for the visual
                
                % Plot
                plot(ax, T(:), w(:), 'Color', [cc(k,:) a]);
                hold(ax, 'on');
                
                % Label spike centroids
                if isShowCentroids && ismember('centCoords', this.spkTb.Properties.VariableNames)
                    y = this.spkTb.centCoords(spkInd, 2);
                    t = zeros(size(y)) + tShift;
                    plot(t, y, '.', 'Color', cc(k,:), 'MarkerSize', 8);
                end
            end
            MPlot.Axes(ax);
            axis(ax,'tight');
        end
        
        function PlotClusterTemplate(this, cid, varargin)
            % Plot templates of a given cluster
            % 
            %   PlotClusterTemplate(cid)
            %   PlotClusterTemplate(cid, spkInd)
            %   PlotClusterTemplate(..., 'Color', [0 0 0])
            %   PlotClusterTemplate(..., 'Axes', gca)
            % 
            % Inputs
            %   cid                 A vector of cluster IDs (zero-based).
            %   spkInd              A vector (for a single cluster) of relative spike indices (logical 
            %                       or numeric) within a cluster, or a cell array of such vectors (for 
            %                       multiple clusters). Empty vector(s) means no masking (default).
            %   'Color'             
            %   'Axes'              
            
            p = inputParser;
            p.KeepUnmatched = true;
            p.addOptional('spkInd', [], @(x) ~ischar(x));
            p.addParameter('Color', [0 0 0], @(x) size(x,2)==3);
            p.addParameter('Axes', gca, @(x) isa(x, 'matlab.graphics.axis.Axes'));
            p.parse(varargin{:});
            I = p.Results.spkInd;
            cc = p.Results.Color;
            ax = p.Results.Axes;
            
            nClus = numel(cid);
            if ~iscell(I)
                I = repmat({I}, [nClus 1]);
            end
            if size(cc,1) < numel(cid)
                cc = repmat(cc, [nClus 1]);
            end
            
            for k = 1 : nClus
                % Find spike indices
                spkMask = this.spkTb.clusId==cid(k);
                spkInd = find(spkMask);
                assert(~isempty(spkInd), 'No spike is found for cluster #%i', cid(k));
                if ~isempty(I{k})
                    spkInd = spkInd(I{k});
                end
                
                % Find unique template IDs
                tpId = this.spkTb.tempId(spkInd);
                [G, tpId] = findgroups(tpId);
                
                % Compute mean template amplitude
                tpAmp = this.spkTb.tempAmp(spkInd);
                tpAmp = splitapply(@median, tpAmp, G);
                
                % Compute transparency by the proportion of each template
                N = splitapply(@numel, spkInd, G);
                a = max(N/sum(N), 0.05); % set a minimal transparency
                
                for i = 1 : numel(tpId)
                    % Get template
%                     m = this.tempTb.tempId == tpId(i);
                    m = tpId(i) + 1;
                    chInd = this.tempTb.centIdx(m) + (-5:4); % take the middle 10 channels
                    y = this.chanTb.ycoords(chInd);
                    tp = this.tempTb.temp{m}(chInd,:)*tpAmp(i)*3 + y; % add depth
                    
                    % Vectorize template
                    tp = tp';
                    tp(end+1,:,:) = NaN;
                    [nTm, nCh] = size(tp);
                    
                    % Make time
                    t = ((0:nTm-1)'-nTm/2)/30; % convert to millisecond
                    tShift = (k-1)*3;
                    t = t + tShift;
                    t = repmat(t, [1 nCh]);
                    
                    % Plot
                    plot(ax, t(:), tp(:), 'Color', [cc(k,:) a(i)]);
                    hold(ax, 'on');
                end
            end
            MPlot.Axes(ax);
            axis(ax,'tight');
        end
        
        function PlotCCG(this, cid, varargin)
            % 
            
            p = inputParser();
            p.addParameter('Color', lines(1), @(x) size(x,2)==3);
            p.addParameter('ShowMetrics', true, @islogical);
            p.addParameter('Figure', gcf, @(x) isa(x, 'matlab.ui.Figure'));
            p.parse(varargin{:});
            cc = p.Results.Color;
            isShowMetrics = p.Results.ShowMetrics;
            fig = p.Results.Figure;
            
            % Not plotting for more than 10 clusters
            N = numel(cid);
            if N > 10
                return
            end
            
            if size(cc,1) == 1
                cc = repmat(cc, [N 1]);
            end
            
            % Find spikes and metrics for each cluster
            tSpk = cell(size(cid));
            for i = 1 : N
                m = this.spkTb.clusId==cid(i);
                tSpk{i} = this.spkTb.timeInd(m) / 30e3;
                tSpk{i} = unique(tSpk{i}); % remove duplicates
            end
            
            % Compute CCG
            tEdge = -0.02 : 0.5e-3 : 0.02;
            ccg = MNeuro.CCG(tEdge, tSpk{:});
            
            % Plot
            tBin = tEdge(2:end) - diff(tEdge)/2;
            for i = 1 : N
                for j = i : N
                    ax = subplot(N, N, (i-1)*N+j, 'Parent', fig);
                    h = bar(ax, tBin*1e3, squeeze(ccg(i,j,:)), 'hist');
                    h.EdgeColor = 'none';
                    h.FaceColor = (cc(i,:)+cc(j,:)) / 2;
                    MPlot.Axes(ax);
                    ax.XLim = tEdge([1 end])*1e3;
                    
                    if i == j
                        c = cid(i);
                        if isShowMetrics && all(ismember({'RPV', 'contam'}, this.clusTb.Properties.VariableNames))
                            m = this.clusTb.clusId==c;
                            rv = this.clusTb.RPV(m);
                            con = this.clusTb.contam(m);
                            ax.Title.String = sprintf('Cluster %i: RPV %.1f / CR %.0f', c, rv, con);
                        else
                            ax.Title.String = ['Cluster ' num2str(c)];
                        end
                    end
                end
            end
        end
        
    end
    
    methods(Static)
        function [W, chWins, tmWins] = LoadWaveform(mmap, wInd, varargin)
            % Load cached waveform from memmapfile
            
            [nCh, nTm, nW] = size(mmap.Data.W);
            
            p = inputParser();
            p.addParameter('NumSamples', nTm, @(x) isscalar(x) && isnumeric(x));
            p.addParameter('NumChannels', nCh, @(x) isscalar(x) && isnumeric(x));
            p.parse(varargin{:});
            nTmOut = min(p.Results.NumSamples, nTm);
            nChOut = min(p.Results.NumChannels, nCh);
            
            % Find indices
            chInd = (1:nChOut) - floor(nChOut/2) + floor(nCh/2);
            tmInd = (1:nTmOut) - floor(nTmOut/2) + floor(nTm/2);
            if nargin < 3 || isempty(wInd)
                wInd = 1 : nW;
            end
            
            % Read data
            W = mmap.Data.W(chInd, tmInd, wInd);
            t = mmap.Data.t(wInd); % spike time sample indices
            c = mmap.Data.c(wInd); % spike template center channel indices
            
            % Make time windows
            tmWins = t(:) + tmInd([1 end]);
            chWins = c(:) + chInd([1 end]);
        end
        
    end
end
