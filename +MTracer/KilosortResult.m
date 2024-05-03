classdef KilosortResult < MTracer.SortingResult
    %KILOSORTRESULT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        ksFolder = '';
        rezLite = struct;
        chanMapFile = '';
        chanMapName = '';
    end
    
    methods
        function this = KilosortResult()
            % do nothing
        end
        
        function kr = Duplicate(this)
            % Make a hard copy of the object
            if isempty(this)
                kr = MTracer.KilosortResult();
                return
            end
            kr = MTracer.KilosortResult();
            p = properties(this);
            for i = 1 : numel(p)
                kr.(p{i}) = this.(p{i});
            end
        end
        
        function ImportData(this, ksFolder, chanMapFile)
            % Read or link data files, and organize results into tables
            
            % Get the Kilosort folder if not provided
            if nargin < 2 || isempty(ksFolder)
                ksFolder = MBrowse.Folder([], 'Select the Kilosort output folder');
                if ~ksFolder
                    return
                end
            end
            
            
            % Load channel info
            if nargin < 3 || isempty(chanMapFile)
                chanMapFile = 'NP1_NHP_HalfCol_kilosortChanMap.mat';
                warning('ChanMapFile not specified. Use default channel map from %s', chanMapFile);
            end
            this.LoadChannelMap(chanMapFile);
            
            
            % Map binary file
            mdat = MKilosort2.MapDatFile(ksFolder);
            
            % Load sampling info
            rezFile = fullfile(ksFolder, 'rez.mat');
            opsFile = fullfile(ksFolder, 'ops.mat');
            paramsFile = fullfile(ksFolder, 'params.py');
            if exist(rezFile, 'file')
                this.LoadRez(rezFile);
                
            elseif exist(opsFile, 'file')
                load(opsFile, 'ops');
                ops.tstart = round(ops.trange(1) * ops.fs);
                if isinf(ops.trange(2))
                    ops.tend = ops.tstart + size(mdat.Data.V, 2);
                else
                    ops.tend = round(ops.trange(2) * ops.fs);
                end
                this.sampleOffset = ops.tstart;
                this.samplingRate = ops.fs;
                this.rezLite.ops = ops;
                
%             elseif exist(paramsFile, 'file')
%                 ops = MKilosort2.ReadParamsPy(paramsFile);
%                 ops.fs = ops.sample_rate;
%                 ops.tstart = ops.offset;
%                 ops.tend = nSample;
%                 ops.trange = [ops.tstart ops.tend];
%                 this.sampleOffset = ops.tstart;
%                 this.samplingRate = ops.fs;
%                 this.rezLite.ops = ops;
                
            else
                error("Cannot find files with necessary metadata.");
            end
            
            
            % Read spike times
            origSpkTimeIndFile = fullfile(ksFolder, 'spike_times_original.npy');
            if exist(origSpkTimeIndFile, 'file')
                tInd = readNPY(origSpkTimeIndFile); % sample indices of all spikes
            else
                tInd = readNPY(fullfile(ksFolder, 'spike_times.npy')); % sample indices of all spikes
            end
            
            % Make spike table
            sTb = table;
            sTb.spkId = (1:numel(tInd))' - 1; % assign spike IDs (zero-based)
            sTb.timeInd = tInd;
            sTb.timeSec = double(tInd) / 30e3;
            
            spkTimeSecFile = fullfile(ksFolder, 'spike_times_in_sec_adj.npy');
            if exist(spkTimeSecFile, 'file')
                sTb.timeSecAdj = readNPY(spkTimeSecFile); % adjusted time of all spikes
            else
                sTb.timeSecAdj = NaN(size(tInd));
            end
            
            sTb.clusId = readNPY(fullfile(ksFolder, 'spike_clusters.npy'));     % cluster ID of all spikes, zero-based
            sTb.tempId = readNPY(fullfile(ksFolder, 'spike_templates.npy'));    % template ID of all spikes, zero-based
            sTb.tempAmp = readNPY(fullfile(ksFolder, 'amplitudes.npy'));        % template amplitude of all spikes
            
            for i = 1 : width(sTb)
                sTb.(i) = double(sTb.(i));
            end
            
            
            % Make template table
            temp = readNPY(fullfile(ksFolder, 'templates.npy'));    % [nTemp, time, chan]
            temp = permute(temp, [3 2 1]);                          % [chan, time, nTemp]
            temp = temp(this.chanTb.sortInd,:,:);                   % apply sorting by depth
            
            pcChanInd = readNPY(fullfile(ksFolder, 'pc_feature_ind.npy'));
            [pcChanInd, pcSortInd] = this.IConvertPcChanInd(pcChanInd);
            
            pcWeights = readNPY(fullfile(ksFolder, 'pc_features.npy'));
            pcWeights = permute(pcWeights, [3 2 1]); % to [chan, pc, spk]
            
            tTb = table;
            for i = size(temp,3) : -1 : 1
                % Take a template
                tp = temp(:,:,i);
                
                % Find channel indices and the index of the template center
                ind = double(pcChanInd(i,:)');
                if isempty(ind)
                    warning('All the values in template #%i (zero-based) are zero. This template has %i spikes.', ...
                        i-1, sum(sTb.tempId==i-1));
                    cIdx = round(mean(ind));
                else
                    [~, k] = MNeuro.ComputeWaveformCenter(tp(ind,:), [], 'power', 'peak');
                    cIdx = ind(k);
                end
                
                % Add a row
                tTb.tempId(i) = i-1;
                tTb.temp{i} = tp;
                tTb.chanInd{i} = ind;
                tTb.chanX{i} = this.chanTb.xcoords(ind);
                tTb.chanY{i} = this.chanTb.ycoords(ind);
                tTb.centIdx(i) = cIdx;
                
                % Sort channels of pcWeights
                m = sTb.tempId == i-1;
                pcWeights(:,:,m) = pcWeights(pcSortInd(i,:), :, m);
            end
            
            pcWeights = mat2cell(pcWeights, size(pcWeights,1), size(pcWeights,2), ones(size(tInd)));
            sTb.pcWeights = pcWeights(:);
            
            
            % Make cluster table
            cTb = table;
            [~, cTb.clusId] = findgroups(sTb.clusId);
            cTb.group = repmat({''}, size(cTb.clusId));
            
            gTb = MKilosort2.ReadTsvFiles(ksFolder, 'cluster_group.tsv'); % KS or manually labeled quality for each cluster
            gTb.Properties.VariableNames{2} = 'group'; % always name the label column 'group'
            for i = 1 : height(gTb)
                m = cTb.clusId == gTb.cluster_id(i);
                if any(m)
                    cTb.group{m} = gTb.group{i};
                else
                    % fprintf('Skip cluster #%i (zero-based) for it is not associated with any spike.\n', gTb.cluster_id(i));
                end
            end
            
            
            % Save variables
            this.ksFolder = ksFolder;
            this.mdat = mdat;
            this.spkTb = sTb;
            this.tempTb = tTb;
            this.clusTb = cTb;
        end
        
        function LoadChannelMap(this, chanMapFile)
            % Load channel configuration to this.chanTb
            % 1) The rows of the table are sorted spatially in the order of increasing shank ID, 
            %    decreasing y coordinate (distance to tip), and increasing x coordinate
            % 2) Channels that are not present in voltage data should not be in the table.
            
            [~, this.chanTb] = MKilosort2.LoadChanMap2Table(chanMapFile);
        end
        
        function rez = LoadRez(this, rezFile)
            % Load a subset of data from rez.mat
            
            disp('Loading rez.mat');
            load(rezFile, 'rez');
            
            this.sampleOffset = rez.ops.tstart;
            this.pcBases = rez.ops.wPCA;
            
            fields2copy = {'ops', 'dshift', 'st0'};
            for i = 1 : numel(fields2copy)
                fn = fields2copy{i};
                if ~ismember(fn, fieldnames(rez))
                    warning("Field '%s' is missing from rez", fn);
                    continue
                end
                this.rezLite.(fn) = rez.(fn);
            end
        end
        
        function ExportData(this, ksFolder)
            % 
            
            if nargin < 2 || isempty(ksFolder)
                ksFolder = this.ksFolder;
            end
            
            % Always backup the files to be modified
            bkFolder = fullfile(ksFolder, ['backup_results_' datestr(now, 'yyyy-mm-dd_HH-MM-SS')]);
            mkdir(bkFolder);
            fileNames = { ...
                'spike_clusters.npy', 'spike_centroids.npy', 'spike_center_channels.npy', 'amplitudes.npy', ...
                'cluster_Amplitude.tsv', 'cluster_ContamPct.tsv', 'cluster_group.tsv', 'cluster_info.tsv', 'cluster_KSLabel.tsv'};
            for i = 1 : numel(fileNames)
                fromPath = fullfile(ksFolder, fileNames{i});
                if exist(fromPath, 'file')
                    movefile(fromPath, fullfile(bkFolder, fileNames{i}));
                end
            end
            
            % Save spike cluster IDs
            writeNPY(uint32(this.spkTb.clusId), fullfile(ksFolder, 'spike_clusters.npy'));
            
            % Save spike coordinates
            writeNPY(this.spkTb.centCoords, fullfile(ksFolder, 'spike_centroids.npy'));
%             writeNPY(this.spkTb.centChanInd, fullfile(ksFolder, 'spike_center_channels.npy'));
            
            % Save spike amplitudes
            writeNPY(this.spkTb.tempAmp, fullfile(ksFolder, 'amplitudes.npy'));
            
            % 
            cTb = this.clusTb;
            
%             fileID = fopen(fullfile(ksFolder, 'cluster_KSLabel.tsv'),'w');
%             fprintf(fileID, 'cluster_id%sKSLabel', char(9));
%             fprintf(fileID, char([13 10]));
%             for j = 1 : height(cTb)
%                 fprintf(fileID, '%d%s%s', cTb.clusId(j), char(9), cTb.group{j});
%                 fprintf(fileID, char([13 10]));
%             end
%             fclose(fileID);
            
            groupFile = fopen(fullfile(ksFolder, 'cluster_group.tsv'),'w');
            fprintf(groupFile, 'cluster_id%sgroup', char(9));
            fprintf(groupFile, char([13 10]));
            for j = 1 : height(cTb)
                fprintf(groupFile, '%d%s%s', cTb.clusId(j), char(9), cTb.group{j});
                fprintf(groupFile, char([13 10]));
            end
            fclose(groupFile);
            
            contamFile = fopen(fullfile(ksFolder, 'cluster_ContamPct.tsv'),'w');
            fprintf(contamFile, 'cluster_id%sContamPct', char(9));
            fprintf(contamFile, char([13 10]));
            for j = 1 : height(cTb)
                fprintf(contamFile, '%d%s%.1f', cTb.clusId(j), char(9), cTb.contam(j));
                fprintf(contamFile, char([13 10]));
            end
            fclose(contamFile);
            
%             fileIDA = fopen(fullfile(ksFolder, 'cluster_Amplitude.tsv'),'w');
%             fprintf(fileIDA, 'cluster_id%sAmplitude', char(9));
%             fprintf(fileIDA, char([13 10]));
%             for j = 1 : height(cTb)
%                 fprintf(fileIDA, '%d%s%.1f', cTb.clusId(j), char(9), cTb.amplitude(j));
%                 fprintf(fileIDA, char([13 10]));
%             end
%             fclose(fileIDA);
        end

        function cvsTb = SaveClusSummaryCSV(this, ksFolder)
            % Save a subset of cluster table as a CSV file
            
            cTb = this.clusTb;
            cvsTb = table;
            for i = 1 : width(cTb)
                col = cTb.(i);
                if isvector(col) && (isnumeric(col) || isstring(col) || iscellstr(col))
                    cvsTb.(cTb.Properties.VariableNames{i}) = cTb.(i);
                end
            end
            
            if nargin > 1
                writetable(cvsTb, fullfile(ksFolder, 'cluster_summary.csv'));
            end
        end
        
        function [pcChInd, pcSortInd] = IConvertPcChanInd(this, pcChInd)
            % Convert channel indices in pc_feature_ind.npy to row indices of this.chanTb
            %   [pcChInd, pcSortInd] = IConvertPcChanInd(pcChInd)
            % Input
            %   pcChInd         Original zero-based channel indices of connected channel(s).
            %                   pcChInd + 1 are the same as this.chanTb.sortInd.
            % Output
            %   pcChInd         Row indices of this.chanTb.
            %   I               Indices that sort pcChInd the same way as this.chanTb
            
            pcChInd = pcChInd + 1;
            pcSortInd = zeros(size(pcChInd));
            for i = 1 : size(pcChInd, 1)
                ind = pcChInd(i,:);
                [ind, pcSortInd(i,:)] = MMath.SortLike(ind, this.chanTb.sortInd, false);
                m = ismember(this.chanTb.sortInd, ind);
                pcChInd(i,:) = find(m);
            end
        end
        
    end
end
