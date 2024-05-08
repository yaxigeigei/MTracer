classdef MKilosort2
    %MKilosort2 Summary of this class goes here
    %   Detailed explanation goes here
    
    methods(Static)
        function varargout = Sort(varargin)
            % Run spike sorting on a binary file. Results will be saved in the folder of .dat file.
            % Auto merge is disabled and no intermediate data is saved. 
            % 
            %   MKilosort2.Sort()
            %   MKilosort2.Sort(binFilePath)
            %   MKilosort2.Sort(binFilePath, outFolder)
            %   MKilosort2.Sort(..., 'ChannelMapFile', [])
            %   MKilosort2.Sort(..., 'ConfigFunc', [])
            %   rez = MKilosort.Sort(...)
            % 
            % Inputs
            %   binFile             A binary file of extracellular recording data. If not specified, a file
            %                       selection window will be prompted. 
            %   outFolder           Kilosort output folder. 
            %   
            %   If either one of 'ChannelMapFunc' or 'ConfigFunc' is not provided, a dialog box will allow you 
            %   to choose a probe type for the preset configuration and channel map. 
            %   
            %   'ChannelMapFile'    A function handle that takes the folder path of the .dat file as input 
            %                       and saves a chanMap.mat file to that folder. Examples can be found in the 
            %                       methods of this class (e.g. @MKilosort.SaveChanMapH3, modified from 
            %                       createChannelMapFile.m in the configFiles folder of Kilosort repository). 
            %   'ConfigFunc'        A function handle that takes the path of the .dat file as input and 
            %                       returns a structure of options. Examples can be found in the methods of this 
            %                       class (e.g. @MKilosort.Config64, modified from StandardConfig_MOVEME.m in 
            %                       the configFiles folder of Kilosort repository).
            % Output
            %   rez                 The raw clustering results that Kilosort outputs. 
            
            % Parse inputs
            p = inputParser();
            p.KeepUnmatched = true;
            p.addOptional('binFile', [], @(x) ischar(x) || isstring(x) || isempty(x));
            p.addOptional('outFolder', [], @(x) ischar(x) || isstring(x) || isempty(x));
            p.addParameter('ChannelMapFile', [], @(x) ischar(x) || isstring(x) || isempty(x));
            p.addParameter('ConfigFunc', [], @(x) isa(x, 'function_handle'));
            p.addParameter('DriftMapOnly', false, @islogical);
            p.parse(varargin{:});
            binFile = char(p.Results.binFile);
            outDir = char(p.Results.outFolder);
            chanMapFile = p.Results.ChannelMapFile;
            fConfig = p.Results.ConfigFunc;
            isMapOnly = p.Results.DriftMapOnly;
            
            % Find the binary file
            if isempty(binFile)
                binFile = MBrowse.File([], 'Select a binary file', {'*.dat', '*.bin'});
            end
            if isempty(binFile)
                return
            end
            
            % Set the output folder
            if isempty(outDir)
                outDir = MBrowse.File([], 'Select the kilosort output folder');
            end
            if isempty(outDir)
                return
            end
            if ~exist(outDir, 'dir')
                mkdir(outDir);
            end
            
            if isempty(chanMapFile)
                % Extract channel map from ap.meta file
                metaFile = strrep(binFile, '.bin', '.meta');
                assert(logical(exist(metaFile, 'file')), "Channel map is not provided and the .meta file is not present with .bin file for extracting channel map info.");
                chanMapFile = fullfile(outDir, "chanMap.mat");
                MSpikeGLX.MetaToChanMap(metaFile, chanMapFile);
            else
                % Make a copy of the channel map file in the output directory
                copyfile(chanMapFile, fullfile(outDir, "chanMap.mat"));
            end
            
            if isempty(fConfig)
                % Use default sorting configuration
                fConfig = @MKilosort2.Config384;
            end
            
            
            % Get default options
            ops.fbinary = binFile;
            ops.fproc   = fullfile(outDir, 'temp_wh.dat'); % proc file on a fast SSD
            ops.chanMap = char(chanMapFile);
            ops.trange  = [0 Inf]; % time range to sort
            ops.NchanTOT = 385; % total number of channels in your recording
            ops = fConfig(ops);
            
            % Overwrite options from function args
            opsNames = intersect(fieldnames(ops), fieldnames(p.Unmatched));
            for k = 1 : numel(opsNames)
                n = opsNames{k};
                ops.(n) = p.Unmatched.(n);
            end
            
            % Initialize GPU (will erase any existing GPU arrays)
            gpuDevice(1);
            
            % Preprocess data and extract spikes for initialization
            ResetGPU();
            rez = preprocessDataSub(ops);
            
            % Save data for raw spike map
            if isMapOnly
                rez = datashift2(rez, 0); % last input is for shifting data
                SaveRez(rez, outDir);
                if nargout > 0
                    varargout{1} = rez;
                end
                return
            else
                % NEW STEP TO DO DATA REGISTRATION
                rez = datashift2(rez, 1); % last input is for shifting data
            end
            
            % ORDER OF BATCHES IS NOW RANDOM, controlled by random number generator
            iseed = 1;
            
            % main tracking and template matching algorithm
            rez = learnAndSolve8b(rez, iseed);
            
            % OPTIONAL: remove double-counted spikes - solves issue in which individual spikes are assigned to multiple templates.
            % See issue 29: https://github.com/MouseLand/Kilosort/issues/29
            %rez = remove_ks2_duplicate_spikes(rez);
            
            % final merges
            rez = find_merges(rez, 1);
            
            % final splits by SVD
            rez = splitAllClusters(rez, 1);
            
            % decide on cutoff
            rez = set_cutoff(rez);
            
            % eliminate widely spread waveforms (likely noise)
            rez.good = get_good_units(rez);
            
            fprintf('found %d good units \n', sum(rez.good > 0))
            
            % write to Phy
            fprintf('Saving results to Phy  \n')
            rezToPhy(rez, outDir);
            
            % discard features in final rez file (too slow to save)
            rez.cProj = [];
            rez.cProjPC = [];
            
            % final time sorting of spikes, for apps that use st3 directly
            [~, isort]   = sortrows(rez.st3);
            rez.st3      = rez.st3(isort, :);
            
            % save rez.mat
            SaveRez(rez, outDir);
            
            % save(fullfile(fpath, 'rez.mat'), 'rez', '-v7.3');
            if nargout > 0
                varargout{1} = rez;
            end
            
            % Helper functions
            function ResetGPU()
                reset(parallel.gpu.GPUDevice.current);
            end
            function SaveRez(rez, outDir)
                % Ensure all GPU arrays are transferred to CPU side before saving to .mat
                rez_fields = fieldnames(rez);
                for i = 1:numel(rez_fields)
                    field_name = rez_fields{i};
                    if(isa(rez.(field_name), 'gpuArray'))
                        rez.(field_name) = gather(rez.(field_name));
                    end
                end
                
                % save final results as rez2
                fprintf('Saving final results in rez  \n')
                fname = fullfile(outDir, 'rez.mat');
                save(fname, 'rez', '-v7.3');
            end
        end
        
        function ops = Config384(ops)
            % 
            
            % sample rate
            ops.fs = 30000;
            
            % frequency for high pass filtering (150)
            ops.fshigh = 300;
            
            % threshold on projections (like in Kilosort1, can be different for last pass like [10 4])
            ops.Th = [10 4];
            
            % how important is the amplitude penalty (like in Kilosort1, 0 means not used, 10 is average, 50 is a lot)
            ops.lam = 10;
            
            % splitting a cluster at the end requires at least this much isolation for each sub-cluster (max = 1)
            ops.AUCsplit = 0.9;
            
            % minimum spike rate (Hz), if a cluster falls below this for too long it gets removed
            ops.minFR = 1/50;
            
            % number of samples to average over (annealed from first to second value)
            ops.momentum = [20 400];
            
            % spatial constant in um for computing residual variance of spike
            ops.sigmaMask = 30;
            
            % threshold crossings for pre-clustering (in PCA projection space)
            ops.ThPre = 8;
            
            % spatial scale for datashift kernel
            ops.sig = 20;
            
            % type of data shifting (0 = none, 1 = rigid, 2 = nonrigid)
            ops.nblocks = 5;
            
            
            % danger, changing these settings can lead to fatal errors
            % options for determining PCs
            ops.spkTh           = -6;      % spike threshold in standard deviations (-6)
            ops.reorder         = 1;       % whether to reorder batches for drift correction.
            ops.nskip           = 25;  % how many batches to skip for determining spike PCs
            
            ops.GPU                 = 1; % has to be 1, no CPU version yet, sorry
            % ops.Nfilt               = 1024; % max number of clusters
            ops.nfilt_factor        = 4; % max number of clusters per good channel (even temporary ones)
            ops.ntbuff              = 64;    % samples of symmetrical buffer for whitening and spike detection
            ops.NT                  = 64*1024+ ops.ntbuff; % must be multiple of 32 + ntbuff. This is the batch size (try decreasing if out of memory).
            ops.whiteningRange      = 32; % number of channels to use for whitening each channel
            ops.nSkipCov            = 25; % compute whitening matrix from every N-th batch
            ops.scaleproc           = 200;   % int16 scaling of whitened data
            ops.nPCs                = 3; % how many PCs to project the spikes into
            ops.useRAM              = 0; % not yet available
            
        end
        
        function [chanTb, chanTbKS] = LoadChanMap2Table(chanMapFile)
            % Load channel info from a channel map .mat file to a table
            % 
            %   [chanTb, chanTbKS] = MKilosort2.LoadChanMap2Table(chanMapFile)
            % 
            % Input:
            %   chanMapFile         Path of a channel map MAT file, e.g. NP1_NHP_HalfCol_kilosortChanMap.mat.
            %                       The following variables are expected from the file:
            %                       chanMap0ind     Zero-based channel indices.
            %                       shankInd        One-based shank indices.
            %                       xcoords         X coordinates of channels in microns.
            %                       ycoords         Y coordinates of channels in microns.
            %                       connected       Binary vector indicating what channels are valid.
            %                                       Reference and bad channels are usually set to 0.
            % Outputs:
            %   chanTb              A table that organizes varibales from chanMapFile in columns.
            %   chanTbKS            Similar to chanTb but only with connected channels. This table is consistent
            %                       with kilosort output where only connected channels are indexed. In addition, 
            %                       the rows are sorted spatially in the order of increasing shank ID, decreasing 
            %                       y coordinate (distance to tip), and increasing x coordinate.
            % 
            
            % Load mat file
            s = load(chanMapFile);
            
            % Make a table of channel info
            chanTb = table;
            chanTb.chanId = s.chanMap0ind(:);
            if isfield(s, 'shankInd')
                chanTb.shankInd = s.shankInd(:);
            else
                chanTb.shankInd = s.kcoords(:);
            end
            chanTb.xcoords = s.xcoords(:);
            chanTb.ycoords = s.ycoords(:);
            chanTb.isConnected = logical(s.connected(:));
            
            % Remove the unconnected channel(s)
            chanTbKS = chanTb;
            chanTbKS = chanTbKS(chanTbKS.isConnected, :);
            
            % Sort table
            [chanTbKS, I] = sortrows(chanTbKS, {'shankInd', 'ycoords', 'xcoords'}, {'ascend', 'descend', 'ascend'});
            chanTbKS.sortInd = I;
        end
        
        function s = ReadParamsPy(filePath)
            % Get parameters from params.py by running lines in MATLAB
            %   these include dat_path, dtype, n_channels_dat, sample_rate ...
            
            fid = fopen(filePath);
            s = struct;
            while ~feof(fid)
                try
                    l = fgetl(fid);
                    eval(['s.' l ';']);
                catch
                    %warning('''%s'' cannot be evaluated.', l);
                end
            end
            fclose(fid);
        end
        
        function varargout = ReadTsvFiles(ksDir, varargin)
            % Read various tsv files generated by kilosort and Phy
            %
            %   tb = MKilosort2.ReadTsvFiles(ksDir, fileName)
            %   [tb1, tb2, ...] = MKilosort2.ReadTsvFiles(ksDir, fileName1, fileName2, ...)
            %
            % Inputs:
            %   ksDir           Kilosort output direcoty where the tsv file(s) live.
            %   fileName        Name of the tsv file to read.
            %                   One can specify multiple files to read in one function call.
            %                   Supported tsv files include:
            %                   cluster_info.tsv, cluster_group.tsv, cluster_KSLabel.tsv
            % Outputs:
            %   tb              Table read from the tsv file. Mutiple tables are returned in 
            %                   separate variables.
            %
            
            varargout = cell(size(varargin));
            
            for k = 1 : numel(varargin)
                % Check file
                fileName = varargin{k};
                filePath = fullfile(ksDir, fileName);
                if ~exist(filePath, 'file')
                    warning('Cannot find %s', fileName);
                    continue
                end
                
                % Specify the pattern to read
                if strcmp(fileName, 'cluster_info.tsv')
                    formatStr = '%f %f %f %s %f %f %f %f %s %f %f';
                    headers = {'cluster_id', 'Amplitude', 'ContamPct', 'KSLabel', 'amp', 'ch', ...
                        'depth', 'fr', 'group', 'n_spikes', 'sh'};
                    
                elseif strcmp(fileName, 'cluster_group.tsv')
                    formatStr = '%f %s';
                    headers = {'cluster_id', 'group'};
                    
                elseif strcmp(fileName, 'cluster_KSLabel.tsv')
                    formatStr = '%f %s';
                    headers = {'cluster_id', 'KSLabel'};
                    
                else
                    warning('Reading %s is not supported', fileName);
                    continue
                end
                
                % Read file
                fid = fopen(filePath);
                C = textscan(fid, formatStr, 'HeaderLines', 1);
                fclose(fid);
                
                % Make a tbale
                tb = table;
                for i = 1 : numel(headers)
                    tb.(headers{i}) = C{i};
                end
                varargout{k} = tb;
            end
        end
        
        function mdat = MapDatFile(ksDir)
            % Make a memory map to temp_wh.dat
            
            % Get parameters from params.py by running lines in MATLAB
            s = MKilosort2.ReadParamsPy(fullfile(ksDir, 'params.py'));
            dat_path = s.dat_path;
            n_channels_dat = s.n_channels_dat;
            dtype = s.dtype;
            
            % Map binary data to memory
            if ~exist(dat_path, 'file')
                [~, datName, datExt] = fileparts(dat_path);
                dat_path = fullfile(ksDir, [datName datExt]);
            end
            if exist(dat_path, 'file')
                mdat = memmapfile(dat_path, 'Format', dtype);
                nSample = numel(mdat.Data) / n_channels_dat;
                mdat = memmapfile(dat_path, 'Format', {dtype, [n_channels_dat nSample], 'V'});
            else
                mdat = [];
            end
        end
        
        function sn = ReadSnippets(mmap, tmInd, tmWin, varargin)
            % Read snippets from mapped binary data array
            %
            %   sn = ReadSnippets(mmap, tmInd, tmWin)
            %   sn = ReadSnippets(mmap, tmInd, tmWin, chInd)
            %   sn = ReadSnippets(mmap, tmInd, tmWin, chInd, chWin)
            %   sn = ReadSnippets(mmap, tmInd, tmWin, ..., 'ChannelOrder', [])
            %   sn = ReadSnippets(mmap, tmInd, tmWin, ..., 'VoltScale', 1)
            %   sn = ReadSnippets(mmap, tmInd, tmWin, ..., 'Filter', [])
            % 
            % Inputs
            %   mmap                A memmapfile object to the binary data with the array size specified.
            %   tmInd               A vector of time indices.
            %   tmWin               The time window. The i-th snippet is read from tmInd(i)+tmWin(1) to tmInd(i)+tmWin(2), inclusively.
            %   chInd               A vector of channel indices. If scalar, the same channel will be extracted for all snippets.
            %                       The default is empty, reading all channels.
            %   chWin               The window of channels. The i-th snippet is read from tmInd(i)+tmWin(1) to tmInd(i)+tmWin(2), inclusively.
            %                       The default is [0 0], reading only the channel specified by chInd. chWin is ignored if chInd is empty.
            %   'ChannelOrder'      A vector of indices that reorder the channels in binary array before extracting snippets.
            %   'VoltScale'         A scaling factor applied to the data.
            %   'Filter'            A digitalFilter object used to filter the snippet data (by each channel in each snippet).
            % Output
            %   sn                  A channel-by-time-by-#snippets array of signal values.
            % 
            
            tmInd = double(tmInd(:));
            nSn = numel(tmInd);
            [nCh, nTm] = size(mmap.Data.V);
            
            p = inputParser();
            p.addOptional('chInd', [], @(x) isnumeric(x) && isvector(x));
            p.addOptional('chWin', [0 0], @(x) isnumeric(x) && numel(x)==2);
            p.addParameter('ChannelOrder', [], @(x) isnumeric(x) && isvector(x));
            p.addParameter('VoltScale', 1, @(x) isnumeric(x) && isvector(x));
            p.addParameter('Filter', [], @(x) isa(x, 'digitalFilter'));
            p.parse(varargin{:});
            chInd = double(p.Results.chInd(:));
            chWin = p.Results.chWin(:)';
            cOrder = p.Results.ChannelOrder;
            vScale = p.Results.VoltScale;
            filtObj = p.Results.Filter;
            
            if isempty(chInd)
                chInd = 1;
                chWin = [0 nCh-1];
            end
            if isscalar(chInd)
                chInd = ones(size(tmInd))*chInd;
            end
            
            if isempty(cOrder)
                cOrder = 1 : nCh;
            end
            
            % Read data
            tmWins = tmInd + tmWin;
            chWins = chInd + chWin;
            for i = numel(tmInd) : -1 : 1
                % Get channel indices
                ind1 = chWins(i,1) : chWins(i,2);
                isChOut = ind1 < 1 | ind1 > nCh;
                ind1 = min(max(ind1, 1), nCh);
                ind1 = cOrder(ind1);
                
                % Get sample indices
                ind2 = tmWins(i,1) : tmWins(i,2);
                ind2 = min(max(ind2, 1), nTm); % propagate end values to fill out of range samples
                
                % Read from mapped data
                sn(:,:,i) = mmap.Data.V(ind1, ind2);
                sn(isChOut,:,i) = 0; % reset values in out of range channels to zeros
            end
            
            if vScale ~= 1
                sn = double(sn) .* vScale(:); % apply voltage scaling
            end
            
            % Highpass filtering
            if ~isempty(filtObj)
                sn = double(sn);
                for i = 1 : nCh
                    for j = 1 : nSn
                        sn(i,:,j) = filtfilt(filtObj, sn(i,:,j));
                    end
                end
            end
        end
        
        function D = MakeHighpassFilter(fs)
            % Prepare parameters for filtering
            %   1) Phy uses bandpass between 500Hz and 14.25kHz(.475*sample_rate) to visulize waveform;
            %      here uses 300Hz highpass to be consistent with Kilosort configuration
            %   2) 82 samples per waveform is consistent with Kilosort template size
            %   3) Consistent with Phy, using samples 3 times of the filter order as margins when filtering
            if nargin < 1
                fs = 30e3;
            end
            D = designfilt('highpassiir', ...
                'PassbandFrequency', 300, ...
                'StopbandFrequency', 250, ...
                'StopbandAttenuation', 60, ...
                'PassbandRipple', 0.1, ...
                'SampleRate', fs, ...
                'DesignMethod', 'ellip'); % order of this filter is 8
        end
        
        function D = MakeLowpassFilter(fs)
            % Prepare a digitalFilter object for lowpass filtering of LFP
            if nargin < 1
                fs = 30e3;
            end
            D = designfilt('lowpassiir', ...
                'PassbandFrequency', 200, ...
                'StopbandFrequency', 250, ...
                'StopbandAttenuation', 60, ...
                'PassbandRipple', 0.1, ...
                'SampleRate', fs, ...
                'DesignMethod', 'ellip');
        end
        
    end
    
end
