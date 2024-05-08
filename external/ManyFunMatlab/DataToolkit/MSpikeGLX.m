classdef MSpikeGLX
    % Methods in this class are mostly copied from DemoReadSGLXData.m and SGLXMetaToCoords.m 
    % with the following changes
    %   1) In ReadBin, convert int16 to single instead of double
    %   2) Add high-level reader methods: ReadNI, ReadLFP, ReadMetaEZ
    %   3) Add helper methods: MetaToChanMap
    
    methods(Static)
        function meta = ReadMetaEZ(metaPath)
            % Read metadata with an optional single path input
            if nargin < 1 || isempty(metaPath)
                metaPath = MBrowse.File();
            end
            metaPath = char(metaPath);
            [folderPath, baseName] = fileparts(metaPath);
            meta = MSpikeGLX.ReadMeta([baseName '.meta'], folderPath);
        end
        
        function [meta, aiArray, diArray, t] = ReadNI(niPath)
            % 
            
            % Read metadata
            niPath = char(niPath);
            [folderPath, baseName] = fileparts(niPath);
            meta = MSpikeGLX.ReadMeta([baseName '.meta'], folderPath);
            
            % Read timeseries
            samp0 = 0;
            nSamp = Inf;
            dataArray = MSpikeGLX.ReadBin(samp0, nSamp, meta, [baseName '.bin'], folderPath)';
            nSamp = size(dataArray,1);
            t = (0 : nSamp-1)' / MSpikeGLX.SampRate(meta);
            
            % Extract analog channels
            if ~isempty(meta.niXAChans1)
                I = eval(meta.niXAChans1) + 1;
                dataArray = MSpikeGLX.GainCorrectNI(dataArray, I, meta); % gain correction
                aiArray = dataArray(:,I);
            else
                aiArray = zeros(nSamp,0);
            end
            
            % Extract digital channels
            if ~isempty(meta.niXDChans1)
                dwReq = 1;
                dLineList = 0:15;
                diArray = MSpikeGLX.ExtractDigital(dataArray, meta, dwReq, dLineList);
                diArray = logical(diArray);
            else
                diArray = false(nSamp,0);
            end
        end
        
        function [meta, lfpArray, t] = ReadLFP(lfPath)
            % 
            
            % Read metadata
            lfPath = char(lfPath);
            [folderPath, baseName] = fileparts(lfPath);
            meta = MSpikeGLX.ReadMeta([baseName '.meta'], folderPath);
            
            % Read timeseries
            samp0 = 0;
            nSamp = Inf;
            dataArray = MSpikeGLX.ReadBin(samp0, nSamp, meta, [baseName '.bin'], folderPath)';
            nSamp = size(dataArray,1);
            t = (0 : nSamp-1)' / MSpikeGLX.SampRate(meta);
            
            % Gain correction
            lfpInd = 1:384; % 385 is sync
            dataArray = MSpikeGLX.GainCorrectIM(dataArray, lfpInd, meta);
            lfpArray = dataArray(:,lfpInd);
        end
        
        function s = MetaToChanMap(metaPath, outPath)
            % Construct and optionally save Kilosort2 channel map from a SGLX .meta file
            % Adapted from Jennifer Colonell's SGLXMetaToCoords function
            % 
            %   s = MSpikeGLX.MetaToChanMap()
            %   s = MSpikeGLX.MetaToChanMap(metaPath)
            %   s = MSpikeGLX.MetaToChanMap(metaPath, outPath)
            % 
            % Inputs
            %   metaPath        Path of a .meta file.
            %   outPath         Path of the output file.
            % Output
            %   s               A struct with fields for KS2 channel map variables.
            % 
            
            if ~exist('metaPath', 'var') || isempty(metaPath)
                % Ask user for metadata file
                [metaName, metaFolder] = uigetfile('*.meta', 'Select Metadata File');
                if ~metaName
                    s = [];
                    return
                end
            else
                [metaFolder, metaBareName, metaExt] = fileparts(char(metaPath));
                metaName = [metaBareName, metaExt];
            end
            
            % Parse in file to get the metadata structure
            meta = MSpikeGLX.ReadMeta(metaName, metaFolder);
            
            % Get coordinates for saved channels from snsGeomMap, if present,
            % otherwise from snsShankMap
            if isfield(meta,'snsGeomMap')
                [nShank, shankWidth, shankPitch, shankInd, xCoord, yCoord, connected] = MSpikeGLX.geomMapToGeom(meta);
            elseif isfield(meta,'snsShankMap')
                [nShank, shankWidth, shankPitch, shankInd, xCoord, yCoord, connected] = MSpikeGLX.shankMapToGeom(meta);
            end
            
            % Build output name and write out file
            [~, metaBareName, ~] = fileparts(metaName);
            
            % Variables for KS2
            s = struct;
            s.name = metaBareName;
            s.chanMap = (1:numel(xCoord))';
            s.chanMap0ind = s.chanMap - 1;
            s.connected = logical(connected);
            s.xcoords = shankInd*shankPitch + xCoord; % KS2 not yet using kcoords, so x coord includes shank sep
            s.ycoords = yCoord;
            s.kcoords = shankInd + 1;
            s.shankInd = s.kcoords;
            
            % Save KS channel map
            if ~exist('outPath', 'var')
                return
            end
            if isempty(outPath)
                outPath = [metaBareName, '_kilosortChanMap.mat'];
            end
            save(outPath, '-struct', 's');
        end
        
        % From SpikeGLX_Datafile_Tools (need to adopt the new update)
        function DemoReadSGLXData()
            % Simple helper functions and MATLAB structures demonstrating
            % how to read and manipulate SpikeGLX meta and binary files.
            %
            % The most important part of the demo is ReadMeta().
            % Please read the comments for that function. Use of
            % the 'meta' structure will make your data handling
            % much easier!
            %
            
            % Ask user for binary file
            [binName,path] = uigetfile('*.bin', 'Select Binary File');
            
            % Parse the corresponding metafile
            meta = MSpikeGLX.ReadMeta(binName, path);
            
            % Get first one second of data
            nSamp = floor(1.0 * MSpikeGLX.SampRate(meta));
            dataArray = ReadBin(0, nSamp, meta, binName, path);
            
            dataType = 'A';         %set to 'A' for analog, 'D' for digital data
            
            % For an analog channel: gain correct saved channel ch (1-based for MATLAB).
            ch = 1;
            
            % For a digital channel: read this digital word dw in the saved file
            % (1-based). For imec data there is never more than one saved digital word.
            dw = 1;
            
            % Read these lines in dw (0-based).
            % For 3B2 imec data: the sync pulse is stored in line 6.
            % May be 1 or more line indices.
            dLineList = [0,1,6];
            
            if dataType == 'A'
                if strcmp(meta.typeThis, 'imec')
                    dataArray = MSpikeGLX.GainCorrectIM(dataArray, [ch], meta);
                else
                    dataArray = MSpikeGLX.GainCorrectNI(dataArray, [ch], meta);
                end
                plot(dataArray(ch,:));
            else
                digArray = MSpikeGLX.ExtractDigital(dataArray, meta, dw, dLineList);
                for i = 1:numel(dLineList)
                    plot(digArray(i,:));
                    hold on
                end
                hold off
            end
        end
        
        function dataArray = ReadBin(samp0, nSamp, meta, binName, path)
            % Read nSamp timepoints from the binary file, starting
            % at timepoint offset samp0. The returned array has
            % dimensions [nChan,nSamp]. Note that nSamp returned
            % is the lesser of: {nSamp, timepoints available}.
            %
            % IMPORTANT: samp0 and nSamp must be integers.
            %
            
            nChan = str2double(meta.nSavedChans);
            
            nFileSamp = str2double(meta.fileSizeBytes) / (2 * nChan);
            samp0 = max(samp0, 0);
            nSamp = min(nSamp, nFileSamp - samp0);
            
            sizeA = [nChan, nSamp];
            
            fid = fopen(fullfile(path, binName), 'rb');
            fseek(fid, samp0 * 2 * nChan, 'bof');
%             dataArray = fread(fid, sizeA, 'int16=>double');
            dataArray = fread(fid, sizeA, 'int16=>single');
            fclose(fid);
        end
        
        function digArray = ExtractDigital(dataArray, meta, dwReq, dLineList)
            % Return an array [timepoints X lines] of logical values for
            % a specified set of digital lines.
            %
            % - dwReq is the one-based index into the saved file of the
            %    16-bit word that contains the digital lines of interest.
            % - dLineList is a zero-based list of one or more lines/bits
            %    to scan from word dwReq.
            
            % Get channel index of requested digital word dwReq
            if strcmp(meta.typeThis, 'imec')
                [AP, LF, SY] = MSpikeGLX.ChannelCountsIM(meta);
                if SY == 0
                    fprintf('No imec sync channel saved\n');
                    digArray = [];
                    return;
                else
                    digCh = AP + LF + dwReq;
                end
            else
                [MN,MA,XA,DW] = MSpikeGLX.ChannelCountsNI(meta);
                if dwReq > DW
                    fprintf('Maximum digital word in file = %d\n', DW);
                    digArray = [];
                    return;
                else
                    digCh = MN + MA + XA + dwReq;
                end
            end
            nSamp = size(dataArray, 1);
            digArray = zeros(nSamp, numel(dLineList), 'logical');
            for i = 1:numel(dLineList)
                digArray(:,i) = bitget(int16(dataArray(:,digCh)), dLineList(i)+1);
            end
        end
        
        function srate = SampRate(meta)
            % Return sample rate as double.
            if strcmp(meta.typeThis, 'imec')
                srate = str2double(meta.imSampRate);
            else
                srate = str2double(meta.niSampRate);
            end
        end
        
        function fI2V = Int2Volts(meta)
            % Return a multiplicative factor for converting 16-bit
            % file data to voltage. This does not take gain into
            % account. The full conversion with gain is:
            %
            %   dataVolts = dataInt * fI2V / gain.
            %
            % Note that each channel may have its own gain.
            %
            if strcmp(meta.typeThis, 'imec')
                if isfield(meta,'imMaxInt')
                    maxInt = str2num(meta.imMaxInt);
                else
                    maxInt = 512;
                end
                fI2V = str2double(meta.imAiRangeMax) / maxInt;
            else
                fI2V = str2double(meta.niAiRangeMax) / 32768;
            end
        end
        
        function chans = OriginalChans(meta)
            % Return array of original channel IDs. As an example,
            % suppose we want the imec gain for the ith channel stored
            % in the binary data. A gain array can be obtained using
            % ChanGainsIM() but we need an original channel index to
            % do the look-up. Because you can selectively save channels
            % the ith channel in the file isn't necessarily the ith
            % acquired channel, so use this function to convert from
            % ith stored to original index.
            %
            % Note: In SpikeGLX channels are 0-based, but MATLAB uses
            % 1-based indexing, so we add 1 to the original IDs here.
            %
            if strcmp(meta.snsSaveChanSubset, 'all')
                chans = (1:str2double(meta.nSavedChans));
            else
                chans = str2num(meta.snsSaveChanSubset);
                chans = chans + 1;
            end
        end
        
        function [AP,LF,SY] = ChannelCountsIM(meta)
            % Return counts of each imec channel type that compose
            % the timepoints stored in binary file.
            %
            M = str2num(meta.snsApLfSy);
            AP = M(1);
            LF = M(2);
            SY = M(3);
        end
        
        function [MN,MA,XA,DW] = ChannelCountsNI(meta)
            % Return counts of each nidq channel type that compose
            % the timepoints stored in binary file.
            %
            M = str2num(meta.snsMnMaXaDw);
            MN = M(1);
            MA = M(2);
            XA = M(3);
            DW = M(4);
        end
        
        function gain = ChanGainNI(ichan, savedMN, savedMA, meta)
            % Return gain for ith channel stored in the nidq file.
            %
            % ichan is a saved channel index, rather than an original
            % (acquired) index.
            %
            if ichan <= savedMN
                gain = str2double(meta.niMNGain);
            elseif ichan <= savedMN + savedMA
                gain = str2double(meta.niMAGain);
            else
                gain = 1;
            end
        end
        
        function [APgain,LFgain] = ChanGainsIM(meta)
            % Return gain arrays for imec channels.
            %
            % Index into these with original (acquired) channel IDs.
            %
            
            if isfield(meta,'imDatPrb_type')
                probeType = str2num(meta.imDatPrb_type);
            else
                probeType = 0;
            end
            if (probeType == 21) || (probeType == 24)
                [AP,LF,~] = MSpikeGLX.ChannelCountsIM(meta);
                % NP 2.0; APgain = 80 for all channels
                APgain = zeros(AP,1,'double');
                APgain = APgain + 80;
                % No LF channels, set gain = 0
                LFgain = zeros(LF,1,'double');
            else
                % 3A or 3B data?
                % 3A metadata has field "typeEnabled" which was replaced
                % with "typeImEnabled" and "typeNiEnabled" in 3B.
                % The 3B imro table has an additional field for the
                % high pass filter enabled/disabled
                if isfield(meta,'typeEnabled')
                    % 3A data
                    C = textscan(meta.imroTbl, '(%*s %*s %*s %d %d', ...
                        'EndOfLine', ')', 'HeaderLines', 1 );
                else
                    % 3B data
                    C = textscan(meta.imroTbl, '(%*s %*s %*s %d %d %*s', ...
                        'EndOfLine', ')', 'HeaderLines', 1 );
                end
                APgain = double(cell2mat(C(1)));
                LFgain = double(cell2mat(C(2)));
            end
        end
        
        function dataArray = GainCorrectNI(dataArray, chanList, meta)
            % Having acquired a block of raw nidq data using ReadBin(),
            % convert values to gain-corrected voltages. The conversion
            % is only applied to the saved-channel indices in chanList.
            % Remember saved-channel indices are in range [1:nSavedChans].
            % The dimensions of the dataArray remain unchanged. ChanList
            % examples:
            %
            %   [1:MN]      % all MN chans (MN from ChannelCountsNI)
            %   [2,6,20]    % just these three channels
            %
            
            [MN,MA] = MSpikeGLX.ChannelCountsNI(meta);
            fI2V = MSpikeGLX.Int2Volts(meta);
            
            for i = 1:length(chanList)
                j = chanList(i);    % index into timepoint
                conv = fI2V / MSpikeGLX.ChanGainNI(j, MN, MA, meta);
                dataArray(j,:) = dataArray(j,:) * conv;
            end
        end
        
        function dataArray = GainCorrectIM(dataArray, chanList, meta)
            % Having acquired a block of raw imec data using ReadBin(),
            % convert values to gain-corrected voltages. The conversion
            % is only applied to the saved-channel indices in chanList.
            % Remember saved-channel indices are in range [1:nSavedChans].
            % The dimensions of the dataArray remain unchanged. ChanList
            % examples:
            %
            %   [1:AP]      % all AP chans (AP from ChannelCountsIM)
            %   [2,6,20]    % just these three channels
            %
            
            % Look up gain with acquired channel ID
            chans = MSpikeGLX.OriginalChans(meta);
            [APgain,LFgain] = MSpikeGLX.ChanGainsIM(meta);
            nAP = length(APgain);
            nNu = nAP * 2;
            
            % Common conversion factor
            fI2V = MSpikeGLX.Int2Volts(meta);
            
            for i = 1:length(chanList)
                j = chanList(i);    % index into timepoint
                k = chans(j);       % acquisition index
                if k <= nAP
                    conv = fI2V / APgain(k);
                elseif k <= nNu
                    conv = fI2V / LFgain(k - nAP);
                else
                    continue;
                end
                dataArray(j,:) = dataArray(j,:) * conv;
            end
        end
        
        % From SGLXMetaToCoords, git 140452d
        function [meta] = ReadMeta(metaName, path)
            % Parse ini file into cell entries C{1}{i} = C{2}{i}
            fid = fopen(fullfile(path, metaName), 'r');
            % -------------------------------------------------------------
            %    Need 'BufSize' adjustment for MATLAB earlier than 2014
            %    C = textscan(fid, '%[^=] = %[^\r\n]', 'BufSize', 32768);
            C = textscan(fid, '%[^=] = %[^\r\n]');
            % -------------------------------------------------------------
            fclose(fid);
            
            % New empty struct
            meta = struct();
            
            % Convert each cell entry into a struct entry
            for i = 1:length(C{1})
                tag = C{1}{i};
                if tag(1) == '~'
                    % remake tag excluding first character
                    tag = sprintf('%s', tag(2:end));
                end
                meta = setfield(meta, tag, C{2}{i});
            end
        end
        
        function [nShank, shankWidth, shankPitch, shankInd, xCoord, yCoord, connected] = geomMapToGeom(meta)
            % Parse snsGeomMap for XY coordinates
            
            C = textscan(meta.snsGeomMap, '(%d:%d:%d:%d', ...
                'EndOfLine', ')', 'HeaderLines', 1 );
            shankInd = double(cell2mat(C(1)));
            xCoord = double(cell2mat(C(2)));
            yCoord = double(cell2mat(C(3)));
            connected = double(cell2mat(C(4)));
            
            % parse header for number of shanks
            geomStr = meta.snsGeomMap;
            headStr = extractBefore(geomStr,')(');
            headParts = split(headStr,',');
            nShank = str2double(headParts{2});
            shankWidth = str2double(headParts{4});
            shankPitch = str2double(headParts{3});
        end
        
        function [nShank, shankWidth, shankPitch, shankInd, xCoord, yCoord, connected] = shankMapToGeom(meta)
            % Get XY coordinates from snsShankMap plus hard coded geom values
            
            % get number of saved AP channels (some early metadata files have a
            % SYNC entry in the snsChanMap
            [nchan,~,~] = MSpikeGLX.ChannelCountsIM(meta);
            
            C = textscan(meta.snsShankMap, '(%d:%d:%d:%d', ...
                'EndOfLine', ')', 'HeaderLines', 1 );
            shankInd = double(cell2mat(C(1)));
            colInd = double(cell2mat(C(2)));
            rowInd = double(cell2mat(C(3)));
            connected = double(cell2mat(C(4)));
            
            % trim these to the number of saved channels
            shankInd = shankInd(1:nchan);
            colInd = colInd(1:nchan);
            rowInd = rowInd(1:nchan);
            connected = connected(1:nchan);
            
            geom = MSpikeGLX.getGeomParams(meta);
            
            oddRows = logical(mod(rowInd,2));
            evenRows = ~oddRows;
            xCoord = colInd*geom.horzPitch;
            xCoord(evenRows) = xCoord(evenRows) + geom.even_xOff ;
            xCoord(oddRows) = xCoord(oddRows) + geom.odd_xOff;
            yCoord = rowInd*geom.vertPitch;
            
            nShank = geom.nShank;
            shankWidth = geom.shankWidth;
            shankPitch = geom.shankPitch;
        end
        
        function geom = getGeomParams(meta)
            % Return geometry paramters for supported probe types
            % These are used to calculate positions from metadata
            % that includes only ~snsShankMap
            %
            
            % create map
            geomTypeMap = MSpikeGLX.makeTypeMap();
            
            % get probe part number; if absent, this is a 3A
            if isfield(meta,'imDatPrb_pn')
                pn = meta.imDatPrb_pn;
            else
                pn = '3A';
            end
            
            if geomTypeMap.isKey(pn)
                geomType = geomTypeMap(pn);
            else
                fprintf('unsupported probe part number\n');
                return;
            end
            
            switch geomType
                case 'np1_stag_70um'
                    geom.nShank = 1;
                    geom.shankWidth = 70;
                    geom.shankPitch = 0;
                    geom.even_xOff = 27;
                    geom.odd_xOff = 11;
                    geom.horzPitch = 32;
                    geom.vertPitch = 20;
                    geom.rowsPerShank = 480;
                    geom.elecPerShank = 960;
                case 'nhp_lin_70um'
                    geom.nShank = 1;
                    geom.shankWidth = 70;
                    geom.shankPitch = 0;
                    geom.even_xOff = 27;
                    geom.odd_xOff = 27;
                    geom.horzPitch = 32;
                    geom.vertPitch = 20;
                    geom.rowsPerShank = 480;
                    geom.elecPerShank = 960;
                case 'nhp_stag_125um_med'
                    geom.nShank = 1;
                    geom.shankWidth = 125;
                    geom.shankPitch = 0;
                    geom.even_xOff = 27;
                    geom.odd_xOff = 11;
                    geom.horzPitch = 87;
                    geom.vertPitch = 20;
                    geom.rowsPerShank = 1368;
                    geom.elecPerShank = 2496;
                case 'nhp_stag_125um_long'
                    geom.nShank = 1;
                    geom.shankWidth = 125;
                    geom.shankPitch = 0;
                    geom.even_xOff = 27;
                    geom.odd_xOff = 11;
                    geom.horzPitch = 87;
                    geom.vertPitch = 20;
                    geom.rowsPerShank = 2208;
                    geom.elecPerShank = 4416;
                case 'nhp_lin_125um_med'
                    geom.nShank = 1;
                    geom.shankWidth = 125;
                    geom.shankPitch = 0;
                    geom.even_xOff = 11;
                    geom.odd_xOff = 11;
                    geom.horzPitch = 103;
                    geom.vertPitch = 20;
                    geom.rowsPerShank = 1368;
                    geom.elecPerShank = 2496;
                case 'nhp_lin_125um_long'
                    geom.nShank = 1;
                    geom.shankWidth = 125;
                    geom.shankPitch = 0;
                    geom.even_xOff = 11;
                    geom.odd_xOff = 11;
                    geom.horzPitch = 103;
                    geom.vertPitch = 20;
                    geom.rowsPerShank = 2208;
                    geom.elecPerShank = 4416;
                case 'uhd_8col_1bank'
                    geom.nShank = 1;
                    geom.shankWidth = 70;
                    geom.shankPitch = 0;
                    geom.even_xOff = 14;
                    geom.odd_xOff = 14;
                    geom.horzPitch = 6;
                    geom.vertPitch = 6;
                    geom.rowsPerShank = 48;
                    geom.elecPerShank = 384;
                case 'uhd_8col_16bank'
                    geom.nShank = 1;
                    geom.shankWidth = 70;
                    geom.shankPitch = 0;
                    geom.even_xOff = 14;
                    geom.odd_xOff = 14;
                    geom.horzPitch = 6;
                    geom.vertPitch = 6;
                    geom.rowsPerShank = 768;
                    geom.elecPerShank = 6144;
                case 'np2_ss'
                    geom.nShank = 1;
                    geom.shankWidth = 70;
                    geom.shankPitch = 0;
                    geom.even_xOff = 27;
                    geom.odd_xOff = 27;
                    geom.horzPitch = 32;
                    geom.vertPitch = 15;
                    geom.rowsPerShank = 640;
                    geom.elecPerShank = 1280;
                case 'np2_4s'
                    geom.nShank = 4;
                    geom.shankWidth = 70;
                    geom.shankPitch = 250;
                    geom.even_xOff = 27;
                    geom.odd_xOff = 27;
                    geom.horzPitch = 32;
                    geom.vertPitch = 15;
                    geom.rowsPerShank = 640;
                    geom.elecPerShank = 1280;
                case 'NP1120'
                    geom.nShank = 1;
                    geom.shankWidth = 70;
                    geom.shankPitch = 0;
                    geom.even_xOff = 6.75;
                    geom.odd_xOff = 6.75;
                    geom.horzPitch = 4.5;
                    geom.vertPitch = 4.5;
                    geom.rowsPerShank = 192;
                    geom.elecPerShank = 384;
                case 'NP1121'
                    geom.nShank = 1;
                    geom.shankWidth = 70;
                    geom.shankPitch = 0;
                    geom.even_xOff = 6.25;
                    geom.odd_xOff = 6.25;
                    geom.horzPitch = 3;
                    geom.vertPitch = 3;
                    geom.rowsPerShank = 384;
                    geom.elecPerShank = 384;
                case 'NP1122'
                    geom.nShank = 1;
                    geom.shankWidth = 70;
                    geom.shankPitch = 0;
                    geom.even_xOff = 12.5;
                    geom.odd_xOff = 12.5;
                    geom.horzPitch = 3;
                    geom.vertPitch = 3;
                    geom.rowsPerShank = 24;
                    geom.elecPerShank = 384;
                case 'NP1123'
                    geom.nShank = 1;
                    geom.shankWidth = 70;
                    geom.shankPitch = 0;
                    geom.even_xOff = 10.25;
                    geom.odd_xOff = 10.25;
                    geom.horzPitch = 4.5;
                    geom.vertPitch = 4.5;
                    geom.rowsPerShank = 32;
                    geom.elecPerShank = 384;
                case 'NP1300'
                    geom.nShank = 1;
                    geom.shankWidth = 70;
                    geom.shankPitch = 0;
                    geom.even_xOff = 11;
                    geom.odd_xOff = 11;
                    geom.horzPitch = 48;
                    geom.vertPitch = 20;
                    geom.rowsPerShank = 480;
                    geom.elecPerShank = 960;
                case 'NP1200'
                    geom.nShank = 1;
                    geom.shankWidth = 70;
                    geom.shankPitch = 0;
                    geom.even_xOff = 27;
                    geom.odd_xOff = 11;
                    geom.horzPitch = 32;
                    geom.vertPitch = 20;
                    geom.rowsPerShank = 64;
                    geom.elecPerShank = 128;
                case 'NXT3000'
                    geom.nShank = 1;
                    geom.shankWidth = 70;
                    geom.shankPitch = 0;
                    geom.even_xOff = 53;
                    geom.odd_xOff = 53;
                    geom.horzPitch = 0;
                    geom.vertPitch = 15;
                    geom.rowsPerShank = 128;
                    geom.elecPerShank = 128;
                otherwise
                    % shouldn't see this case
                    fprintf('unsupported probe part number\n');
                    return;
            end
        end
        
        function M = makeTypeMap()
            % Return geometry paramters for supported probe types
            % Note that geom only contains enough info to calculate
            % positions for the electrodes listed in snsShankMap
            %
            
            % many part numbers have the same geometry parameters ;
            % make a map that pairs geometry type (value) with probe part number (key)
            M = containers.Map('KeyType','char','ValueType','char');
            
            M('3A') = 'np1_stag_70um';
            M('PRB_1_4_0480_1') = 'np1_stag_70um';
            M('PRB_1_4_0480_1_C') = 'np1_stag_70um';
            M('NP1010') = 'np1_stag_70um';
            M('NP1011') = 'np1_stag_70um';
            M('NP1012') = 'np1_stag_70um';
            M('NP1013') = 'np1_stag_70um';
            
            M('NP1015') = 'nhp_lin_70um';
            M('NP1015') = 'nhp_lin_70um';
            M('NP1016') = 'nhp_lin_70um';
            M('NP1017') = 'nhp_lin_70um';
            
            M('NP1020') = 'nhp_stag_125um_med';
            M('NP1021') = 'nhp_stag_125um_med';
            M('NP1030') = 'nhp_stag_125um_long';
            M('NP1031') = 'nhp_stag_125um_long';
            
            M('NP1022') = 'nhp_lin_125um_med';
            M('NP1032') = 'nhp_lin_125um_long';
            
            M('NP1100') = 'uhd_8col_1bank';
            M('NP1110') = 'uhd_8col_16bank';
            
            M('PRB2_1_2_0640_0') = 'np2_ss';
            M('PRB2_1_4_0480_1') = 'np2_ss';
            M('NP2000') = 'np2_ss';
            M('NP2003') = 'np2_ss';
            M('NP2004') = 'np2_ss';
            
            M('PRB2_4_2_0640_0') = 'np2_4s';
            M('PRB2_4_4_0480_1') = 'np2_4s';
            M('NP2010') = 'np2_4s';
            M('NP2013') = 'np2_4s';
            M('NP2014') = 'np2_4s';
            
            M('NP1120') = 'NP1120';
            M('NP1121') = 'NP1121';
            M('NP1122') = 'NP1122';
            M('NP1123') = 'NP1123';
            M('NP1300') = 'NP1300';
            
            M('NP1200') = 'NP1200';
            M('NXT3000') = 'NXT3000';
        end
        
        % From Neuropixel-utils by Daniel J. O'Shea et al.
        function [pathRoot, fileStem, type, imecNumber] = ParseImecFileName(file)
            if iscell(file)
                [pathRoot, fileStem, type] = cellfun(@MSpikeGLX.ParseImecFileName, file, 'UniformOutput', false);
                return;
            end
            file = char(file);
            
            [pathRoot, f, e] = fileparts(file);
            if isempty(e)
                error('No file extension specified on Imec file name');
            end
            file = [f, e];
            
            match = regexp(file, '(?<stem>[\w\-\.]+).imec(?<imecNumber>\d*).(?<type>\w+).bin', 'names', 'once');
            if ~isempty(match)
                type = match.type;
                fileStem = match.stem;
                if isempty(match.imecNumber)
                    imecNumber = NaN;
                else
                    imecNumber = str2double(match.imecNumber);
                end
                return;
            end
            
            fileStem = file;
            type = '';
            imecNumber = NaN;
        end
        
    end
    
end

