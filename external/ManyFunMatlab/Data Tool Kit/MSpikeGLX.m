classdef MSpikeGLX
    % Methods in this class are mostly copied from DemoReadSGLXData.m and SGLXMetaToCoords.m 
    % with the following changes
    %   1) In ReadBin, convert int16 to single instead of double
    %   2) Add high-level reader methods: ReadNI, ReadLFP, ReadMetaEZ
    %   3) Add helper methods: GetCoordNP10
    
    methods(Static)
        function meta = ReadMetaEZ(metaPath)
            % Read metadata with an optional single path input
            if nargin < 1 || isempty(metaPath)
                metaPath = MBrowse.File();
            end
            [folderPath, baseName] = fileparts(metaPath);
            meta = MSpikeGLX.ReadMeta([baseName '.meta'], folderPath);
        end
        
        function [meta, aiArray, diArray, t] = ReadNI(metaPath)
            % 
            
            % Read metadata
            [folderPath, baseName] = fileparts(metaPath);
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
        
        function [meta, lfpArray, t] = ReadLFP(metaPath)
            % 
            
            % Read metadata
            [folderPath, baseName] = fileparts(metaPath);
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
        
        function tb = GetCoordNP10(apMeta)
            % Make a table of positional info
            tb = table;
            tb.chanId = (1:384)';
            [tb.elecInd, tb.connected] = MSpikeGLX.NP10_ElecInd(apMeta);
            [tb.xCoord, tb.yCoord] = MSpikeGLX.XYCoord10(apMeta, tb.elecInd);
        end
        
        % From DemoReadSGLXData.m
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
        
        % From SGLXMetaToCoords.m
        function SGLXMetaToCoords()
            % Write out coordinates for a Neuropixels 3A, 1.0 or 2.0 metadata file.
            % Format selected with the outType variable.
            % Jennifer Colonell, Janelia Research Campus
            %
            % Output selection: 0 for text coordinate file;
            %                   1 for Kilosort or Kilosort2 channel map file;
            %                   2 for strings to paste into JRClust .prm file
            outType = 1;
            
            % Ask user for metadata file
            [metaName,path] = uigetfile('*.meta', 'Select Metadata File');
            
            % Shank separation for multishank
            shankSep = 250;
            
            % Parse in file to get the metadata structure
            meta = ReadMeta(metaName, path);
            
            if isfield(meta, 'imDatPrb_type')
                pType = str2num(meta.imDatPrb_type);
            else
                pType = 0; %3A probe
            end
            
            if pType <= 1
                
                %Neuropixels 1.0 or 3A probe
                [elecInd, connected] = MSpikeGLX.NP10_ElecInd(meta);
                
                % Get saved channels
                chans = MSpikeGLX.OriginalChans(meta);
                [AP,LF,SY] = MSpikeGLX.ChannelCountsIM(meta);
                chans = chans(1:AP);        %1-based channel numbers
                
                % Trim elecInd and shankInd to include only saved channels
                elecInd = elecInd(chans);
                shankind = zeros(size(elecInd));
                
                % Get XY coords for saved channels
                [xcoords, ycoords] = MSpikeGLX.XYCoord10(meta, elecInd);
                
            else
                
                % Parse imro table for shank and electrode indicies
                [elecInd, shankind, bankMask, connected] = MSpikeGLX.NP20_ElecInd(meta);
                
                % Get saved channels
                chans = MSpikeGLX.OriginalChans(meta);
                [AP,LF,SY] = MSpikeGLX.ChannelCountsIM(meta);
                chans = chans(1:AP);        %1-based channel numbers
                
                % Trim elecInd and shankInd to include only saved channels
                elecInd = elecInd(chans);
                shankind = shankind(chans);
                
                % Get XY coords for saved channels
                [xcoords, ycoords] = MSpikeGLX.XYCoord20(meta, elecInd, bankMask, shankind);
                
            end
            
            % Build output name and write out file
            [~,fname,~] = fileparts(metaName);
            
            switch outType
                case 0      %tab delimited, chan, x, y, shank
                    newName = [fname,'-siteCoords.txt'];
                    fid = fopen( newName, 'w');
                    for i = 1:numel(elecInd)
                        currX = shankind(i)*shankSep + xcoords(i);
                        fprintf( fid, '%d\t%d\t%d\t%d\n', chans(i)-1, currX, ycoords(i), shankind(i));
                    end
                    fclose(fid);
                    
                case 1     %KS2 *.mat
                    newName = [fname,'_kilosortChanMap.mat'];
                    chanMap = (1:numel(chans))';
                    chanMap0ind = chanMap - 1;
                    connected = logical(connected);
                    xcoords = shankind*shankSep + xcoords;   %KS2 not yet using kcoords, so x coord includes shank sep
                    kcoords = shankind + 1;     %KS1 uses kcoords to force templates to be on one shank
                    name = fname;
                    save( newName, 'chanMap', 'chanMap0ind', 'connected', 'name', 'xcoords', 'ycoords', 'kcoords' );
                    
                case 2  %strings to copy into JRC prm file
                    newName = [fname,'_forJRCprm.txt'];
                    nchan = numel(chans);
                    fid = fopen( newName, 'w' );
                    fprintf( fid, 'shankMap = [' );
                    for i = 1:nchan-1
                        fprintf( fid, '%d,', shankind(i) + 1 ); % switch to 1-based for MATLAB
                    end
                    fprintf( fid, '%d];\n',shankind(nchan) + 1 );
                    
                    xcoords = shankind*shankSep + xcoords;
                    
                    fprintf( fid, 'siteLoc = [' );
                    for i = 1:nchan-1
                        fprintf(fid, '%d,%d;', xcoords(i), ycoords(i));
                    end
                    fprintf( fid, '%d,%d];\n', xcoords(nchan), ycoords(nchan) );
                    
                    fprintf( fid, 'siteMap = [' );
                    for i = 1:nchan-1
                        fprintf( fid, '%d,', chans(i) );
                    end
                    fprintf( fid, '%d];\n', chans(nchan) );
                    fclose(fid);
            end
        end
        
        function [meta] = ReadMeta(metaName, path)
            % Parse ini file returning a structure whose field names
            % are the metadata left-hand-side tags, and whose right-
            % hand-side values are MATLAB strings. We remove any
            % leading '~' characters from tags because MATLAB uses
            % '~' as an operator.
            %
            % If you're unfamiliar with structures, the benefit
            % is that after calling this function you can refer
            % to metafile items by name. For example:
            %
            %   meta.fileCreateTime  // file create date and time
            %   meta.nSavedChans     // channels per timepoint
            %
            % All of the values are MATLAB strings, but you can
            % obtain a numeric value using str2double(meta.nSavedChans).
            % More complicated parsing of values is demonstrated in the
            % utility functions below.
            %
            
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
        
        function [elecInd, shankInd, bankMask, connected] = NP20_ElecInd(meta)
            % Return shank and electrode number for NP2.0.
            %
            % Index into these with original (acquired) channel IDs.
            %
            pType = str2num(meta.imDatPrb_type);
            if pType == 21
                % Single shank probe
                % imro table entries: (channel, bank, refType, electrode #)
                C = textscan(meta.imroTbl, '(%*s %d %*s %d', ...
                    'EndOfLine', ')', 'HeaderLines', 1 );
                elecInd = int32(cell2mat(C(2)));
                bankMask = int32(cell2mat(C(1)));
                shankInd = zeros(size(elecInd));
                connected = ones(size(elecInd));
                exChan = MSpikeGLX.FindDisabled(meta);
                for i = 1:numel(exChan)
                    connected(elecInd == exChan(i)) = 0;
                end
                
            else
                % 4 shank probe
                % imro table entries: (channel, shank, bank, refType, electrode #)
                C = textscan(meta.imroTbl, '(%d %d %d %*s %d', ...
                    'EndOfLine', ')', 'HeaderLines', 1 );
                chan = double(cell2mat(C(1)));
                elecInd = int32(cell2mat(C(4)));
                bankMask = int32(cell2mat(C(3)));
                shankInd = double(cell2mat(C(2)));
                connected = ones(size(chan));
                exChan = MSpikeGLX.FindDisabled(meta);
                %exChan = [127];
                for i = 1:numel(exChan)
                    connected(chan == exChan(i)) = 0;
                end
            end
        end
        
        function [elecInd, connected] = NP10_ElecInd(meta)
            % Return shank and electrode number for NP1.0.
            %
            % Index into these with original (acquired) channel IDs.
            %
            % 3A or 3B data?
            % 3A metadata has field "typeEnabled" which was replaced
            % with "typeImEnabled" and "typeNiEnabled" in 3B.
            % The 3B imro table has an additional field for the
            % high pass filter enabled/disabled
            % Note that the textscan funtion places line breaks at each
            % instance of the 'EndofLine' character -- here, ')'
            % 'HeaderLines' = 1 skips the initial entry in the table with
            % the probe type and number of entries.
            
            if isfield(meta,'typeEnabled')
                % 3A data
                C = textscan(meta.imroTbl, '(%d %d %*s %*s %*s', ...
                    'EndOfLine', ')', 'HeaderLines', 1 );
                exChan = MSpikeGLX.FindDisabled(meta);
                %exChan = [36, 75, 112, 151, 188, 227, 264, 303, 340, 373];
            else
                % 3B data
                C = textscan(meta.imroTbl, '(%d %d %*s %*s %*s %*s', ...
                    'EndOfLine', ')', 'HeaderLines', 1 );
                exChan = MSpikeGLX.FindDisabled(meta);
                %exChan = [191];
            end
            chan = double(cell2mat(C(1)));
            bank = double(cell2mat(C(2)));
            elecInd = bank*384 + chan;
            connected = ones(size(chan));
            for i = 1:numel(exChan)
                connected(chan == exChan(i)) = 0;
            end
        end
        
        function [exChan] = FindDisabled(meta)
            % Read shank map for any probe type and return list
            % of channels that are disabled. This will include the
            % reference channels
            %
            % Note that the textscan funtion places line breaks at each
            % instance of the 'EndofLine' character -- here, ')'
            % 'HeaderLines' = 1 skips the initial entry in the table with
            % the number of shanks, columns, and rows
            
            % read in the shank map
            C = textscan(meta.snsShankMap, '(%d:%d:%d:%d', ...
                'EndOfLine', ')', 'HeaderLines', 1 );
            enabled = double(cell2mat(C(4)));
            % There's an entry in the shank map for each saved channel.
            % Get the array of saved channels:
            chan = MSpikeGLX.OriginalChans(meta);
            % Find out how many are non-SY chans
            [AP,~,~] = MSpikeGLX.ChannelCountsIM(meta);
            exChan = [];
            for i = 1:AP
                if enabled(i) == 0
                    exChan = [exChan, chan(i)];
                end
            end
        end
        
        function [xCoord, yCoord] = XYCoord20(meta, elecInd, bankMask, shankind)
            % Return x y coords for electrode index for 2.0 probes
            
            pType = str2num(meta.imDatPrb_type);
            
            nElec = 1280;   %per shank; pattern repeats for the four shanks
            vSep = 15;   % in um
            hSep = 32;
            
            elecPos = zeros(nElec, 2);
            
            elecPos(1:2:end,1) = 0;         %sites 0,2,4...
            elecPos(2:2:end,1) = hSep;      %sites 1,3,5...
            
            % fill in y values
            viHalf = (0:(nElec/2-1))';                %row numbers
            elecPos(1:2:end,2) = viHalf * vSep;       %sites 0,2,4...
            elecPos(2:2:end,2) = elecPos(1:2:end,2);  %sites 1,3,5...
            
            xCoord = elecPos(elecInd+1,1);
            yCoord = elecPos(elecInd+1,2);
            
            if pType == 21
                % single shank probe. Plot only lowest selected electrode
                figure(1)
                % plot all positions
                scatter( elecPos(:,1), elecPos(:,2), 150, 'k', 'square' ); hold on;
                scatter( xCoord, yCoord, 100, 'b', 'square', 'filled' );hold on;
                xlim([-16,64]);
                ylim([-10,10000]);
                title('NP 2.0 single shank view');
                hold off;
            else
                % four shank probe, no multiple connections
                figure(1)
                shankSep = 250;
                for sI = 0:3
                    cc = find(shankind == sI);
                    scatter( shankSep*sI + elecPos(:,1), elecPos(:,2), 30, 'k', 'square' ); hold on;
                    scatter( shankSep*sI + xCoord(cc), yCoord(cc), 20, 'b', 'square', 'filled' ); hold on;
                end
                xlim([-16,3*shankSep+64]);
                ylim([-10,10000]);
                title('NP2.0 MS shank view');
                hold off;
            end
        end
        
        function [xCoord, yCoord] = XYCoord10(meta, elecInd)
            % Return x y coords for electrode index for 1.0 probes
            
            nElec = 960;   %per shank; pattern repeats for the four shanks
            vSep = 20;   % in um
            hSep = 32;
            
            elecPos = zeros(nElec, 2);
            
            elecPos(1:4:end,1) = hSep/2;            %sites 0,4,8...
            elecPos(2:4:end,1) =  (3/2)*hSep;       %sites 1,5,9...
            elecPos(3:4:end,1) = 0;                 %sites 2,6,10...
            elecPos(4:4:end,1) =  hSep;             %sites 3,7,11...
            elecPos(:,1) = elecPos(:,1) + 11;       %x offset on the shank
            
            % fill in y values
            viHalf = (0:(nElec/2-1))';                %row numbers
            elecPos(1:2:end,2) = viHalf * vSep;       %sites 0,2,4...
            elecPos(2:2:end,2) = elecPos(1:2:end,2);  %sites 1,3,5...
            
            xCoord = elecPos(elecInd+1,1);
            yCoord = elecPos(elecInd+1,2);
            
            % single shank probe. Plot only lowest selected electrode
            figure(1)
            % plot all positions
            scatter( elecPos(:,1), elecPos(:,2), 150, 'k', 'square' ); hold on;
            scatter( xCoord, yCoord, 100, 'b', 'square', 'filled' );hold on;
            xlim([0,70]);
            ylim([-10,8000]);
            title('NP 1.0 single shank view');
            hold off;
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

