classdef MIntan
    %MDATA Summary of this class goes here
    %   Detailed explanation goes here
    
    methods(Static)
        function ops = GetOptions(signalList)
            % Get option structures to use with MIntan.ReadPhdFiles
            %
            %   ops = MIntan.GetOptions()
            %   ops = MIntan.GetOptions(signalList)
            %
            % Input
            %   signalList          A string or a cell array of strings specifying the signalName field (see below) of 
            %                       each output struct. Possible options are 'amplifier', 'aux_in', 'adc', 'dig_in'. 
            %                       The default is {'amplifier', 'aux_in', 'adc', 'dig_in'} which creates a four-element
            %                       structure array. 
            % Output
            %   ops                 Structure(s) of options used to specify the processig in MIntan.ReadRhdFiles. Each 
            %                       structure is independent from any others.
            %     signalName        The name of signal source to read data from.
            %     signalFunc        A handle to a function that operates on the signal read from each file. The input to
            %                       this function is an [#channel,#sample] array. The output size should be consistent
            %                       with timestamps. The default is empty which does nothing.
            %     downsampleFactor  A factor to downsample the signal with median filter. For example, a factor of 30
            %                       can downsample a 30kHz signal to 1kHz using a 1-by-30 median filter. Downsampling
            %                       is performed across RHD files when samlpes are incomplete at file edges. The default
            %                       is 1 which does nothing.
            %     isReturn          Whether or not this signal should be returned in the output of MIntan.ReadRhdFiles.
            %                       The default is true.
            %     varBaseName       Variable base name used to store processed signal in the output of MIntan.ReadRhdFiles.
            %                       Timestamps and data will be specified by '_time' and '_data' suffix to the base name.
            %                       The default is simply the signalName.
            %     binFilePath       If specified, the signal data will be saved as a binary file. The data type will
            %                       inherit from the data array. However, you can use signalFunc to change data type.
            %                       Timestamps are not saved.
            %
            
            % Handle user inputs
            if nargin < 1
                signalList = {'amplifier', 'aux_in', 'adc', 'dig_in'};
            end
            signalList = cellstr(signalList);
            
            % Create option structs
            for i = numel(signalList) : -1 : 1
                ops(i) = signalOptions(signalList{i});
            end
            
            function s = signalOptions(signalName)
                % Signal processing
                s.signalName = signalName;
                s.signalFunc = [];
                s.downsampleFactor = 1;
                
                % Returning data in function output
                s.isReturn = true;
                s.varBaseName = signalName;
                
                % Generation of binary files
                s.binFilePath = '';
            end
        end
        
        function result = ReadRhdFiles(rhdFilePaths, ops)
            % Read Intan RHD2000 files (a fancy wrapper of Intan.read_Intan_RHD2000_file_noUI)
            %
            %   result = MIntan.ReadRhdFiles()
            %   result = MIntan.ReadRhdFiles(rhdFilePaths)
            %   result = MIntan.ReadRhdFiles(rhdFilePaths, ops)
            %
            % Inputs
            %   rhdFilePaths        A string or a cell array of strings of RHD2000 file paths. If left empty,
            %                       a file selection window will show up.
            %   ops                 Structure(s) returned from MIntan.GetOptions method and optionally
            %                       modified by user to customize loading. See MIntan.GetOptions for details.
            %                       The default is extracting all data with no custom preprocessing.
            % Output
            %   result              A structure with read/processed data and metadata.
            %
            
            warning('off', 'backtrace');
            
            % Handles user inputs
            if nargin < 2
                ops = MIntan.GetOptions();
            end
            if nargin < 1 || isempty(rhdFilePaths)
                rhdFilePaths = MBrowse.Files();
                if isempty(rhdFilePaths)
                    result = [];
                    return;
                end
            end
            rhdFilePaths = cellstr(rhdFilePaths);
            
            % Preallocation
            timeCells = cell(numel(rhdFilePaths), numel(ops));
            dataCells = cell(numel(rhdFilePaths), numel(ops));
            carryoverTime = cell(1, numel(ops));
            carryoverData = cell(1, numel(ops));
            
            for k = 1 : numel(rhdFilePaths)
                % Load data from a Intan file
                [notes, frequency_parameters, ...
                    amplifier_channels, amplifier_data, amplifier_time, ...
                    aux_in_channels, aux_in_data, aux_in_time, ...
                    adc_channels, adc_data, adc_time, ...
                    dig_in_channels, dig_in_data, dig_in_time] = ...
                    Intan.read_Intan_RHD2000_file_noUI(rhdFilePaths{k});
                
                amplifier_data = single(amplifier_data);
                aux_in_data = single(aux_in_data);
                adc_data = single(adc_data);
                dig_in_data = uint8(dig_in_data);
                
                fprintf('\n');
                
                % Process requested signals
                for i = 1 : numel(ops)
                    % Check and cache variables
                    sigName = ops(i).signalName;
                    assert(ismember(sigName, {'amplifier', 'adc', 'dig_in', 'aux_in'}), ...
                        '%s (case-sensitive) is not a valid signal name', sigName);
                    
                    sigData = eval([sigName '_data;']);
                    sigTime = eval([sigName '_time;']);
                    
                    if isempty(sigData)
                        warning('%s data is not available', sigName);
                        continue;
                    else
                        fprintf('Processing %s data\n', sigName);
                    end
                    
                    % Apply user functions
                    if isa(ops(i).signalFunc, 'function_handle')
                        fprintf('- applying user''s signal function\n');
                        sigData = ops(i).signalFunc(sigData);
                    end
                    
                    % Downsampling
                    dsF = ops(i).downsampleFactor;
                    if dsF > 1
                        fprintf('- downsampling %d times with 1-by-%d median filter\n', dsF, dsF);
                        
                        % Combine carry-over data
                        sigData = [carryoverData{i} sigData];
                        sigTime = [carryoverTime{i} sigTime];
                        
                        % Carry remainder data over
                        nSp = numel(sigTime);
                        r = mod(nSp, dsF);
                        carryoverInd = nSp-r+1 : nSp;
                        carryoverData{i} = sigData(:,carryoverInd);
                        carryoverTime{i} = sigTime(carryoverInd);
                        
                        if ~isempty(carryoverInd)
                            if k < numel(rhdFilePaths)
                                % Trim off the carry-over data
                                sigData(:,carryoverInd) = [];
                                sigTime(carryoverInd) = [];
                            else
                                % Pad remainder data for computing the last median
                                sigPad = repmat(median(carryoverData{i},2), [1 dsF-r]);
                                sigData = [sigData sigPad];
                            end
                        end
                        
                        % Downsample data
                        sigData = reshape(sigData, [size(sigData,1), dsF, size(sigData,2)/dsF]);
                        sigData = median(sigData, 2);
                        sigData = permute(sigData, [1 3 2]);
                        sigTime = downsample(sigTime, dsF);
                    end
                    
                    % Store processed data
                    if ops(i).isReturn
                        dataCells{k,i} = sigData';
                        timeCells{k,i} = sigTime';
                    end
                    
                    % Save processed data to binary file
                    bPath = ops(i).binFilePath;
                    if isempty(bPath)
                        continue;
                    end
                    
                    bDir = fileparts(bPath);
                    if ~isempty(bDir) && ~exist(bDir, 'dir')
                        mkdir(bDir);
                    end
                    
                    if k == 1 && exist(bPath, 'file')
                        warning('The binary file already exists and will not be generated again. \n%s', bPath);
                        ops(i).binFilePath = [];
                        continue;
                    end
                    
                    fprintf('- writing data to binary file\n');
                    fid = fopen(bPath, 'a');
                    fwrite(fid, sigData, class(sigData));
                    fclose(fid);
                end
            end
            
            % Combine data chunks
            for i = 1 : numel(ops)
                if ops(i).isReturn
                    result.([ops(i).varBaseName '_data']) = cell2mat(dataCells(:,i));
                    result.([ops(i).varBaseName '_time']) = cell2mat(timeCells(:,i));
                end
            end
            
            % Add metadata
            result.info.file_paths = rhdFilePaths;
            result.info.frequency_parameters = frequency_parameters;
            result.info.amplifier_channels = amplifier_channels;
            result.info.aux_in_channels = aux_in_channels;
            result.info.adc_channels = adc_channels;
            result.info.dig_in_channels = dig_in_channels;
            result.info.ops = ops;
            
            warning('on', 'backtrace');
        end
    end
end


