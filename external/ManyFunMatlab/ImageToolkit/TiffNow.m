classdef TiffNow < matlab.mixin.Copyable
    %TiffNow Summary of this class goes here
    %   
    % Loading TIFF file in a dynamic way allowing for instant retreival and memory management
    %
    %   tiffNowObject = TiffNow(path)
    %   
    % Inputs
    %    path: The path of the TIFF file. File browsing UI will prompt out if passing empty value.
    %   'dataType': The data type of the image. (default 'uint8' or 'uint16', otherwise 'double')
    %   'maxBytes': Maximal memory usage in byte. (default is the current system available memory)
    %   'framePerBlock': Number of frames per loading block. The smaller the faster response but longer loading time. (default 10)
    %   'loadingInterval': Interval (in second) between loading two blocks. The longer the faster response but longer loading time. (default 0.05)
    %   'verbose': A boolean value. Whether to show internal message. 
    %
    % Output
    %   TiffNow object: can almost be used and indexed as a numeric array (keyword 'end' is not supported)
    %   
    %   
    
    properties(SetAccess = private)
        filePath;
        imgSize;
    end
    
    properties(Access = private)
        tifObj;
        dataType;
        
        blocks;
        numBlocks;
        frPerBlock;
        maxNumBlocks;
        
        timerObj;
        queue;
        
        verbose;
    end
    
    methods
        % Constructor
        function this = TiffNow(path, varargin)
            p = inputParser;
            p.addParameter('dataType', '', @ischar);
            p.addParameter('maxBytes', inf, @isnumeric);
            p.addParameter('loadingInterval', 0.05, @isnumeric);
            p.addParameter('framePerBlock', 10, @isnumeric);
            p.addParameter('verbose', false, @islogical);
            p.parse(varargin{:});
            this.dataType = p.Results.dataType;
            maxBytes = p.Results.maxBytes;
            loadingInterval = p.Results.loadingInterval;
            this.frPerBlock = p.Results.framePerBlock;
            this.verbose = p.Results.verbose;
            
            if isempty(path)
                path = Browse.File();
            end
            this.filePath = path;
            
            try
                % Creates tiflib object
                [ ~, ~, ext ] = fileparts(this.filePath);
                if ~strcmp(ext, '.tif')
                    this.filePath = [ this.filePath, '.tif' ];
                end
                this.tifObj = Tiff(this.filePath, 'r');
                
                % Gets image dimensions
                [ this.imgSize(1), this.imgSize(2), this.imgSize(3) ] = Img23.GetTiffImageSize(this.tifObj);
                
                % Gets datatype if not specified by user
                if isempty(this.dataType)
                    pxlBytes = round(this.tifObj.getTag('BitsPerSample') / 8);
                    switch pxlBytes
                        case 1
                            this.dataType = 'uint8';
                        case 2
                            this.dataType = 'uint16';
                        otherwise
                            this.dataType = 'double';
                    end
                end
            catch
                this.imgSize = [];
            end
            
            if isempty(this.imgSize)
                error('bad TIFF file...');
            else
                % Initializes image blocks
                % Gets all frame indices
                allInd = (1 : this.imgSize(3))';
                
                % Creates the modulus block
                modNum = mod(this.imgSize(3), this.frPerBlock);
                if modNum
                    modBlockInd = { allInd(end-modNum+1:end) };
                else
                    modBlockInd = {};
                end
                
                % Creates full blocks
                if this.imgSize(3) >= this.frPerBlock
                    fullBlockInd = allInd(1:end-modNum);
                    div = repmat(this.frPerBlock, length(fullBlockInd)/this.frPerBlock, 1);
                    fullBlockInd = mat2cell(fullBlockInd, div);
                    this.blocks = [ fullBlockInd; modBlockInd ];
                else
                    this.blocks = modBlockInd;
                end
                this.numBlocks = length(this.blocks);
                
                % Makes cell array 3 columns - 'indices', 'images', 'flags'
                this.blocks = [ ...
                    this.blocks, ...
                    cell(size(this.blocks)), ...
                    num2cell(false(this.numBlocks, 1)) ];
                
                % Calculates the maxium number of blocks
                user = memory;  % check system memory
                maxBytes = min(user.MaxPossibleArrayBytes, maxBytes); % confine user input
                a = ones(1, this.dataType);
                s = whos('a');
                bytePerBlock = prod(this.imgSize(1:2)) * this.frPerBlock * s.bytes;
                this.maxNumBlocks = max(1, floor(maxBytes / bytePerBlock));
                
                % Setup timer
                this.timerObj = timer;
                this.timerObj.TimerFcn = @(~,~) LoadBlock(this);
                this.timerObj.StartDelay = loadingInterval;
                this.timerObj.Period = loadingInterval;
                if this.maxNumBlocks < this.numBlocks
                    this.timerObj.TasksToExecute = inf;
                else
                    this.timerObj.TasksToExecute = this.numBlocks + 1;  % plus 1 to show completion message
                end
                this.timerObj.ExecutionMode = 'fixedSpacing';
                
                % TIFF NOW!
                this.LineUpFrom(1);
                start(this.timerObj);
            end
        end
        
        % User gets image data with conventional indexing (but no 'end' keyword)
        function img = img(this, dim1, dim2, dim3)
            % Records performance in verbose mode
            if this.verbose
                tic;
            end
            
            % Translates indexing
            if strcmp(dim1, ':')
                dim1 = 1 : this.imgSize(1);
            end
            if strcmp(dim2, ':')
                dim2 = 1 : this.imgSize(2);
            end
            if nargin < 4
                dim3 = ':';
            end
            if strcmp(dim3, ':')
                dim3 = 1 : this.imgSize(3);
            end
            
            % Loads image frame by frame
            img = zeros(length(dim1), length(dim2), length(dim3), this.dataType);
            for i = 1 : length(dim3)
                img(:,:,i) = this.GetFrame(dim3(i), dim1, dim2);
            end
            
            % Reports performance in verbose mode
            if this.verbose
                fprintf('image retreived: ');
                toc;
            end
        end
        
        % Kills the timer object
        function delete(this)
            if ~isempty(this.timerObj)
                stop(this.timerObj);
                delete(this.timerObj);
            end
            if ~isempty(this.tifObj)
                close(this.tifObj);
            end
        end
    end
    
    methods(Access = private)
        % Sets the queue from the block of access
        function LineUpFrom(this, accessIdx)
            blockInd = 1 : this.numBlocks;                          % indices of blocks
            flags = cell2mat(this.blocks(:,3));                     % logical mask of loading status
            unloadedInd = blockInd(~flags);                         % keeps unloaded block indices
            dists = abs(unloadedInd - accessIdx);                   % finds distances to the block of interst
            [ ~, indOrder ] = sort(dists);                          % sorts distances in ascending order
            maxLength = min(length(indOrder), this.maxNumBlocks);   % imposes the memory limit
            this.queue = unloadedInd(indOrder(1:maxLength));        % queues unloaded indices with sorted order
        end
        
        % Checks and manages memory
        function ManageMemory(this, accessIdx)
            blockInd = 1 : this.numBlocks;                          % indices of blocks
            flags = cell2mat(this.blocks(:,3));                     % logical mask of loading status
            loadedInd = blockInd(flags);                            % keeps loaded block indices
            dists = abs(loadedInd - accessIdx);                     % finds distances to the block of interst
            [ ~, indOrder ] = sort(dists, 'descend');               % sorts distances in ascending order
            antiQueue = loadedInd(indOrder);                        % queues loaded indices with sorted order
            
            numBlockRemove = length(antiQueue) + 1 - this.maxNumBlocks;     % plus 1 for loading the next block
            if numBlockRemove > 0
                for i = antiQueue(1:min(length(antiQueue),numBlockRemove))
                    this.blocks{i,2} = [];
                    this.blocks{i,3} = false;
                end
            end
        end
        
        % Load next block in the queue
        function LoadBlock(this)
            if ~isempty(this.queue)
                % Gets the index of next block
                blockIdx = this.queue(1);
                if this.verbose
                    tic;
                end
                
                % Memory check and management
                this.ManageMemory(blockIdx);
                
                % Preallocates space for image block
                bSize = length(this.blocks{blockIdx,1});
                this.blocks{blockIdx,2} = zeros([ this.imgSize(1:2) bSize ], this.dataType);
                
                % Reads from file
                for i = 1 : bSize
                    frIdx = this.blocks{blockIdx,1}(i);
                    this.tifObj.setDirectory(frIdx);
                    this.blocks{blockIdx,2}(:,:,i) = this.tifObj.read();
                end
                
                % Updates load/unload flag and chronicle
                this.blocks{blockIdx,3} = true;
                this.queue(1) = [];
                
                if this.verbose
                    if this.numBlocks > 150
                        fprintf('block %d loaded - ', blockIdx);
                    else
                        flags = cell2mat(this.blocks(:,3))';
                        fprintf(char([ flags + '0', 1 ]));
                    end
                    toc;
                    if isempty(this.queue)
                        disp('end of the queue');
                    end
                end
            end
        end
        
        % Retrieves a single frame and (re)directs the reading process (with optional size)
        function fr = GetFrame(this, frIdx, dim1, dim2)
            % Translates the indexing
            blockIdx = ceil(frIdx / this.frPerBlock);
            subIdx = mod(frIdx, this.frPerBlock);
            if subIdx == 0
                subIdx = this.frPerBlock;
            end
            
            if ~this.blocks{blockIdx,3}
                % Directs TimerFnc to the block of interest
                this.LineUpFrom(blockIdx);
                
                % Waits for the block to be available if not already
                ready = false;
                while ~ready
                    pause(0.01);
                    ready = this.blocks{blockIdx,3};
                end
            end
            
            % Returns the frame
            fr = this.blocks{blockIdx,2}(dim1,dim2,subIdx);
        end
    end
    
end

