classdef MImgBaseClass
    %MImgBaseClass Summary of this class goes here
    %   Detailed explanation goes here
    
    methods(Static)
        function Export(path, img, indexFormat)
            % Export image array as a TIFF file (will overwrite the file with the same name)
            %
            %   ManyFuncImgBaseClass.Export(path, img)
            %   ManyFuncImgBaseClass.Export(path, img, indexFormat)
            %
            % Inputs
            %   path            The path (directory + file name without suffix) to save
            %   img             A two, three or four dimensional array, or a cell array of images(or stacks)
            %   indexFormat     A string for formated indexing of individual stacks
            %                   (default is '_%04d', i.e. an underscore followed by four digits integer)
            %
            
            if nargin < 3
                indexFormat = '_%04d';
            end
            
            img = Img34.Array2Cell(img);
            
            for i = 1 : length(img)
                a = img{i}(1);
                s = whos('a');
                bps = s.bytes * 8;
                
                if length(img) > 1
                    counting = sprintf(indexFormat, i);
                    tObj = Tiff([ path, counting '.tif' ], 'w');
                else
                    [ ~, ~, ext ] = fileparts(path);
                    if ~strcmp(ext, '.tif')
                        path = [ path, '.tif' ];
                    end
                    tObj = Tiff(path, 'w');
                end
                
                tagstruct.ImageLength = size(img{i}, 1);
                tagstruct.ImageWidth = size(img{i}, 2);
                tagstruct.BitsPerSample = bps;
                tagstruct.Compression = Tiff.Compression.None;
                tagstruct.SamplesPerPixel = 1;
                tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
                tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
                
                tObj.setTag(tagstruct);
                tObj.write(img{i}(:,:,1));
                
                for j = 2 : size(img{i}, 3)
                    tObj.writeDirectory();
                    tObj.setTag(tagstruct);
                    tObj.write(img{i}(:,:,j))
                end
                
                tObj.close();
            end
        end
        
        function [ imgHeight, imgWidth, imgSlices ] = GetTiffImageSize( tifObj )
            %GETTIFFOBJFRAMENUM Summary of this function goes here
            %   Detailed explanation goes here
            
            % Convert input to a Tiff object
            if ischar(tifObj)
                tifObj = Tiff(tifObj, 'r');
            end
            
            % Get the height and width
            imgWidth = tifObj.getTag('ImageWidth');
            imgHeight = tifObj.getTag('ImageLength');
            
            % Determine the search range
            probeIncrement = 1000;
            upperLim = 1;
            
            while true
                try
                    tifObj.setDirectory(upperLim);
                    lowerLim = upperLim;
                    upperLim = upperLim + probeIncrement;
                catch
                    break;
                end
            end
            
            % Dichotomizing the search range
            while upperLim - lowerLim > 1
                midDir = floor((upperLim + lowerLim)/2);
                try
                    tifObj.setDirectory(midDir);
                    lowerLim = midDir;
                catch
                    upperLim = midDir;
                end
            end
            
            try
                tifObj.setDirectory(upperLim);
                imgSlices = upperLim;
            catch
                imgSlices = lowerLim;
            end
        end

        function img = Import(path, varargin)
            % Imports a TIFF file (conventional non-dynamic loading)
            %
            %   img = ManyFuncImgBaseClass.Import(path)
            %   img = ManyFuncImgBaseClass.Import(path, ..., 'range', [1 Inf])
            %   img = ManyFuncImgBaseClass.Import(path, ..., 'interleave', [1 1])
            %   img = ManyFuncImgBaseClass.Import(path, ..., 'dataType', '')
            %
            % Inputs
            %    path           The path of the TIFF file. File browsing UI will prompt if passing empty value.
            %   'range'         Two-element vector indicating the begining and ending index of the range of interest.
            %                   (default is the entire tif)
            %   'interleave'    Two-element vector. The first specifies total number of channels; the second specifies
            %                   the channel you want to keep. (default [1,1])
            %   'dataType'      The data type of the image. (default 'uint8' for 8-bit TIFF, 'uint16' for 16-bit,
            %                   'double' for other depths)
            %
            % Output
            %   img             A two or three dimensional array
            
            p = inputParser;
            p.addParameter('range', [ 1 inf ], @isnumeric);
            p.addParameter('interleave', [ 1 1 ], @isnumeric);
            p.addParameter('dataType', '', @ischar);
            p.parse(varargin{:});
            frameBegin = p.Results.range(1);
            frameEnd = p.Results.range(2);
            channelTotal = p.Results.interleave(1);
            channelSelect = p.Results.interleave(2);
            datatype = p.Results.dataType;
            
            if nargin < 1 || isempty(path)
                path = MBrowse.File();
            end
            
            warning('off', 'MATLAB:imagesci:Tiff:libraryWarning');
            
            [ ~, ~, ext ] = fileparts(path);
            if ~strcmp(ext, '.tif')
                path = [ path, '.tif' ];
            end
            tObj = Tiff(path, 'r');
            imgWidth = tObj.getTag('ImageWidth');
            imgHeight = tObj.getTag('ImageLength');
            
            if isempty(datatype)
                pxlBytes = round(tObj.getTag('BitsPerSample') / 8);
                switch pxlBytes
                    case 1
                        datatype = 'uint8';
                    case 2
                        datatype = 'uint16';
                    otherwise
                        datatype = 'double';
                end
            end
            
            frameObj = 1;
            while ~tObj.lastDirectory()
                frameObj = frameObj + 1;
                tObj.nextDirectory();
            end
            
            frameEnd = floor(min([ frameObj / channelTotal, frameEnd ]));
            
            img = zeros(imgHeight, imgWidth, frameEnd-frameBegin+1, datatype);
            for i = frameBegin : frameEnd
                tObj.setDirectory((i-1) * channelTotal + channelSelect);
                img(:,:,i-frameBegin+1) = tObj.read();
            end
            
            tObj.close();
        end
        
        function stack = Import2Stack(dirPath, key)
            % Combine multiple TIFF files into a stack or stacks
            %
            %   stack = ManyFuncImgBaseClass.Import2Stack()
            %   stack = ManyFuncImgBaseClass.Import2Stack(dirPath)
            %   stack = ManyFuncImgBaseClass.Import2Stack(dirPath, key)
            %
            % Inputs
            %   dirPath     The directory of TIFF files. By default prompting browser to select files
            %   key         A string of keyword with wildcards to specify a group of files (e.g. 'ChanA_????_*.tif'). Files
            %               within a group must have the same size in every dimensions. You can also use formated number
            %               (e.g. '*_%04d.tif') to incrementally specify different groups of files and import them into
            %               respective stacks. Images across stacks are not required to have the same size. (default is '*.tif')
            % Output
            %   stacks      A three dimensional array or a cell array containing multiple 3D arrays, depending on the
            %               grouping specified by the keyword
            %
            
            if nargin < 1
                % Uses UI to select files
                [ filePaths, dirPath ] = MBrowse.Files();
                
                % Groups them into one stack
                list{1} = filePaths;
            else
                % Default keyword finds all TIFF files in the folder
                if nargin < 2
                    key = '*.tif';
                end
                
                % Initializes grouping
                list = cell(1);
                listEnd = false;
                i = 1;
                
                while(~listEnd)
                    % Applies indexing, if any, to the keyword
                    fileName = sprintf(key, i);
                    
                    % Finds all files satisfying the keyword
                    tempFilePaths = dir(fullfile(dirPath, fileName));
                    tempFilePaths = struct2cell(tempFilePaths);
                    tempFilePaths = tempFilePaths(1,:)';
                    
                    % Saves above file paths for current image stack
                    if ~isempty(tempFilePaths)
                        list{i} = tempFilePaths;
                    else
                        % Stops grouping if nothing more to do
                        listEnd = true;
                    end
                    
                    % Increments index for next stack
                    i = i + 1;
                    
                    % Stops grouping if the keyword is not incremental
                    if strcmp(fileName, sprintf(key, i))
                        listEnd = true;
                    end
                end
            end
            
            % Loading images to stack(s)
            stack = cell(length(list), 1);
            if ~isempty(list{1})
                for i = 1 : length(list)
                    % Gets the dimension info for preallocation
                    firstImg = ManyFuncImgBaseClass.Import(fullfile(dirPath, list{i}{1}));
                    [ imgHeight, imgWidth, imgSlices ] = size(firstImg);
                    stack{i} = zeros([ imgHeight, imgWidth, imgSlices * length(list{i}) ], 'like', firstImg);
                    
                    % Loads files into the current stack
                    stack{i}(:,:,1:imgSlices) = firstImg;
                    for j = 2 : length(list{i})
                        stack{i}(:, :, (j-1)*imgSlices+1 : j*imgSlices) = ManyFuncImgBaseClass.Import(fullfile(dirPath, list{i}{j}));
                    end
                end
            end
            
            if length(stack) == 1
                stack = stack{1};
            end
        end
        
        function proj = ProjCustom(img, hFnc)
            %ProjCustom Summary of this function goes here
            %   Detailed explanation goes here
            
            img = Img34.Array2Cell(img);
            
            proj = zeros(size(img{1},1), size(img{1},2), length(img), 'like', img{1});
            for i = 1 : length(img)
                proj(:,:,i) = hFnc(img{i});
            end
        end
        
        function proj = ProjMax(img)
            %ProjMax Summary of this function goes here
            %   Detailed explanation goes here
            
            img = Img34.Array2Cell(img);
            
            proj = zeros(size(img{1},1), size(img{1},2), length(img), 'like', img{1});
            for i = 1 : length(img)
                proj(:,:,i) = max(img{i}, [ ], 3);
            end
        end
        
        function proj = ProjMean(img)
            %ProjAverage Summary of this function goes here
            %   Detailed explanation goes here
            
            img = Img34.Array2Cell(img);
            
            proj = zeros(size(img{1},1), size(img{1},2), length(img), 'like', img{1});
            for i = 1 : length(img)
                proj(:,:,i) = mean(img{i}, 3);
            end
        end
        
        function [ comboStack, xProj, yProj, zProj ] = ProjXYZ(img, zFactor)
            %ProjCombo Summary of this function goes here
            %   Detailed explanation goes here
            
            if nargin < 2
                zFactor = 1;
            end
            
            img = Img34.Array2Cell(img);
            numStacks = length(img);
            [ imgHeight, imgWidth, imgSlices ] = size(img{1});
            
            xProj = zeros(imgHeight, imgSlices, numStacks, 'like', img{1});
            yProj = zeros(imgSlices, imgWidth, numStacks, 'like', img{1});
            zProj = zeros(imgHeight, imgWidth, numStacks, 'like', img{1});
            
            for i = 1 : numStacks
                xProj(:,:,i) = squeeze(max(img{i}, [ ], 2));
                yProj(:,:,i) = squeeze(max(img{i}, [ ], 1))';
                zProj(:,:,i) = max(img{i}, [ ], 3);
            end
            
            xProj = Img23.Scale(xProj, zFactor, 1, 1);
            yProj = Img23.Scale(yProj, 1, zFactor, 1);
            
            imgSlices = size(yProj,1);
            comboStack = [ yProj, ones(imgSlices, imgSlices, numStacks); zProj, xProj ];
        end
        
        function Viewer(img, varargin)
            % Display image, image stack or stacks in a figure-based viewer for easy browsing
            %
            %   MImgBaseClass.Viewer(img)
            %   MImgBaseClass.Viewer(..., 'delay', 0.03)
            %   MImgBaseClass.Viewer(..., 'userFunc', {@(s,f,d) [], []})
            %
            % Inputs
            %   img                     A two, three or four dimensional array, or a cell array of images (or stacks).
            %   'delay'                 Amount of time to pause after each image update. It allows unfinished plotting to
            %                           finish and affects the maximum speed of browsing. 
            %   'userFunc'              You can make automatic custom plots after each update by providing your function 
            %                           handle and data in a cell array. 
            %                           The first element of this cell array is a handle of a function that takes at least 
            %                           two input arguments (stack index, frame index, ...) and returns handle(s) of 
            %                           plotted objects (otherwise old plots will not be cleared). 
            %                           The later elements in this cell array are input variable(s) expected by the user 
            %                           function in addition to stack index and frame index. 
            % Controls
            %   Mouse
            %     scroll wheel          Go back or forth within the current stack
            %   Keyboard
            %     left/right arrow      Go back or forth within the current stack. Hold Shift key at the same time 
            %                           to jump every 5 frames, Ctrl key every 25 frames. 
            %     up/down arrow         Go to previous or next stack. 
            %     number 0 to 9         Go to 0% to 90% position in the current stack, respectively.
            %     
            
            % Handle user input
            if ~iscell(img)
                img = {img};
            end
            
            if ~all(cellfun(@(x) isnumeric(x) || islogical(x), img))
                error('All image data needs to be numeric or logical array.');
            end
            
            p = inputParser();
            p.addParameter('delay', 0.03, @isscalar);
            p.addParameter('userFunc', {@(s,f,d) [], []}, @iscell);
            p.parse(varargin{:});
            delayTime = p.Results.delay;
            userFuncHandle = p.Results.userFunc{1};
            userFuncArg = p.Results.userFunc(2:end);
            
            % Instruction
            eval('help MImgBaseClass.Viewer');
            
            % Construct UI
            clf;
            fh = gcf;
            set(fh, ...
                'WindowScrollWheelFcn', @mouseScroll, ...
                'WindowKeyPressFcn', @keyPress, ...
                'WindowKeyReleaseFcn', @keyRelease);
            
            hAx = axes();
            hImg = imagesc(); hold on
            axis ij equal tight
            colormap gray
            
            % Initialize display
            sIdx = 1;
            fIdx = 1;
            userPlotHandles = [];
            updateDisplay();
            caxis manual
            
            % Utilities
            function updateDisplay()
                % Image
                sIdx = MMath.Bound(sIdx, [1, length(img)]);
                
                if ndims(img{sIdx}) < 4
                    fIdx = MMath.Bound(fIdx, [1, size(img{sIdx},3)]);
                    hImg.CData = img{sIdx}(:,:,fIdx);
                else
                    fIdx = MMath.Bound(fIdx, [1, size(img{sIdx},4)]);
                    hImg.CData = img{sIdx}(:,:,:,fIdx);
                end
                
                fh.Name = sprintf('stack #%i, frame #%i', sIdx, fIdx);
                
                % User function
                try
                    delete(userPlotHandles);
                    userPlotHandles = userFuncHandle(sIdx, fIdx, userFuncArg{:});
                catch
                    warning('Error occured when processing user function');
                end
                
                pause(delayTime);
            end
            
            % Callbacks
            function keyPress(src, eventdata)
                isUpdate = true;
                
                if any(strcmp(eventdata.Modifier, 'control'))
                    frChange = 25;
                elseif any(strcmp(eventdata.Modifier, 'shift'))
                    frChange = 5;
                else
                    frChange = 1;
                end
                
                switch eventdata.Key
                    case 'rightarrow'
                        fIdx = fIdx + frChange;
                    case 'leftarrow'
                        fIdx = fIdx - frChange;
                    case 'uparrow'
                        sIdx = sIdx - 1;
                    case 'downarrow'
                        sIdx = sIdx + 1;
                    otherwise
                        isUpdate = false;
                end
                
                if isUpdate
                    updateDisplay();
                end
            end
            
            function mouseScroll(src, eventdata)
                signChange = sign(eventdata.VerticalScrollCount);
                frChange = min(abs(eventdata.VerticalScrollCount), 3);
                fIdx = fIdx + signChange*frChange;
                updateDisplay();
            end
            
            function keyRelease(src, eventdata)
                val = str2double(eventdata.Key(end));
                if ~isnan(val)
                    if ndims(img{sIdx}) < 4
                        fIdx = floor(size(img{sIdx},3) * val / 10);
                    else
                        fIdx = floor(size(img{sIdx},4) * val / 10);
                    end
                    updateDisplay();
                end
            end
            
        end
    end
    
end




