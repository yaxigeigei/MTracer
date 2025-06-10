classdef Img23 < MImgBaseClass
    %IMAGE23 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods(Static)
        function [ img, ijRange ] = Crop(img, ijRange)
            %Crops a two or three dimensional array with specified range
            %
            % Interactively:    [ img, ijRange ] = Img23.Crop(img)
            % Programatically:  [ img, ijRange ] = Img23.Crop(img, ijRange)
            % 
            %Input discription
            %
            % img: a two or three dimensional array to be cropped
            % ijRange(optional): could be an 2-by-2 array of [ iBegin, iEnd; jBegin, jEnd ], specifying the 
            % desired row and column (instead of x,y) ranges. You can also specify different ranges for each 
            % frames/slices in the 3rd dimension and the size of ijRange will be 2-by-2-by-size(img,3). 
            % 
            %Output discription
            %
            % img: the resulting array
            % ijRange: the row and column range(s) used for cropping
            
            if iscell(img)
                img = img{1};
            end
            [ imgHeight, imgWidth, imgSlices ] = size(img);
            
            % Do manual cropping given minimal input
            if nargin < 2
                figure('Name', 'Crop the image if you want (otherwise just close the window)')
                imagesc(Img23.ProjMean(img));
                axis equal tight
                colormap gray
                
                try
                    hCrop = imrect;
                    ijRange = wait(hCrop);                          % [ columnBegin, rowBegin, width, height ]
                    ijRange(3:4) = ijRange(1:2) +ijRange(3:4);      % [ columnBegin, rowBegin, columnEnd, rowEnd ]
                    ijRange = round(ijRange);
                    ijRange = reshape(ijRange, 2, 2);               % [ columnBegin, columnEnd; rowBegin, rowEnd ]
                    ijRange = ijRange([2,1],:);                     % [ rowBegin, rowEnd; columnBegin, columnEnd ]
                catch
                    ijRange = [ 1, imgHeight; 1, imgWidth ];
                end
            end
            
            % Checks (and fixes) boundaries
            ijRange(1,:,:) = MMath.Bound(ijRange(1,:,:), [ 1 imgHeight ]);
            ijRange(2,:,:) = MMath.Bound(ijRange(2,:,:), [ 1 imgWidth ]);
            
            % Expands reduced input
            if size(ijRange,3) == 1
                ijRange = repmat(ijRange, 1, 1, imgSlices);
            end
            
            % Applys cropping
            for i = imgSlices : -1 : 1
                imgCut(:,:,i) = img(ijRange(1,1,i):ijRange(1,2,i), ijRange(2,1,i):ijRange(2,2,i), i);
            end
            img = imgCut;
        end
        
        function [ img, ijRange ] = CropByShifts(img, centerShifts)
            %Crops a three dimensional array by shifting image centers (e.g. for aligning frames/slices or
            %tracking an object in a video);
            %
            % [ img, ijRange ] = Img23.CropByShifts(img, centerShifts)
            % 
            %Input discription
            %
            % img: a three dimensional array to be cropped. (If the input is a 2D array, the output will be the
            % same array since it will always overlap with itself no matter what shifting vector you give.)
            % shiftVects: an size(img,3)-by-2 array where each row specifies [ rowShift, columnShift ] of the 
            % center of corresponding frame/slice in the image stack.
            %
            %Output discription
            %
            % img: the resulting array
            % ijRange: the row and column range(s) used for cropping
            
            if iscell(img)
                img = img{1};
            end
            [ imgHeight, imgWidth, imgSlices ] = size(img);
            
            iMax = max(centerShifts(:,1));
            iMin = min(centerShifts(:,1));
            jMax = max(centerShifts(:,2));
            jMin = min(centerShifts(:,2));
            
            ijRange = zeros(2,2,imgSlices);
            
            ijRange(1,1,:) = 1 +            (centerShifts(:,1) - iMin);
            ijRange(1,2,:) = imgHeight +    (centerShifts(:,1) - iMax);
            ijRange(2,1,:) = 1 +            (centerShifts(:,2) - jMin);
            ijRange(2,2,:) = imgWidth +     (centerShifts(:,2) - jMax);
            
            img = Img23.Crop(img, ijRange);
        end
        
        function img = Filt2(img, h)
            %Performs 2D imfilter() along the 3rd dimension
            %
            % img = Img23.Filt2(img, h)
            % 
            %Input discription
            %
            % img: a two or three dimensional array to be filtered
            % h: a fspecial() filter object
            % 
            %Output discription
            %
            % img: the resulting image array
            
            for i = 1 : size(img, 3)
                img(:,:,i) = imfilter(img(:,:,i), h);
            end
        end
        
        function img = Filt2Average(img, radius)
            %Performs 2D averaging along the 3rd dimension
            %
            % img = Img23.Filt2Average(img, radius)
            % 
            %Input discription
            %
            % img: a two or three dimensional array to be filtered
            % radius: the radius of the spatial filter in pixels
            % 
            %Output discription
            %
            % img: the resulting image array
            
            if nargin < 2
                radius = 3;
            end
            
            h = fspecial('average', radius);
            
            for i = 1 : size(img, 3)
                img(:,:,i) = imfilter(img(:,:,i), h);
            end
        end
        
        function img = Filt2Gaussian(img, radius, sigma)
            %Performs 2D averaging along the 3rd dimension
            %
            % img = Img23.Filt2Average(img, radius, sigma)
            % 
            %Input discription
            %
            % img: a two or three dimensional array to be filtered
            % radius: the radius of the spatial filter in pixels
            % sigma: the spread of the spatial filter in pixels
            % 
            %Output discription
            %
            % img: the resulting image array
            
            if nargin < 2
                radius = 3;
                sigma = 0.5;
            end
            
            h = fspecial('gaussian', radius, sigma);
            
            for i = 1 : size(img, 3)
                img(:,:,i) = imfilter(img(:,:,i), h);
            end
        end
        
        function img = Filt2Median(img, filterSize)
            %Performs 2D averaging along the 3rd dimension
            %
            % img = Img23.Filt2Average(img, filterSize)
            % 
            %Input discription
            %
            % img: a two or three dimensional array to be filtered
            % filterSize: the length of the edge of the spatial filter in pixels
            % 
            %Output discription
            %
            % img: the resulting image array
            
            if nargin < 2 
                filterSize = [3, 3];
            end
            
            for i = 1 : size(img, 3)
                img(:,:,i) = medfilt2(img(:,:,i), filterSize);
            end
        end
        
        function img = Filt3Average(img, radius)
            %Filt3Average Summary of this function goes here
            %   Detailed explanation goes here
            
            if iscell(img)
                img = img{1};
            end
            
            imgFrames = size(img, 3);
            
            imgAvg = zeros(size(img), 'like', img);
            for i = 1 : imgFrames
                frameBound(1) = i - radius;
                frameBound(2) = i + radius;
                frameBound = MMath.Bound(frameBound, [ 1 imgFrames ]);
                imgAvg(:,:,i) = mean(img(:,:,frameBound(1):frameBound(2)), 3);
            end
        end
        
        function Play(img, interval)
            %Play Summary of this function goes here
            %   Detailed explanation goes here
            
            if nargin < 2
                interval = max(10/size(img,3), 0.03);
            end
            
            figure
            h = imagesc(img(:,:,1));
            colormap gray
            axis equal tight
            caxis([min(img(:)), max(img(:))]);
            
            try
                for i = 2 : size(img,3)
                    set(h, 'CData', img(:,:,i));
                    pause(interval);
                end
            catch
            end
        end
        
        function Record(path, img, varargin)
            %Play Summary of this function goes here
            %   Detailed explanation goes here
            
            [ imgHeight, imgWidth, imgSlices ] = size(img);
            
            p = inputParser();
            p.addParameter('frameRate', max(imgSlices/10,2), @isnumeric);
            p.parse(varargin{:});
            frameRate = p.Results.frameRate;
            
            hFig = figure('Name', 'You can make adjustment to the window then click on the image to start');
            figPos = get(hFig, 'Position');
            set(hFig, 'Position', [ figPos(1:2), imgWidth imgHeight ]);
            axes('Parent', hFig, 'Units', 'pixels', 'Position', [ 0 1 imgWidth imgHeight ]);
            
            hold on
            ih = imagesc(img(:,:,1));
            xlim([ 1 size(img,2) ]);
            ylim([ 1 size(img,1) ]);
            axis ij equal off
            colormap gray
            
            try
                ginput(1);
                
                vidObj = VideoWriter(path);
                vidObj.FrameRate = frameRate;
                open(vidObj);
                
                for i = 1 : imgSlices
                    set(ih, 'CData', img(:,:,i));
                    drawnow;
                    frameObj = getframe;
                    writeVideo(vidObj, frameObj.cdata);
                end

                close(vidObj);
            catch
            end
        end
        
        function [ stack, centers, shifts, ijRange ] = RegManual(stack)
            %RegManual Summary of this function goes here
            %   Detailed explanation goes here
            
            figure('Name', 'Click Frame by Frame at a Reference Point')
            
            centers = zeros(size(stack,3), 2);
            shifts = zeros(size(stack,3), 2);
            
            try
                for i = 1 : size(stack, 3)
                    imagesc(imadjust(stack(:,:,i)));
                    axis equal tight
                    colormap gray
                    h = impoint;
                    centers(i, :) = getPosition(h);
                    shifts(i, :) = centers(i,:) - centers(1,:);
                end

                shifts = round(shifts);

                [ stack, ijRange ] = Img23.CropByShifts(stack, shifts);
            catch
                stack = [];
                centers = [];
                shifts = [];
                ijRange = [];
            end
        end
        
        function img = Scale(img, xFactor, yFactor, zFactor)
            if nargin < 4
                zFactor = 1;
            end
            
            xInterval = 1 / xFactor;
            yInterval = 1 / yFactor;
            zInterval = 1 / zFactor;
            
            dataType = class(img);
            img = single(img);
            
            if ndims(img) == 3
                [ XI, YI, ZI ] = meshgrid(1:xInterval:size(img, 2), 1:yInterval:size(img, 1), 1:zInterval:size(img, 3));
                XI = single(XI);
                YI = single(YI);
                ZI = single(ZI);
                img = interp3(img, XI, YI, ZI);
            else
                [ X, Y ] = meshgrid(1:size(img, 2), 1:size(img, 1));
                [ XI, YI ] = meshgrid(1:xInterval:size(img, 2), 1:yInterval:size(img, 1));
                img = interp2(X, Y, img, XI, YI);
            end
            
            img = cast(img, dataType);
        end
    end
    
end

