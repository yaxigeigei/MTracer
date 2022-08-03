classdef Img34 < MImgBaseClass
    %IMAGE234 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods(Static)
        function imgCells = Array2Cell(img)
            %Array2Cell Summary of this function goes here
            %   Detailed explanation goes here
            
            if ~iscell(img)
                imgCells = cell(size(img, 4), 1);
                for i = 1 : size(img, 4)
                    imgCells{i} = img(:,:,:,i);
                end
            else
                imgCells = img;
            end
        end
        
        function [ img, ijRange ] = Crop(img, ijRange, zRange)
            %Crops image stack or stacks with specified range
            %
            % Interactively:    [ img, ijRange ] = Img34.Crop(img)
            % Programatically:  [ img, ijRange ] = Img34.Crop(img, ijRange)
            %                   [ img, ijRange ] = Img34.Crop(img, ijRange, zRange)
            %
            %Input discription
            %
            % img: a 4D array or 3D arrays organized in a 1D cell array
            % ijRange(optional): could be a 2-by-2 array of [ iBegin, iEnd; jBegin, jEnd ], specifying the 
            % desired row and column (instead of x,y) ranges. You can also specify different ranges for each 
            % image stack in the 3rd dimension with a 2-by-2-by-#stacks array.
            % zRange(optional): could be a 1-by-2 array of [ zBegin, zEnd ], specifying the desired range of 
            % depth. You can also specify different ranges for each image stack which requires a
            % #stacks-by-2 array. The default is full range. 
            % 
            %Output discription
            %
            % img: the cropped image organized as a cell array
            % ijRange: the row and column range(s) used for cropping
            
            img = Img34.Array2Cell(img);
            [ imgHeight, imgWidth, imgSlices ] = size(img{1});
            
            if nargin < 3
                zRange = [ 1, imgSlices ];
            end
            
            % Do manual cropping given minimal input
            if nargin < 2
                figure('Name', 'Crop the image if you want (otherwise just close the window)')
                imagesc(Img23.ProjMean(img{1}));
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
            zRange = MMath.Bound(zRange, [ 1 imgSlices ]);
            
            % Expands reduced input
            if size(ijRange,3) == 1
                ijRange = repmat(ijRange, 1, 1, length(img));
            end
            if size(zRange,1) == 1
                zRange = repmat(zRange, length(img), 1);
            end
            
            % Applys cropping
            for i = length(img) : -1 : 1
                imgCut{i} = img{i}( ...
                    ijRange(1,1,i):ijRange(1,2,i), ...
                    ijRange(2,1,i):ijRange(2,2,i), ...
                    zRange(i,1):zRange(i,2));
            end
            img = imgCut;
        end
        
        function [ img, ijRange, zRange ] = CropByShifts(img, centerShifts)
            %Crops image stacks by shifting image centers (e.g. for aligning stacks or
            %tracking an region in volumn imaging result);
            %
            % [ img, ijRange, zRange ] = Img34.CropByShifts(img, centerShifts)
            % 
            %Input discription
            %
            % img: a three dimensional array to be cropped. (If the input is a 2D array, the output will be the
            % same array since it will always overlap with itself no matter what shifting vector you give.)
            % shiftVects: an #stacks-by-3 array where each row specifies [ rowShift, columnShift, depthShift ] of
            % the center of corresponding frame/slice in the image stack.
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
            zMax = max(centerShifts(:,3));
            zMin = min(centerShifts(:,3));
            
            ijRange = zeros(2,2,imgSlices);
            zRange = zeros(length(img),2);
            
            ijRange(1,1,:) = (centerShifts(:,1) - iMin) + 1;
            ijRange(1,2,:) = (centerShifts(:,1) - iMax) + imgHeight;
            ijRange(2,1,:) = (centerShifts(:,2) - jMin) + 1;
            ijRange(2,2,:) = (centerShifts(:,2) - jMax) + imgWidth;
            zRange(:,1)    = (centerShifts(:,3) - zMin) + 1;
            zRange(:,2)    = (centerShifts(:,3) - zMax) + imgSlices;
            
            img = Img34.Crop(img, ijRange, zRange);
        end
        
        function [ stacks, centers, shifts, ijRange ] = RegManual(stacks)
            %RegManual Summary of this function goes here
            %   Detailed explanation goes here
            
            stacks = Img34.Array2Cell(stacks);
            
            figure('Name', 'Click Frame by Frame at a Reference Point')
            
            centers = zeros(length(stacks), 2);
            shifts = zeros(length(stacks), 2);
            
            try
                for i = 1 : length(stacks)
                    imagesc(imadjust(max(stacks{i}, [ ], 3)));
                    axis equal tight
                    colormap gray
                    h = impoint;
                    centers(i, :) = getPosition(h);
                    shifts(i, :) = centers(i,:) - centers(1,:);  % Vector of column shifts and row shifts
                end

                shifts = round([ shifts(:,2), shifts(:,1) ]); % Transform into the order of row shifts and column shifts
                stacks = StacksCrop(stacks, resultSize, shifts);
            catch
                stacks = [];
                centers = [];
                shifts = [];
                ijRange = [];
            end
        end
    end
    
end

