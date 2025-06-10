classdef MNN
    %MNN Summary of this class goes here
    %   Detailed explanation goes here
    
    methods(Static)
        function [vidTf, tform] = RoiTransform(vid, roiTemplate, tfStr)
            % Transform (crop) a video to a region of interest that matches the template image
            % 
            %   [vidTf, tform] = MNN.RoiTransform(vid, roiTemplate, tfStr)
            % 
            % 
            
            % Handle user inputs
            if nargin < 3
                tfStr = 'XY';
            end
            tfStr = upper(tfStr);
            
            % Compute an average frame
            if ndims(vid) < 4
                frMean = mean(vid, 3);
            else
                vidMono = squeeze(mean(vid, 3)); % convert RGB to mono
                frMean = mean(vidMono, 3);
            end
            frMean = cast(frMean, 'like', roiTemplate); % make sure the data types are the same
            
            % Find transformation by image registration to the template
            [optimizer, metric] = imregconfig('monomodal');
            tform = imregtform(frMean, roiTemplate, 'translation', optimizer, metric);
            
            % Cancel translation in certain axis
            tfCancel = setxor(tfStr, 'XY');
            for i = 1 : numel(tfCancel)
                switch tfCancel(i)
                    case 'X'
                        tform.T(3,1) = 0;
                    case 'Y'
                        tform.T(3,2) = 0;
                    otherwise
                        warning('Canceling transformation in %s is not defined or not supported', tfCancel(i));
                end
            end
            
            % Crop out ROI
            vidTf = imwarp(vid, tform, 'OutputView', imref2d(size(roiTemplate)));
        end
        
        function [vid, t, vObj] = ReadVideo(vidFilePaths, varargin)
            % Read video file(s)
            % 
            %   [vid, t, vObj] = MNN.ReadVideo(vidFilePaths)
            %   [vid, t, vObj] = MNN.ReadVideo(..., 'FrameFunc', @(x) x)
            % 
            % Inputs
            %   vidFilePaths        A video file path or a cell array of paths for multiple videos. 
            %   'FrameFunc'         The handle of a function that processes each frame upon reading. For example, 
            %                       you can directly convert RGB videos to Mono with @rgb2gray. The default does
            %                       nothing. 
            % Outputs
            %   vid                 An array of size [height,width,frames] or [height,width,colors,frames], or a 
            %                       cell array of them for multiple videos. 
            %   t                   A [frames,1] vector of frame time in second, or a cell array of them for 
            %                       multiple videos.
            %   vObj                A VideoReader object, or a cell array of them for multiple videos. 
            % 
            
            % Handle user inputs
            p = inputParser();
            p.addParameter('FrameFunc', @(x) x);
            p.parse(varargin{:});
            frameFunc = p.Results.FrameFunc;
            
            if nargin < 1 || isempty(vidFilePaths)
                vidFilePaths = MBrowse.Files([], 'Please select video file(s) to read');
            end
            vidFilePaths = cellstr(vidFilePaths);
            
            % Read each video file
            vid = cell(numel(vidFilePaths), 1);
            t = cell(numel(vidFilePaths), 1);
            vObj = cell(numel(vidFilePaths), 1);
            
            for i = 1 : numel(vid)
                % Construct a VideoReader
                vObj{i} = VideoReader(vidFilePaths{i});
                
                % Read all frames
                s = struct('cdata', [], 'time', 0);
                vObj{i}.CurrentTime = 0;
                k = 0;
                while hasFrame(vObj{i})
                    k = k + 1;
                    s(k).time = vObj{i}.CurrentTime;
                    s(k).cdata = frameFunc(readFrame(vObj{i}));
                end
                vid{i} = cat(ndims(s(k).cdata)+1, s.cdata);
                t{i} = cat(1, s.time);
            end
            
            % Denest if there is only one video
            if numel(vid) == 1
                vid = vid{1};
                t = t{1};
                vObj = vObj{1};
            end
        end
        
        function img = ReadImage2RGB(filePath)
            % Read an image file and output the image with RGB layers
            
            img = imread(filePath);
            
            switch size(img,3)
                case 1
                    img = repmat(img, [1 1 3]);
                case 3
                    % do nothing
                otherwise
                    error('Cannot convert image to RGB with %d layers', size(img,3));
            end
        end
        
        function tbOut = DenestTable(tbIn)
            % Denest table columns where each element is a cell array
            
            varNames = tbIn.Properties.VariableNames;
            tbOut = table();
            for i = 1 : width(tbIn)
                tbOut.(varNames{i}) = cat(1, tbIn.(varNames{i}){:});
            end
        end
        
        function [imgTf, ptsTf, randParams] = AugmentImage(augObj, img, pts)
            % Manually augment images
            
            edgeLen = size(img);
            edgeLen = edgeLen(1:2);
            halfEdge = (edgeLen-1)/2;
            
            a = getRandFrom(augObj.RandRotation);
            sx = getRandFrom(augObj.RandXScale);
            sy = getRandFrom(augObj.RandYScale);
            tx = getRandFrom(augObj.RandXTranslation);
            ty = getRandFrom(augObj.RandYTranslation);
            
            if augObj.RandXReflection && rand(1) > 0.5
                sx = -sx;
            end
            if augObj.RandYReflection && rand(1) > 0.5
                sy = -sy;
            end
            
            centerMat = [1 0 0; 0 1 0; -halfEdge 1];
            rotateMat = [cosd(a) sind(a) 0; -sind(a) cosd(a) 0; 0 0 1];
            scaleMat = [sx 0 0; 0 sy 0; 0 0 1];
            transMat = [1 0 0; 0 1 0; tx ty 1];
            resetMat = [1 0 0; 0 1 0; halfEdge 1];
            affObj = affine2d(centerMat*rotateMat*scaleMat*transMat*resetMat);
            
            refObj = imref2d(edgeLen);
            imgTf = imwarp(img, affObj, 'OutputView', refObj);
            
            ptsTf = transformPointsForward(affObj, pts);
            ptsTf(isnan(ptsTf)) = 0;
            
            randParams = [a sx sy tx ty];
            
            function randVal = getRandFrom(valRange)
                randVal = valRange(1) + rand(1) * diff(valRange);
            end
        end
        
    end
end

