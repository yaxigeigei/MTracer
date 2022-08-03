classdef MPlot
    %MPlot A collection of functions useful for plotting
    
    properties(Constant)
        isSave = true;
    end
    
    methods(Static)
        function varargout = Axes(varargin)
            % Change Axes object with custom defaults
            switch numel(varargin)
                case 0
                    ax = axes();
                case 1
                    h = varargin{1};
                    if isa(h, 'matlab.graphics.axis.Axes')
                        ax = h;
                    elseif isa(h, 'matlab.ui.Figure')
                        ax = axes(h);
                    else
                        error('The input must be an axes or a figure handle.');
                    end
                case 3
                    ax = subplot(varargin{:});
            end
            ax.Box = 'off';
            ax.TickDir = 'out';
            if nargout > 0
                varargout{1} = ax;
            end
        end
        
        function varargout = Blocks(xRange, yRange, varargin)
            % Plot rectangles based on X anf Y ranges
            % 
            %   MPlot.Blocks(xRange, yRange)
            %   MPlot.Blocks(xRange, yRange, color)
            %   MPlot.Blocks(..., 'PatchPropertyName', value)
            %   p = MPlot.Blocks(...)
            %
            % Inputs
            %   xRange      Boundary coordinates of block(s) along X-axis in a n-by-2 matrix where n is 
            %               the number of block(s) to plot. Alternatively, you can provide logical vector 
            %               for X-axis where region(s) containing block(s) are set to true. If n is 1, 
            %               xRange will be expanded to match rows in yRange. 
            %   yRange      Same as xRange but for Y-axis. 
            %   color       A 1-by-3 matrix of RGB color for uniform coloring or a 3-by-3 matrix of row 
            %               vectors for a gradient of colors. This parameter applies to all blocks. The 
            %               default color is uniform gray ([.9 .9 .9]).
            %   Any PropertyName-Value pair of MATLAB patch function except for 'EdgeColor'.
            %   
            % Output
            %   p           A Patch object, which can consist of one or more polygons. Use p to query or 
            %               change properties of the patch object after it is created.
            % 
            % See also patch
            
            % Handles user inputs
            if isempty(varargin) || ischar(varargin{1})
                % Default color is gray
                varargin = [{[.9 .9 .9]}, varargin];
            end
            
            % Convert logical mask to boundaries
            if islogical(xRange)
                xRange = MMath.Logical2Bounds(xRange);
            end
            if islogical(yRange)
                yRange = MMath.Logical2Bounds(yRange);
            end
            
            % Apply yRange to all xRanges if necessary
            if numel(yRange) == 2
                yRange = repmat(yRange, size(xRange,1), 1);
            end
            if numel(xRange) == 2
                xRange = repmat(xRange, size(yRange,1), 1);
            end
            
            % Plots blocks
            xx = xRange(:,[1 2 2 1])';
            yy = yRange(:,[1 1 2 2])';
            p = patch(xx, yy, varargin{:}, 'EdgeColor', 'none');
            if nargout > 0
                varargout{1} = p;
            end
        end
        
        function h = Circle(x, y, r, c)
            % Plot a circle
            % 
            %   h = MPlot.Circle(x, y, r, c)
            %
            % Inputs:
            %   x       X-coordinate of the center
            %   y       Y-coordinate of the center
            %   r       Radius of the circle
            %   c       Color of the circle
            % Output:
            %   h       Object handle of the circle shape. By nature, it is a rectangle with rounded corners. 
            
            d = r*2;
            px = x-r;
            py = y-r;
            h = rectangle('Position', [px py d d], 'Curvature', [1 1], 'FaceColor', c, 'LineStyle', 'none');
            daspect([1 1 1]);
        end
        
        function cc = Color2Str(cc, txt)
            % Convert color array to strings in the form '%f,%f,%f' or '\color[rgb]{%f,%f,%f}%s'
            % 
            %   cc = Color2Str(cc)
            %   cc = Color2Str(cc, txt)
            %
            % Inputs
            %   cc      n-by-3 RGB or n-by-4 RGBA color array.
            %   txt     n strings to format.
            % Output
            %   cc      n formated strings. If txt is not provided, each output is formated as 
            %           '%f,%f,%f' or '%f,%f,%f,%f', otherwise as '\color[rgb]{%f,%f,%f}%s' 
            %           where %s is a string in txt.
            cc = mat2str(cc);
            cc = strrep(cc(2:end-1), ' ', ',');
            cc = strsplit(cc, ';')';
            if nargin > 1
                txt = cellstr(txt);
                cc = cellfun(@(c,x) ['\color[rgb]{' c '}' x], ...
                    cc, txt, 'Uni', false);
            end
        end
        
        function ErrorShade(varargin)
            % Plot error as shading
            % 
            %   MPlot.ErrorShade(y, err)
            %   MPlot.ErrorShade(x, y, err)
            %   MPlot.ErrorShade(x, y, errPos, errNeg)
            %   MPlot.ErrorShade(..., 'IsRelative', true)
            %   MPlot.ErrorShade(..., 'Orientation', 'vertical')
            %   MPlot.ErrorShade(..., 'Color', 'k')
            %   MPlot.ErrorShade(..., 'Alpha', 0.15)
            %   MPlot.ErrorShade(..., 'Parent', gca)
            %
            % Inputs:
            %   x               X-coordinates. Default is indices of elements in y. 
            %   y               Y-coordinates. 
            %   err             Errors, which apply to both sides. 
            %   errPos          Errors on the positive side. 
            %   errNeg          Errors on the negative side.
            %   'IsRelative'    Logical variable indicate whether (default) or not error inputs are relative to y. 
            %   'Orientation'   Orientation along which err is applied. 'vertical' (default) for Y-axis 
            %                   and 'horizontal' for X-axis.
            %   'Color'         Color of the shade. Default is black. 
            %   'Alpha'         Transparancy of the shade. Default 0.3. 
            %   'Parent'        The axes to plot in. Default is the current axes returned by gca.
            %
            
            % Handles user inputs
            p = inputParser();
            p.addRequired('arg1');
            p.addRequired('arg2');
            p.addOptional('arg3', []);
            p.addOptional('arg4', []);
            p.addParameter('IsRelative', true, @islogical);
            p.addParameter('Color', 'k');
            p.addParameter('Alpha', 0.15, @isnumeric);
            p.addParameter('Orientation', 'vertical', @(x) any(strcmp(x, {'vertical', 'horizontal'})));
            p.addParameter('Parent', [], @(x) isa(x, 'matlab.graphics.axis.Axes'));
            
            p.parse(varargin{:});
            arg1 = p.Results.arg1;
            arg2 = p.Results.arg2;
            arg3 = p.Results.arg3;
            arg4 = p.Results.arg4;
            if ~isempty(arg4)
                x = arg1;
                y = arg2;
                errPos = arg3;
                errNeg = arg4;
            elseif ~isempty(arg3)
                x = arg1;
                y = arg2;
                errPos = arg3;
                errNeg = arg3;
            else
                y = arg1;
                x = cumsum(ones(size(y)));
                errPos = arg2;
                errNeg = arg2;
            end
            
            if isvector(y)
                x = x(:);
                y = y(:);
                errPos = errPos(:);
                errNeg = errNeg(:);
            end
            
            isRelative = p.Results.IsRelative;
            color = p.Results.Color;
            faceAlpha = p.Results.Alpha;
            ori = p.Results.Orientation;
            ax = p.Results.Parent;
            
            if isempty(ax)
                ax = gca;
            end
            
            % Ploting
            for k = 1 : size(y,2)
                x = [x(:,k); flip(x(:,k))];
                if isRelative
                    err = [y(:,k)+errPos(:,k); flip(y(:,k)-errNeg(:,k))];
                else
                    err = [errPos(:,k); flip(errNeg(:,k))];
                end
                if strcmp(ori, 'vertical')
                    patch(ax, x, err, color, 'FaceAlpha', faceAlpha, 'LineStyle', 'none');
                else
                    patch(ax, err, x, color, 'FaceAlpha', faceAlpha, 'LineStyle', 'none');
                end
            end
        end
        
        function varargout = GroupRibbon(xRange, yRange, varargin)
            % Plot a segmented ribbon to label groups
            % 
            %   MPlot.GroupRibbon(xRange, yRange)
            %   MPlot.GroupRibbon(xRange, yRange, colors)
            %   MPlot.GroupRibbon(xRange, yRange, ..., 'IndVal', [])
            %   MPlot.GroupRibbon(xRange, yRange, ..., 'Groups', [])
            %   MPlot.GroupRibbon(xRange, yRange, ..., 'Style', 'patch')
            %   MPlot.GroupRibbon(xRange, yRange, ..., 'PlotArgs', {})
            %   [xCenter, yCenter] = MPlot.GroupRibbon(xRange, yRange, ...)
            % 
            % Inputs
            %   xRange          If the robbon goes along the y-axis, xRange is a 2-element numeric vector 
            %                   indicating where the width of the ribbon begins and ends. 
            %                   If the ribbon goes along the x-axis, xRange is
            %                   1) A vector with discrete values (i.e. input to MMath.ValueBonds function).
            %                   2) Edges of segments in an n-by-2 matrix where n is the number of segments.
            %                   3) An n-element cell array where each element is a m-by-2 matrix of segment 
            %                      edges. All segments in one matrix will be plotted using the same color.
            %   yRange          Same as xRange but for the other axis.
            %   colors          1) An n-by-3 or n-by-4 matrix of RGB or RGBA colors for the n groups.
            %                   2) A function handle of colormap. Defualt is the MATLAB lines function.
            %   'Groups'        If xRange or yRange is the 1) type. 'Groups' is useful to specify a subset 
            %                   or a superset of group values of interest. The default is the unique values 
            %                   present in the original vector. 
            %   'IndVal'        
            %   'Style'         Style of the ribbon, 'patch' or 'line'.
            %   'PlotArgs'      Arguments in a cell array for MATLAB patch or line function, depending on 
            %                   the choice of 'Style'.
            % Outputs
            %   xCenters, yCenters      Center coordinates of segments in cell array where each element 
            %                           corresponds to one group.
            
            % Parse user inputs
            p = inputParser;
            p.addOptional('Colors', [], @(x) isnumeric(x) || isa(x, 'function_handle'));
            p.addParameter('Groups', []);
            p.addParameter('IndVal', []);
            p.addParameter('Style', 'patch', @(x) ismember(lower(x), {'line', 'patch'}));
            p.addParameter('PlotArgs', {}, @iscell);
            p.parse(varargin{:});
            colors = p.Results.Colors;
            style = p.Results.Style;
            groups = p.Results.Groups;
            indVal = p.Results.IndVal;
            plotArgs = p.Results.PlotArgs;
            
            % Standardize ranges
            if isnumeric(xRange) && numel(xRange) == 2
                isVertical = true;
                yRange = findRange(yRange, groups, indVal);
                xRange = repmat({xRange(:)'}, size(yRange));
                [yCenters, xCenters] = findCenter(yRange, xRange);
            elseif isnumeric(yRange) && numel(yRange) == 2
                isVertical = false;
                xRange = findRange(xRange, groups, indVal);
                yRange = repmat({yRange(:)'}, size(xRange));
                [xCenters, yCenters] = findCenter(xRange, yRange);
            else
                error('Either xRange or yRange must be a 2-element numeric vector specifying the width of the ribbon');
            end
            
            function r = findRange(r, g, iv)
                if isvector(r)
                    % Find ranges indices from a vector of discrete values
                    r = MMath.ValueBounds(r, g, 'Uni', false);
                    for k = 1 : numel(r)
                        r{k} = r{k} + [-.5 .5];
                    end
                elseif ~iscell(r)
                    % Treat each row of range as a group
                    r = num2cell(r,2);
                end
                if ~isempty(iv)
                    % Convert indices to axis values
                    iv = iv(:);
                    ind = (1:numel(iv))';
                    for k = 1 : numel(r)
                        r{k} = interp1(ind, iv, r{k}, 'linear', 'extrap');
                    end
                end
            end
            
            function [cSeg, cWidth] = findCenter(rSeg, rWidth)
                cSeg = cell(size(rSeg));
                cWidth = cell(size(rWidth));
                for k = 1 : numel(rSeg)
                    cSeg{k} = mean(rSeg{k}, 2);
                    cWidth{k} = repmat(mean(rWidth{k}), size(cSeg{k}));
                    if isempty(cSeg{k})
                        cSeg{k} = NaN;
                        cWidth{k} = NaN;
                    end
                end
                if all(cellfun(@isscalar, cSeg))
                    cSeg = cell2mat(cSeg);
                    cWidth = cell2mat(cWidth);
                end
            end
            
            % Specify colors
            if isempty(colors)
                colors = lines(numel(xRange));
            elseif isa(colors, 'function_handle')
                colors = colors(numel(xRange));
            end
            
            % Plot the ribbon
            hold on
            for i = 1 : numel(xRange)
                if strcmpi(style, 'patch')
                    MPlot.Blocks(xRange{i}, yRange{i}, colors(i,:), plotArgs{:});
                elseif isVertical
                    w = diff(xRange{i});
                    y = [yRange{i} NaN(size(yRange{i},1),1)];
                    x = repmat(mean(xRange{i}), size(y));
                    line(x', y', 'Color', colors(i,:), 'LineWidth', w, plotArgs{:});
                else
                    w = diff(yRange{i});
                    x = [xRange{i} NaN(size(xRange{i},1),1)];
                    y = repmat(mean(yRange{i}), size(x));
                    line(x', y', 'Color', colors(i,:), 'LineWidth', w, plotArgs{:});
                end
            end
            
            % Output segment centers
            if nargout > 0
                varargout{1} = xCenters;
            end
            if nargout > 1
                varargout{2} = yCenters;
            end
        end
        
        function varargout = Figure(varargin)
            % Make figure with white background
            if isempty(varargin)
                f = gcf;
            else
                f = figure(varargin{:});
            end
            f.Color = 'w';
            if nargout > 0
                varargout{1} = f;
            end
        end
        
        function Paperize(varargin)
            % Make axes comply with conventions of publication
            %
            %   MPlot.Paperize(h)
            %   MPlot.Paperize(..., 'FontSize', 6)
            %   MPlot.Paperize(..., 'ColumnsWide', [])
            %   MPlot.Paperize(..., 'ColumnsHigh', [])
            %   MPlot.Paperize(..., 'AspectRatio', [])
            %   MPlot.Paperize(..., 'Zoom', 2)
            %   MPlot.Paperize(..., 'JournalStyle', 'Cell')
            %
            % Inputs:
            %   h           Array of Axes or Figure handle(s). Default is the current figure.
            %               If h is empty, all existing Axes will be operated on. 
            %   fontSize    Font size. Default 6.
            
            p = inputParser();
            p.addOptional('h', gcf, @ishandle);
            p.addParameter('FontSize', 6, @isscalar);
            p.addParameter('FontName', 'arial', @ischar);
            p.addParameter('Zoom', 2, @isscalar);
            p.addParameter('ColumnsWide', [], @isscalar);
            p.addParameter('ColumnsHigh', [], @isscalar);
            p.addParameter('AspectRatio', [], @isscalar);
            p.addParameter('JournalStyle', 'cell', @(x) any(strcmpi(x, {'nature', 'cell'})));
            
            p.parse(varargin{:});
            h = p.Results.h;
            fontSize = p.Results.FontSize;
            fontName = p.Results.FontName;
            z = p.Results.Zoom;
            colsWide = p.Results.ColumnsWide;
            colsHigh = p.Results.ColumnsHigh;
            aRatio = p.Results.AspectRatio;
            journalStyle = lower(p.Results.JournalStyle);
            
            switch journalStyle
                case 'nature'
                    widthSet = [8.9 12 18.3];
                case 'cell'
                    widthSet = [8.5 11.4 17.4];
            end
            
            % Resolve figure width
            if ~isempty(colsWide)
                % Calculate width by fold of cols
                figWidth = widthSet(1) * colsWide;
                
                % Overwrite if at specific #cols
                colOpts = [1 1.5 2];
                optIdx = colsWide == colOpts;
                if any(optIdx)
                    figWidth = widthSet(colsWide == colOpts);
                end
            end
            
            % 
            if isempty(h)
                h = findobj('Type', 'Axes');
            end
            
            for i = 1 : numel(h)
                if isa(h(i), 'matlab.ui.Figure')
                    if ~isempty(colsWide)
                        h(i).Color = 'w';
                        h(i).Units = 'centimeter';
                        h(i).Position(3) = figWidth * z;
                        if ~isempty(colsHigh)
                            h(i).Position(4) = figWidth / colsWide * colsHigh * z;
                        elseif ~isempty(aRatio)
                            h(i).Position(4) = figWidth * aRatio * z;
                        end
                    end
                    ax = findobj(h, 'Type', 'Axes');
                else
                    ax = h(i);
                end
                
                for j = 1 : numel(ax)
                    set(ax(j), ...
                        'TickDir', 'out', ...
                        'FontName', fontName, ...
                        'FontSize', fontSize*z, ...
                        'LabelFontSizeMultiplier', 1, ...
                        'TitleFontSizeMultiplier', 1);
                end
            end
        end
        
        function varargout = PlotPointAsLine(x, y, d, varargin)
            % Plot points as lines with specified length. It also accepts arguments of Line properties 
            % like 'line' function. 
            % 
            %   MPlot.PlotPointAsLine(x, y, d)
            %   MPlot.PlotPointAsLine(..., 'Orientation', 'vertical')
            %   MPlot.PlotPointAsLine(..., 'LinePropertyName', value)
            %   hh = MPlot.PlotPointAsLine(...)
            % 
            % Inputs
            %   x                   x coordinates of line centers. 
            %   y                   y coordinates of line centers. 
            %   d                   Length of lines. It can be a scalar or an array for each line. 
            %   'Orientation'       The orientation of lines. Default is 'vertical'.
            %   Any PropertyName-Value pair for MATLAB built-in 'line' function. 
            % 
            % Output
            %   hh                  Handles of Line objects
            
            p = inputParser();
            p.KeepUnmatched = true;
            p.addParameter('Orientation', 'vertical', @(x) any(strcmpi(x, {'vertical', 'horizontal'})));
            p.parse(varargin{:});
            orientationStr = p.Results.Orientation;
            varargin = p.Unmatched;
            
            x = x(:);
            y = y(:);
            d = d(:);
            
            if strcmpi(orientationStr, 'horizontal')
                xx = [x-d/2, x+d/2, NaN(size(x))]';
                yy = [y, y, NaN(size(y))]';
            else
                xx = [x, x, NaN(size(x))]';
                yy = [y-d/2, y+d/2, NaN(size(y))]';
            end
            
            xx = xx(:);
            yy = yy(:);
            
            hh = plot(xx, yy, varargin);
            
            if nargout == 1
                varargout{1} = hh;
            end
        end
        
        function varargout = PlotRaster(xx, varargin)
            % Plot raster
            % 
            %   MPlot.PlotRaster(xx)
            %   MPlot.PlotRaster(xx, yPos)
            %   MPlot.PlotRaster(xx, yPos, d)
            %   MPlot.PlotRaster(..., 'ColorArray', [])
            %   MPlot.PlotRaster(..., 'LinePropertyName', value)
            %   hh = MPlot.PlotRaster(...)
            % 
            % Inputs
            %   xx                  1) A numeric vector of x coordinates that plots one row.
            %                       3) A cell array of 1).
            %                       2) A matrix whose columns are 1).
            %   yPos                Y position of each row. 
            %   d                   Length of line segments. It can be a scalar or an array for each line. 
            %                       Will plot dots if set to zero (default). 
            %   'ColorArray'        1) An n-by-3 array of RGB colors. n is the number of rows. 
            %                       2) An n-by-4 array of RGBA colors. 'A'(alpha) controls transparency. 
            %                       3) An n-element char vector of colors. (e.g. 'k', 'r', 'm')
            %   Any PropertyName-Value pair for MATLAB built-in 'line' function. 
            % 
            % Output
            %   hh                  Handles of plotted line objects. 
            
            % Parse inputs
            p = inputParser();
            p.KeepUnmatched = true;
            p.PartialMatching = false;
            
            p.addRequired('xx');
            p.addOptional('yPos', []);
            p.addOptional('d', 0, @isnumeric);
            p.addParameter('ColorArray', [], @(x) isnumeric(x) || ischar(x));
            
            p.parse(xx, varargin{:});
            yPos = p.Results.yPos;
            d = p.Results.d;
            colorArray = p.Results.ColorArray;
            varargin = p.Unmatched;
            
            % Format xx
            if ~iscell(xx)
                if isvector(xx)
                    xx = xx(:);
                end
                xx = num2cell(xx, 1);
            end
            
            % Format yPos
            if isempty(yPos)
                yPos = 0 : numel(xx)-1;
            end
            yPos = yPos(:)';
            
            % Plot rasters
            for i = numel(xx) : -1 : 1
                x = xx{i};
                y = repmat(yPos(i), size(x));
                if d
                    hh(i) = MPlot.PlotPointAsLine(x, y, d, varargin);
                else
                    hh(i) = plot(x, y, '.', varargin);
                end
                if i == numel(xx)
                    hold on;
                end
            end
            
            % Apply colors
            if ~isempty(colorArray)
                if ischar(colorArray)
                    colorArray = colorArray(:);
                elseif isvector(colorArray)
                    colorArray = repmat(colorArray(:)', [numel(hh) 1]);
                end
                for i = 1 : numel(hh)
                    hh(i).Color = colorArray(i,:);
                end
            end
            
            % Output
            if nargout == 1
                varargout{1} = hh;
            end
        end
        
        function varargout = PlotTraceLadder(varargin)
            % Plot traces as a ladder
            % 
            %   MPlot.PlotTraceLadder(yy)
            %   MPlot.PlotTraceLadder(xx, yy)
            %   MPlot.PlotTraceLadder(xx, yy, yPos)
            %   MPlot.PlotTraceLadder(..., 'ColorArray', [])
            %   MPlot.PlotTraceLadder(..., 'LinePropertyName', value)
            %   hh = MPlot.PlotTraceLadder(...)
            % 
            % Inputs
            %   yy                  1) A numeric vector of y coordinates that plots one trace.
            %                       2) A matrix whose columns are 1).
            %                       3) A cell array of 1). Vectors do not need to have the same length.
            %   xx                  1) A vector of x coordinates that applies to all series in yy.
            %                       2) A matrix of 1) as columns for individual series in yy.
            %                       3) A cell array of 1) for individual series in yy.
            %                       4) An empty array []. Use sample indices as x coordinates. 
            %   yPos                Y position of each trace's zero after shifting them into a ladder. 
            %   'ColorArray'        1) An n-by-3 array of RGB colors. 
            %                       2) An n-by-4 array of RGBA colors. 'A'(alpha) controls transparency. 
            %                       3) An n-element char vector of colors. (e.g. ['k', 'r', 'm']')
            %                       When n equals the number of traces, colors are applied respectively. 
            %                       When n equals 1, the color applies to all traces (same as using the 
            %                       'Color' property). 
            %   Any PropertyName-Value pair for MATLAB line function. 
            % 
            % Output
            %   hh                  Handles of plotted line objects. 
            
            % Parse inputs
            p = inputParser();
            p.KeepUnmatched = true;
            p.PartialMatching = false;
            
            p.addRequired('arg1');
            p.addOptional('arg2', []);
            p.addOptional('yPos', []);
            p.addParameter('ColorArray', [], @(x) isnumeric(x) || ischar(x));
            p.addParameter('Scalar', 1, @isnumeric);
            
            p.parse(varargin{:});
            arg1 = p.Results.arg1;
            arg2 = p.Results.arg2;
            yPos = p.Results.yPos;
            colorArray = p.Results.ColorArray;
            r = p.Results.Scalar;
            varargin = p.Unmatched;
            
            if ~isempty(arg2)
                xx = arg1;
                yy = arg2;
            else
                xx = [];
                yy = arg1;
            end
            
            % Format yy
            if iscell(yy)
                cellfun(@(x) assert(isvector(x), 'Element of the cell array must be numeric vector'), yy);
                yy = PadNaN(yy);
            end
            if isvector(yy)
                yy = yy(:);
            end
            yy = yy .* r(:)';
            
            % Format xx
            if isempty(xx)
                xx = 1 : size(yy,1);
            end
            if iscell(xx)
                cellfun(@(x) assert(isvector(x), 'Element of the cell array must be numeric vector'), xx);
                xx = PadNaN(xx);
            end
            if isvector(xx)
                xx = repmat(xx(:), [1 size(yy,2)]);
            end
            
            % Format yPos
            if isempty(yPos)
                yPos = cumsum(-min(yy) + [0, max(yy(:,1:end-1))]);
            end
            yPos = yPos(:)';
            
            % Plot traces
            hh = plot(xx, yy+yPos, varargin);
            
            % Apply colors
            if ~isempty(colorArray)
                if (ischar(colorArray) && isscalar(colorArray)) || (isnumeric(colorArray) && isvector(colorArray))
                    colorArray = repmat(colorArray(:)', [numel(hh) 1]);
                end
                for i = 1 : numel(hh)
                    hh(i).Color = colorArray(i,:);
                end
            end
            
            % Output
            if nargout == 1
                varargout{1} = hh;
            end
            
            % Helper function
            function vOut = PadNaN(vIn)
                L = cellfun(@numel, vIn);
                for c = numel(vIn) : -1 : 1
                    vOut{c} = NaN(max(L), 1);
                    vOut{c}(1:L(c)) = vIn{c};
                end
                vOut = cell2mat(vOut);
            end
        end
        
        function cmap = PolarMap(varargin)
            % Colormap with zero-center white shading
            % Adapted from https://www.mathworks.com/matlabcentral/fileexchange/37099-polarmap-polarized-colormap
            % 
            % POLARMAP Polarized color map
            %   POLARMAP applies a "polarized" blue-white-red colormap to current figure,
            %	and adjusts the color axis limits to be centered to zero.
            %
            %	POLARMAP(M) fixes the number of colors to M (default is 64).
            %
            %	POLARMAP(MAP) applies linear shading to white to the center of colormap
            %	MAP which can be any of existing colormaps (an Mx3 matrix of RGB).
            %
            %	POLARMAP(MAP,C) uses exponent C to modify the shading contrast. Default
            %	is C = 1 for linear shading. Use C = 2 to strengthen the shading, or
            %	C = 0.5 to attenuate it.
            %
            %	C=POLARMAP(...) returns an M-by-3 matrix containing the colormap, that
            %	can be used with COLORMAP function like other colormaps.
            %
            %	Examples:
            %		pcolor(peaks), shading interp
            %		polarmap, colorbar
            %
            %	then try the following
            %		polarmap(jet,0.5)
            %
            %	Note the polar shading has no real interest with colormaps that include
            %	white color as one of the extremes (like GRAY, BONE, HOT, ...)
            
            % default parameters
            m = 64;	% number of colors
            c = 1;	% exponent of shading factor (1 = linear)
            
            if nargin > 0
                if ~isnumeric(varargin{1}) || (size(varargin{1},2) ~= 3 && ~isscalar(varargin{1}))
                    error('First argument must be numeric: scalar M or Mx3 color matrix');
                end
                if isscalar(varargin{1})
                    m = varargin{1};
                end
            end
            if nargin > 0 && size(varargin{1},2) == 3
                map = varargin{1};
                m = size(map,1);
            else
                map = bluered(m);
            end
            
            if nargin > 1 && isscalar(varargin{2})
                c = varargin{2};
            end
            
            % linear shading from min/max (colormap value) to center (white)
            r = repmat(abs(linspace(1,-1,m)).^c,[3,1])';
            map = map.*r + 1 - r;
            
            if nargout > 0
                cmap = map;
            else
                colormap(map)
                caxis([-1,1]*max(abs(caxis)))
                % Note: this fixes color axis to manual mode...
            end
            
            function map = bluered(m)
                
                if mod(m,2)
                    z = [0,0,0];
                    m2 = floor(m/2);
                else
                    z = zeros([0,3]);
                    m2 = m/2;
                end
                map = [repmat([0,0,1],[m2,1]);z;repmat([1,0,0],[m2,1])];
            end
        end
        
        function cc = Rainbow(numColors)
            % Returns a color map based on rainbow spectrum
            %
            %   cc = MPlot.Rainbow(numColors)
            %
            % Input:
            %   numColors       The number of colors to return
            % Output:
            %   cc              A series of color in numColors-by-3 matrix ranging from red to violet. 
            
            cBase = { [ 1 1 0 0 0 1 ]; [ 0 1 1 1 0 0 ]; [ 0 0 0 1 1 1 ] };
            cc = cellfun(@(x) interp1(1:6, x, linspace(1,6,numColors)), cBase, 'UniformOutput', false);
            cc = cell2mat(cc)';
        end
        
        function SavePDF(f, path, flag, res)
            % Saves figure as a pdf with the paper size clipped to the size of the
            % 
            %   SavePDF(f, path)
            %   SavePDF(f, path, flag)
            %   SavePDF(f, path, flag, res)
            % 
            %   f       Figure handle of figure to be saved
            % 	path    Name of file to be saved
            %   flag    flag = 0 => save as a regular pdf (default)
            %           flag = 1 => save with -zbuffer and -r flags
            %           Setting flag = 1 can help deal with blurry rendering of pdfs on mac
            %   res     (Only if flag = 1, default 300) resolution of picture
            % 
            % Robert Wilson
            % 18-Mar-2010
            
            if ~MPlot.isSave
                return
            end
            
            if ~exist('flag', 'var')
                flag = 0;
            end
            
            if ~exist('res', 'var')
                res = 300;
            end
            
            set(f, 'windowstyle', 'normal')
            set(f, 'paperpositionmode', 'auto')
            
            pp = get(f, 'paperposition');
            wp = pp(3);
            hp = pp(4);
            set(f, 'papersize', [wp hp])
            
            if flag
                print('-painters', f, '-dpdf', ['-r' num2str(res)], '-zbuffer', path);
            else
                print('-painters', f, '-dpdf', path);
            end
        end
        
        function SavePNG(f, path)
            % Save a figure as PNG file
            if ~MPlot.isSave
                return
            end
            print(f, path, '-dpng', '-r0');
        end
        
        function barLength = ScaleBar(data, heightRatio)
            % Return the length of a scale bar with appropriate length
            %
            %   barLength = MPlot.ScaleBar(data, heightRatio)
            %
            % Inputs:
            %   data                Numeric array. 
            %   heightRatio         The desired ratio of the length of scale bar to the span of data. 
            %                       Default ratio is 0.5.
            % Output:
            %   barLength           The length of scale bar in the same unit as data. 
            %
            
            if nargin < 2
                heightRatio = 0.5;
            end
            
            if isvector(data)
                data = data(:);
            end
            
            barLength = (max(data) - min(data)) * heightRatio;
            barLength = round(barLength, 1, 'significant');
        end
        
        function varargout = Violin(pos, bb, nn, varargin)
            % Plot histograms as violins
            % 
            %   MPlot.Violin(pos, bb, nn)
            %   MPlot.Violin(..., 'Color', 'k')
            %   MPlot.Violin(..., 'Style', 'patch')
            %   MPlot.Violin(..., 'Alpha', 0.3)
            %   MPlot.Violin(..., 'Orientation', 'vertical')
            %   MPlot.Violin(..., 'Alignment', 'center')
            %   MPlot.Violin(..., 'Percentiles', [])
            %   MPlot.Violin(..., 'PrctLineArgs', {'Color', 'r'})
            %   h = MPlot.Violin(...)
            %
            % Inputs:
            %   pos             An n-element vector indicating the position of n violin plots. 
            %   bb              A matrix where each column contains bin centers or edges of a violin.
            %   nn              A matrix where each column contains widths of a violin.
            %   'Color'         Face color of the violin plots. Default is black.
            %   'Style'         'patch' (default) or 'contour'.
            %   'Alpha'         Transparancy of the shade. Default 0.3.
            %   'Orientation'   Orientation of violins. Default is 'vertical'.
            %   'Alignment'     The side of violins aligned to straight line. Options include 'low', 
            %                   'high' and 'center' (default).
            %   'Percentiles'   Mark percentiles by providing a vector of percentages.
            %   'PrctLineArgs'  A cell array of addtional arguments for plotting percentiles.
            % Output: 
            %   hh              Handles of Patch object
            
            % Handles user inputs
            p = inputParser();
            p.addParameter('Color', 'k');
            p.addParameter('Style', 'patch', @(x) ismember(x, {'patch', 'contour'}));
            p.addParameter('Alpha', 1, @isnumeric);
            p.addParameter('Orientation', 'vertical', @(x) ismember(x, {'vertical', 'horizontal'}));
            p.addParameter('Alignment', 'center', @(x) ismember(x, {'low', 'high', 'center'}));
            p.addParameter('Percentiles', [], @isnumeric);
            p.addParameter('PrctLineArgs', {'Color', 'r'}, @iscell);
            
            p.parse(varargin{:});
            color = p.Results.Color;
            style = p.Results.Style;
            faceAlpha = p.Results.Alpha;
            ori = p.Results.Orientation;
            alignType = p.Results.Alignment;
            prct = p.Results.Percentiles;
            prctArgs = p.Results.PrctLineArgs;
            
            if size(bb,1) == size(nn,1)+1
                % Convert bin edges to bin centers
                bb = bb(1:end-1,:) + diff(bb,1,1);
            end
            
            if ~isempty(prct)
                % Compute histogram height at percentiles
                nnEps = nn;
                nnEps(~nnEps) = 10*eps; % add eps for monotonic CDFs
                nnCDF = cumsum(nnEps);
                nnCDF = nnCDF ./ nnCDF(end,:);
                bbPrct = zeros(numel(prct), size(bb,2));
                nnPrct = bbPrct;
                for i = 1 : size(bb,2)
                    bbPrct(:,i) = interp1(nnCDF(:,i), bb(:,i), prct/100);
                    nnPrct(:,i) = interp1(bb(:,i), nn(:,i), bbPrct(:,i));
                end
            else
                bbPrct = NaN(1, size(bb,2));
                nnPrct = bbPrct;
            end
            
            hold on;
            for k = size(nn,2) : -1 : 1
                b = bb(:,k);
                n = nn(:,k);
                
                % Trim zeros at the ends (necessary for 'contour' style)
                isTrim = cumsum(n, 'omitnan') == 0 | cumsum(n, 'omitnan', 'reverse') == 0;
                b(isTrim) = [];
                n(isTrim) = [];
                
                % Close the loop
                b = [b; flip(b)];
                bPrct = bbPrct(:,[k k]);
                nPrct = nnPrct(:,k);
                switch alignType
                    case 'center'
                        n = pos(k) + [-n/2; flip(n/2)];
                        nPrct = pos(k) + [-nPrct nPrct]/2;
                    case 'low'
                        n = pos(k) + [zeros(size(n)); flip(n)];
                        nPrct = pos(k) + [zeros(size(nPrct)) nPrct];
                    case 'high'
                        n = pos(k) + [-n; zeros(size(n))];
                        nPrct = pos(k) + [-nPrct zeros(size(nPrct))];
                end
                
                % Choose style
                if strcmp(style, 'patch')
                    violin = @(x,y) patch(x, y, color, 'FaceAlpha', faceAlpha, 'LineStyle', 'none');
                else
                    violin = @(x,y) line(x, y, 'Color', color);
                end
                
                % Plot a violin
                if strcmp(ori, 'vertical')
                    h(k) = violin(n, b);
                    line(nPrct', bPrct', prctArgs{:});
                else
                    h(k) = violin(b, n);
                    line(bPrct', nPrct', prctArgs{:});
                end
            end
            
            if nargout > 0
                varargout{1} = h;
            end
        end
    end
    
end

