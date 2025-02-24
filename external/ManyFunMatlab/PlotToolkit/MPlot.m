classdef MPlot
    %MPlot A collection of functions useful for plotting
    
    properties(Constant)
        isSave = true;
    end
    
    methods(Static)
        function alpha = AlphaForOverlap(N)
            % Find a good alpha transparency based on the number of overlapping traces
            % 
            %   alpha = MPlot.AlphaForOverlap(N)
            % 
            % Inputs
            %   N           The number of overlapping traces (or other objects).
            % Output
            %   alpha       An optimal alpha transparency.
            % 
            alpha = MMath.Bound(1/log2(N), [0 1]);
        end
        
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
                % Default face color is gray
                varargin = [{[.9 .9 .9]}, varargin];
            end
            if ~any(cellfun(@(x) (ischar(x) || isstring(x)) && startsWith('EdgeColor', x), varargin))
                % Default No edge
                varargin = [varargin {'EdgeColor', 'none'}];
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
            p = patch(xx, yy, varargin{:});
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
        
        function spArgs = FindSubplotInd(rows, cols, rowRange, colRange)
            % Return subplot indices based on the requested layout
            % 
            %   spArgs = MPlot.FindSubplotInd(nRows, nCols, rowRange, colRange)
            %   spArgs = MPlot.FindSubplotInd(rowDist, colDist, rowRange, colRange)
            % 
            % Inputs
            %   nRows, nCols        Same as the first two inputs to the subplot function.
            %                       Each row and column will be partitioned as a block.
            %   rowDist, colDist    Numeric vector of integer numbers specifying how to 
            %                       partition rows and columns into blocks following the 
            %                       same convention as the mat2cell function.
            %   rowRange, colRange  Vector of block row and column indiced indicating 
            %                       where the plot will be placed. Note that the rows and 
            %                       columns here refer to the partitioned blocks.
            % Output
            %   spArgs              A 1-by-3 cell array containing the three inputs to 
            %                       the subplot function, e.g. ax = subplot(spArgs{:})
            % 
            % See also subplot, mat2cell
            
            if isscalar(rows)
                rowDist = ones(rows, 1);
            else
                rowDist = rows;
            end
            if isscalar(cols)
                colDist = ones(cols, 1);
            else
                colDist = cols;
            end
            nRow = sum(rowDist);
            nCol = sum(colDist);
            
            if isempty(rowRange)
                rowRange = 1 : numel(rowDist);
            end
            if isempty(colRange)
                colRange = 1 : numel(colDist);
            end
            
            grid = mat2cell(zeros(nRow, nCol), rowDist, colDist);
            for i = rowRange(:)'
                for j  = colRange(:)'
                    grid{i,j}(:) = 1;
                end
            end
            grid = cell2mat(grid);
            spArgs = {nRow, nCol, find(grid')};
        end
        
        function ntArgs = FindTileInd(rows, cols, rowRange, colRange)
            % Return tile indices based on the requested layout
            % 
            %   ntArgs = MPlot.FindTileInd(nRows, nCols, rowRange, colRange)
            %   ntArgs = MPlot.FindTileInd(rowDist, colDist, rowRange, colRange)
            % 
            % Inputs
            %   nRows, nCols        Same as the first two inputs to the subplot function.
            %                       Each row and column will be partitioned as a block.
            %   rowDist, colDist    Numeric vector of integer numbers specifying how to 
            %                       partition rows and columns into blocks following the 
            %                       same convention as the mat2cell function.
            %   rowRange, colRange  Vector of block row and column indiced indicating 
            %                       where the plot will be placed. Note that the rows and 
            %                       columns here refer to the partitioned blocks.
            % Output
            %   ntArgs              A 1-by-2 cell array containing the two inputs to 
            %                       the nexttile function, e.g. ax = nexttile(ntArgs{:})
            % 
            % See also nexttile, mat2cell
            
            if isscalar(rows)
                rowDist = ones(rows, 1);
            else
                rowDist = rows;
            end
            if isscalar(cols)
                colDist = ones(cols, 1);
            else
                colDist = cols;
            end
            nRow = sum(rowDist);
            nCol = sum(colDist);
            
            if isempty(rowRange)
                rowRange = 1 : numel(rowDist);
            end
            if isempty(colRange)
                colRange = 1 : numel(colDist);
            end
            
            grid = mat2cell(zeros(nRow, nCol), rowDist, colDist);
            for i = rowRange(:)'
                for j  = colRange(:)'
                    grid{i,j}(:) = 1;
                end
            end
            grid = cell2mat(grid);
            
            k = find(grid', 1);
            cSpan = sum(any(grid, 1));
            rSpan = sum(any(grid, 2));
            
            ntArgs = {k, [rSpan cSpan]};
        end
        
        function Paperize(varargin)
            % Make axes comply with conventions of publication
            %
            %   MPlot.Paperize(h)
            %   MPlot.Paperize(h, ColumnsWide)
            %   MPlot.Paperize(h, ColumnsWide, ColumnsHigh)
            %   MPlot.Paperize(..., 'FontSize', 6)
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
            p.addOptional('colsWide', [], @isnumeric);
            p.addOptional('colsHigh', [], @isnumeric);
            p.addParameter('FontSize', 6, @isscalar);
            p.addParameter('FontName', 'arial', @ischar);
            p.addParameter('Zoom', 2, @isscalar);
            p.addParameter('ColumnsWide', [], @isscalar);
            p.addParameter('ColumnsHigh', [], @isscalar);
            p.addParameter('AspectRatio', [], @isscalar);
            p.addParameter('JournalStyle', 'cell', @(x) any(strcmpi(x, {'nature', 'cell'})));
            
            p.parse(varargin{:});
            h = p.Results.h;
            colsWide = p.Results.colsWide;
            colsHigh = p.Results.colsHigh;
            if isempty(colsWide)
                colsWide = p.Results.ColumnsWide;
            end
            if isempty(colsHigh)
                colsHigh = p.Results.ColumnsHigh;
            end
            fontSize = p.Results.FontSize;
            fontName = p.Results.FontName;
            z = p.Results.Zoom;
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
                
                % % Overwrite if at specific #cols
                % colOpts = [1 1.5 2];
                % optIdx = colsWide == colOpts;
                % if any(optIdx)
                %     figWidth = widthSet(colsWide == colOpts);
                % end
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
                    hold(hh(i).Parent, 'on');
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
        
        function PlotRasterStack(spk, varargin)
            % Plot rasters from a spikeTime table (or its cell array) in one axes with stacking units
            % 
            %   MPlot.PlotRasterStack(spk)
            %   MPlot.PlotRasterStack(spk, Y)
            %   MPlot.PlotRasterStack(..., 'HeightScale', 0.8)
            %   MPlot.PlotRasterStack(..., 'Color', [])
            %   MPlot.PlotRasterStack(..., 'LineWidth', .5)
            %   MPlot.PlotRasterStack(..., 'Parent', [])
            %   MPlot.PlotRasterStack(..., 'MarkUnits', [])
            %   MPlot.PlotRasterStack(..., 'MarkTrials', [])
            % 
            % Inputs
            %   spk             1) A trial-by-unit table or cell array of spike time vectors.
            %                   2) A unit-element cell array where each element is a nested trial-element 
            %                      cell array of spike time vectors. This allows different units to 
            %                      have different number of trials.
            %   Y               A numeric vector for each raster's middle Y position.
            %   'HeightScale'   The height of each raster set.
            %   'Color'         A unit-by-3 RGB or unit-by-4 RGBA array. If empty [], units alternate 
            %                   colors between [0 0 0 .7] and [.3 .3 .3 .7].
            %   'Parent'        Axes object to plot in.
            %   'LineWidth'     LineWidth of the spikes.
            %   'MarkUnits'     A vector of indices for units to be marked by dark red [.8 0 0 1].
            %   'MarkTrials'    A vector of indices for trials to be marked by dark red [.8 0 0 1].
            
            p = inputParser();
            p.addOptional('Y', 1:size(spk,2), @(x) isnumeric(x));
            p.addParameter('HeightScale', 0.8, @(x) isnumeric(x) && isscalar(x));
            p.addParameter('Color', [], @(x) true);
            p.addParameter('LineWidth', .5, @isnumeric);
            p.addParameter('Parent', [], @(x) isa(x, 'matlab.graphics.axis.Axes'));
            p.addParameter('MarkUnits', [], @isnumeric);
            p.addParameter('MarkTrials', [], @isnumeric);
            p.parse(varargin{:});
            Y = p.Results.Y;
            height_scale = p.Results.HeightScale;
            unit_c = p.Results.Color;
            indMU = p.Results.MarkUnits;
            indMT = p.Results.MarkTrials;
            lw = p.Results.LineWidth;
            ax = p.Results.Parent;
            
            % Standardize input to nested cell array
            if istable(spk)
                spk = spk{:,:};
            end
            if ~iscell(spk{1}) && ~isempty(spk{1})
                [n_trials, n_units] = size(spk);
                spk = mat2cell(spk, n_trials, ones(1,n_units));
            end
            n_units = numel(spk);
            
            % Determine colors
            if isempty(unit_c)
                unit_c = repmat([0 0 0 .7; .3 .3 .3 .7], [ceil(n_units/2) 1]);
            end
            if ~isempty(indMU)
                unit_c(indMU,:) = [.8 0 0 1]; % hightlight selected histograms in red
            end
            
            % Setup axes
            if isempty(ax)
                ax = gca;
            end
            hold(ax, 'on');
            
            for i = 1 : n_units
                u_spk = spk{i};
                n_trials = numel(u_spk);
                trial_height = 1/n_trials * height_scale;
                y = Y(i) - trial_height*n_trials/2 + trial_height/2;
                
                for j = 1 : n_trials
                    spk_t = u_spk{j};
                    spk_y = repelem(y, length(spk_t));
                    spk_h = trial_height * .8;
                    
                    if ~ismember(j, indMT)
                        MPlot.PlotPointAsLine(spk_t, spk_y, spk_h, 'Color', unit_c(i,:), 'LineWidth', lw, 'Parent', ax);
                    else
                        MPlot.PlotPointAsLine(spk_t, spk_y, spk_h, 'Color', [.8 0 0 1], 'LineWidth', lw, 'Parent', ax);
                    end
                    
                    y = y + trial_height;
                end
            end
            
            ax.YLim = [.5, n_units+.5];
            ax.YTick = 1 : n_units;
            ax.YDir = 'reverse';
            ax.XGrid = 'on';
        end
        
        function PlotHistStack(tt, hh, ee, varargin)
            % Plot a stack of histograms in one axes
            % 
            %   MPlot.PlotHistStack(tt, hh, ee)
            %   MPlot.PlotHistStack(..., 'Scaling', 1)
            %   MPlot.PlotHistStack(..., 'Style', 'trace')
            %   MPlot.PlotHistStack(..., 'Color', [])
            %   MPlot.PlotHistStack(..., 'Parent', [])
            % 
            % Inputs
            %   tt              1) m-element vector of timestamps, shared by all histograms.
            %                   2) m-by-n array of timestamps for the n different histograms.
            %   hh              m-by-n array for the bin heights of the n histograms.
            %   ee              m-by-n array for the error of the n histograms.
            %   'Scaling'       A scaling factor to adjust the height of the histograms.
            %   'Style'         'trace' or 'bar'.
            %   'Color'         n-by-3 RGB or n-by-4 RGBA array. If empty [], units alternate 
            %                   colors between [0 0 0] and [.3 .3 .3].
            %   'Parent'        Axes object to plot in.
            % 
            
            p = inputParser();
            p.KeepUnmatched = true;
            p.addParameter('Color', [], @(x) true);
            p.addParameter('Scaling', 1, @(x) isnumeric(x) && isscalar(x));
            p.addParameter('Style', 'trace', @(x) ismember(x, {'trace', 'bar'}));
            p.addParameter('MarkUnits', [], @isnumeric);
            p.addParameter('Parent', [], @(x) isa(x, 'matlab.graphics.axis.Axes'));
            p.parse(varargin{:});
            unit_c = p.Results.Color;
            alpha = 1;
            frac = p.Results.Scaling;
            style = p.Results.Style;
            indMU = p.Results.MarkUnits;
            ax = p.Results.Parent;
            
            if isempty(ax)
                ax = gca;
            end
            
            % Replicate timestamps for each histogram
            if isvector(tt)
                tt = tt(:);
            end
            if size(tt,2) == 1
                tt = repmat(tt, [1 size(hh,2)]);
            end
            
            % Scale heights
            hh = hh .* frac;
            ee = ee .* frac;
            
            % Determine colors
            n_units = size(hh, 2);
            if isempty(unit_c)
                unit_c = repmat([0 0 0; .3 .3 .3], [ceil(n_units/2) 1]);
            end
            if size(unit_c, 1) == 1
                unit_c = repmat(unit_c, [n_units 1]);
            end
            if ~isempty(indMU)
                unit_c(indMU,:) = [.8 0 0]; % hightlight selected histograms in red
            end
            
            % Plotting
            m = 0;
            y = frac/2; %0.5;
            hold(ax, 'on');
            
            for i = 1 : n_units
                t = tt(:,i);
                m = m + 1;
                switch style
                    case 'bar'
                        binEdges = MMath.BinCenters2Edges(t);
                        px = repelem(binEdges, 2);
                        py = [0 repelem(hh(:,i)',2) 0];
                        patch(ax, px, -py+y+1, unit_c(m,:), 'FaceAlpha', .1, 'EdgeColor', 'none', p.Unmatched);
                    case 'trace'
                        MPlot.ErrorShade(t, -hh(:,i)+y+1, ee(:,i), 'Color', unit_c(m,:), 'Alpha', 0.1, 'Parent', ax);
                        plot(ax, t, -hh(:,i)+y+1, 'Color', [unit_c(m,:) alpha], p.Unmatched);
                    otherwise
                        error('%s is not a supported style', style);
                end
                y = y + 1;
            end
            
            ax.YLim = [.5, n_units+.5];
            ax.YTick = 1 : n_units;
            ax.YDir = 'reverse';
            ax.XGrid = 'on';
        end
        
        function PlotHeatmapStack(T, MM, varargin)
            % Plot a stack of heatmaps
            % 
            %   MPlot.PlotHeatmapStack(T, MM)
            %   MPlot.PlotHeatmapStack(T, MM, ..., 'Scaling', 1)
            %   MPlot.PlotHeatmapStack(T, MM, ..., 'Parent', [])
            % 
            % Inputs
            %   T               s
            %   'Scaling'       s
            %   'Parent'        Axes object to plot in.
            % 
            
            % Handles inputs
            if isnumeric(MM)
                MM = mat2cell(MM, size(MM,1), size(MM,2), ones(size(MM,3),1));
            end
            
            n_maps = numel(MM);
            
            if ~isempty(T) && isnumeric(T)
                if isvector(T)
                    T = {T(:)};
                else
                    T = mat2cell(T, size(T,1), ones(size(T,2),1));
                end
            end
            if isscalar(T)
                T = repelem(T, n_maps);
            end
            
            p = inputParser();
            p.addParameter('Scaling', 1, @(x) isnumeric(x) && isscalar(x));
            p.addParameter('Parent', [], @(x) isa(x, 'matlab.graphics.axis.Axes'));
            p.parse(varargin{:});
            scaling = p.Results.Scaling;
            ax = p.Results.Parent;
            
            % Plot heatmaps
            hold(ax, 'on');
            for i = 1 : n_maps
                M = MM{i};
                [n_h, n_w] = size(M);
                
                if isempty(T)
                    t = (1:n_w)';
                else
                    t = T{i};
                end
                
                yBin = (0:n_h-1) / (n_h-1) - 0.5;
                yBin = yBin*scaling + i;
                
                imagesc(ax, t, yBin, flip(M,1));
            end
%             ax.CLim = [0 2];
%             ax.Colormap = brewermap([], 'Greys');
            ax.YLim = [.5, n_maps+.5];
            ax.YTick = 1 : n_maps;
            ax.YDir = 'reverse';
%             ax.XGrid = 'on';
            
%             ax.XLim = [min(tWin(:,1)) max(tWin(:,2))];
%             ax.XLabel.String = 'Time (s)';
            ax.FontSize = 8;
            MPlot.Axes(ax);
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
            %   yy                  1) A numeric vector of y coordinates for a single trace.
            %                       2) A matrix where each column is (1).
            %                       3) A cell array of (1). Vectors do not need to have the same length.
            %   xx                  1) A vector of x coordinates that applies to all series in yy.
            %                       2) A matrix of (1) as columns for individual series in yy.
            %                       3) A cell array of (1) for individual series in yy.
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
                dyPos = -min(yy) + [0, max(yy(:,1:end-1))];
                dyPos(isnan(dyPos)) = median(dyPos, 'omitnan');
                yPos = cumsum(dyPos);
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
        
        function lb = StaggerLabels(lb, offsets)
            % Satgger labels with added spaces
            %
            %   lb = StaggerLabels(lb)
            %   lb = StaggerLabels(lb, offsets)
            % 
            % Inputs
            %   lb          Label strings.
            %   offsets     A vector of character offsets.
            % Output
            %   lb          Staggered label strings.
            % 
            lb = cellstr(lb(:));
            if nargin < 2
                offsets = -[0 median(cellfun(@length, lb))];
            end
            if isscalar(offsets)
                offsets = [0 offsets];
            end
            offsets = round(offsets);
            offsets = repmat(offsets(:), ceil(numel(lb)/numel(offsets)), 1);
            offsets = offsets(1:numel(lb));
            sp = arrayfun(@(x) repelem(' ', x), abs(offsets), 'Uni', false);
            for i = 1 : numel(lb)
                if offsets(i) > 0
                    lb{i} = [sp{i} lb{i}];
                else
                    lb{i} = [lb{i} sp{i}];
                end
            end
            lb = string(lb);
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
        
        function h = ViolinScatter(pos0, val, varargin)
            % Plot data points arranged in the shape of a violin plot
            % 
            %   h = MPlot.ViolinScatter(pos0, val)
            %   h = MPlot.ViolinScatter(pos0, val, nbins)
            %   h = MPlot.ViolinScatter(pos0, val, edges)
            %   h = MPlot.ViolinScatter(pos0, val, ..., 'Span', 0.8)
            %   h = MPlot.ViolinScatter(pos0, val, ..., LineArgs)
            %   pos = MPlot.ViolinScatter(pos0, val, ..., 'IsPlot', false)
            % 
            % Inputs
            %   pos0            Center position of the plot.
            %   val             Values of the data points.
            %   nbins, edges    histcounts parameters that control the binning. For example, less bins create 
            %                   horizontally more spreadout and vertically more interspaced plot.
            %   'Span'          The maximum horizontal span.
            %   LineArgs        Any name-value arguments accepted by the MATLAB built-in plot or line function.
            % Output
            %   h               Handle of the plot object.
            % 
            
            % Parse user inputs
            p = inputParser;
            p.KeepUnmatched = true;
            p.addOptional('histParam', [], @(x) isnumeric(x));
            p.addParameter('Span', 0.8, @(x) isnumeric(x));
            p.addParameter('IsPlot', true, @islogical);
            p.parse(varargin{:});
            histParam = p.Results.histParam;
            sp = p.Results.Span;
            isPlot = p.Results.IsPlot;
            
            % Compute histogram
            if isempty(histParam)
                [N, ~, bins] = histcounts(val, 'Normalization', 'count');
            else
                [N, ~, bins] = histcounts(val, histParam, 'Normalization', 'count');
            end
            
            % Position data points
            ud = 1/max(N)*sp;
            pos = NaN(size(val));
            for i = 1 : numel(N)
                isIn = bins == i;
                binVal = val(isIn);
                [binVal, I] = sort(binVal);
                
                dPos = (1:ceil(numel(binVal)/2)) * ud;
                dPos = [0, repelem(dPos,2)];
                dPos(1:2:end) = -dPos(1:2:end);
                dPos = dPos(1:numel(binVal));
                
                dPos(I) = dPos;
                pos(isIn) = pos0 + dPos;
            end
            
            % Plot scatter
            if ~isPlot
                h = pos;
                return
            end
            h = plot(pos, val, '.', p.Unmatched);
        end
        
        
    end
    
end