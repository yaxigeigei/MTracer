classdef Polygon < handle
    %POLYGON Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        mt              % handle of the main app view-model
        coords          % a n-by-2 numeric array that stores the (t,y) coordinates of the cutting polygon
        hAxes           % handle of the current axes to plot in
        hPlot           % handle of the polygon plot
    end
    properties(Dependent)
        numPoints
        hasPlot
        target
    end
    methods
        function val = get.numPoints(this)
            val = size(this.coords, 1);
        end
        function tf = get.hasPlot(this)
            tf = ~isempty(this.hPlot) && isvalid(this.hPlot);
        end
        function val = get.target(this)
            h = this.hAxes;
            if ~isempty(h) && isvalid(h)
                val = h.Tag;
            else
                val = '';
            end
        end
    end
    
    methods
        function this = Polygon(mt)
            %POLYGON Construct an instance of this class
            this.mt = mt;
        end
        
        function delete(this)
            this.ClearPolygon();
        end
        
        function AddPoint(this, ax, tPt, yPt)
            % Add a new point
            
            % Clear the existing polygon if users start setting point in a different window
            if isempty(this.hAxes) || ax ~= this.hAxes
                this.ClearPolygon();
                this.hAxes = ax;
            end
            
            % Append the point
            this.coords = [this.coords; tPt yPt];
            
            % Make or update plot
            this.PlotPolygon();
        end
        
        function RemoveLastPoint(this)
            % Remove the last point
            if size(this.coords, 1) > 0
                this.coords(end,:) = [];
                this.PlotPolygon();
            end
        end
        
        function PlotPolygon(this)
            % Plot or update polygon
            
            if size(this.coords, 1) < 1
                return
            end
            
            t = this.coords([1:end 1],1);
            y = this.coords([1:end 1],2);
            
            if ~this.hasPlot
                this.hPlot = plot(this.hAxes, t, y, '-', 'LineWidth', 1, 'Color', [0 0 1]);
            else
                set(this.hPlot, 'XData', t, 'YData', y);
            end
        end
        
        function ClearPolygon(this)
            % Clear the vertices coordinates and the polygon plot
            this.coords = [];
            delete(this.hPlot);
        end
        
        function inMask = IsInPolygon(this, xq, yq)
            % Return a logical mask with Trues indicating the input points encircled by the polygon
            if size(this.coords, 1) < 3
                inMask = false(size(xq));
                return
            end
            xv = this.coords(:,1);
            yv = this.coords(:,2);
            inMask = inpolygon(xq, yq, xv, yv);
        end
        
    end
end

