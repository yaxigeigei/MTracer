classdef LayeredFigure < handle
    % 
    
    properties
        hFig                    % handle to the figure window
        hAxes                   % handle to the axes in the figure window
        layers                  % stores plot elements in the axes, each field is a layer
        
        name = ''               % name of the object
        figPos                  % last changed hFig.Position
        axPos                   % last changed hAxes.Position
    end
    
    properties(Dependent)
        isOpen                  % whether the current window is open
        layerNames              % names of the layers
        positions               % last changed figure and axes positions in a struct
    end
    methods
        function val = get.isOpen(this)
            h = this.hFig;
            val = ~isempty(h) && ishandle(h) && isvalid(h);
        end
        function val = get.layerNames(this)
            if ~isempty(this.layers)
                val = fieldnames(this.layers);
            else
                val = {};
            end
        end
        function s = get.positions(this)
            s.figPos = this.figPos;
            s.axPos = this.axPos;
        end
        function set.positions(this, s)
            this.figPos = s.figPos;
            this.axPos = s.axPos;
            this.ApplyLayout();
        end
    end
    
    methods
        % Construction
        function this = LayeredFigure(layerNames)
            % Constructor
            if nargin > 0
                for i = 1 : numel(layerNames)
                    this.layers.(layerNames{i}) = [];
                end
            end
        end
        
        function f = Duplicate(this)
            f = MTracer.LayeredFigure(this.layerNames);
            f.name = this.name;
            f.figPos = this.figPos;
            f.axPos = this.axPos;
        end
        
        function delete(this)
            arrayfun(@(x) delete(x.hFig), this);
        end
        
        % Display
        function [f, ax] = Open(this, varargin)
            % Create figure and axes
            if ~this.isOpen
                f = MPlot.Figure(varargin{:});
                f.SizeChangedFcn = @this.FigureSizeChanged;
                ax = axes(f);
                this.hFig = f;
                this.hAxes = ax;
                this.ApplyLayout();
            else
                f = this.hFig;
                ax = this.hAxes;
            end
        end
        
        function ApplyLayout(this)
            % Apply figure and axes layout
            if ~this.isOpen
                return
            end
            if ~isempty(this.axPos)
                this.hAxes.Position = this.axPos;
            end
            if ~isempty(this.figPos)
                this.hFig.Position = this.figPos;
            end
        end
        
        function FigureSizeChanged(this, src, eventdata)
            % Keep track of the figure and axes positions
            this.figPos = src.Position;
            this.axPos = src.Children.Position;
        end
        
        function ToggleLayer(this, num)
            % Toggle layer on and off
            if num < 1 || num > numel(this.layerNames)
                fprintf("Cannot toggle layer #%s because it does not exist.\n", num);
                return
            end
            hh = this.layers.(this.layerNames{num});
            for i = 1 : numel(hh)
                hh(i).Visible = ~hh(i).Visible;
            end
        end
        
        function AddHandle2Layer(this, h, layerName)
            % Add plot object handle(s) to a layer
            
            if ~isfield(this.layers, layerName) || isempty(this.layers.(layerName))
                % Add object to the new or empty layer
                this.layers.(layerName) = h;
            else
                % Append object handle to the end
                this.layers.(layerName)(end+1:end+numel(h), 1) = h;
                this.ClearInvalidHandles();
            end
        end
        
        function ClearInvalidHandles(this)
            % Remove invalid object handles in every layers
            n = fieldnames(this.layers);
            for i = 1 : numel(n)
                if ~isempty(this.layers.(n{i}))
                    isValid = isvalid(this.layers.(n{i}));
                    this.layers.(n{i})(~isValid) = [];
                end
            end
        end
        
        function DeleteLayers(this, varargin)
            % Delete the specified layers and clear the plotted objects
            for i = 1 : numel(varargin)
                layerName = varargin{i};
                if isfield(this.layers, layerName)
                    delete(this.layers.(layerName));
                    this.layers = rmfield(this.layers, layerName);
                end
            end
        end
        
    end
    
end