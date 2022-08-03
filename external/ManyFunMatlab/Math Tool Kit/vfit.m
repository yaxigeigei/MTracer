classdef vfit < handle
    %VFIT Volume Fit Object
    % 
    %   Detailed explanation goes here
    
    properties
        dataVolume = [];
        dataPoints = [];
        axisVals = {};
        method = 'nearest';
    end
    
    methods
        function this = vfit(v, varargin)
            % Handles user inputs
            this.dataVolume = v;
            
            p = inputParser();
            p.addOptional('axisVals', {});
            p.addParameter('method', 'nearest', @(x) any(strcmp(x, {'linear', 'nearest'})));
            p.parse(varargin{:});
            this.axisVals = p.Results.axisVals;
            this.method = p.Results.method;
            
            % Converts data volume to data points
            numDims = sum(size(v) > 1);
            ndgridStr = 'ndgrid(';
            for i = 1 : numDims
                ndgridStr = [ ndgridStr, '1:', num2str(size(v,i)), ', ' ];
            end
            ndgridStr = [ ndgridStr(1:end-2), ')' ];
            
            coor = cell(1, numDims);
            [ coor{:} ] = eval(ndgridStr);
            if isempty(this.axisVals)
                this.axisVals = arrayfun(@(x) (1:x)', size(v), 'Uni', false);
            else
                coor = cellfun(@(x,i) x(i), this.axisVals(1:length(coor)), coor, 'Uni', false);
            end
            coor = cell2mat(cellfun(@(x) x(:), coor, 'UniformOutput', false));
            
            v = [coor, v(:)];
            this.dataPoints = v;
        end
        
        function yi = feval(this, varargin)
            % Handles user inputs
            p = inputParser();
            p.addRequired('xi');
            p.addParameter('axisVals', this.axisVals);
            p.parse(varargin{:});
            xi = p.Results.xi;
            this.axisVals = p.Results.axisVals;
            
            if strcmpi(this.method, 'nearest') && ~isempty(this.axisVals)
                for i = length(this.axisVals)-1 : -1 : 1
                    dists = arrayfun(@(x) abs(this.axisVals{i}-x), xi(:,i), 'Uni', false);
                    [~, nearestInd{i}] = cellfun(@min, dists);
                end
                nearestInd = sub2ind(cellfun(@length, this.axisVals), nearestInd{:});
                yi = this.dataVolume(nearestInd);
            else
                yi = griddataVolumen(this.dataPoints(:,1:end-1), this.dataPoints(:,end), xi, this.method);
            end
        end
        
        function plot(this)
            scatter3(this.dataPoints(:,1), this.dataPoints(:,2), this.dataPoints(:,3), ...
                (MMath.Normalize(this.dataPoints(:,4)) * 64) + 1, this.dataPoints(:,4), 'fill');
        end
    end
    
end

