classdef TracerBaseClass < handle
    %CNN Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        vm MTracerVM;
        mapSize= [4 500];       % time in seconds and distance in microns
        imgSize = [100 100];    % width and height
        dtPred = 0.2;           % amount of time in sec from focus to predict
    end
    
    properties(Dependent)
        imgRes;
    end
    
    methods
        function val = get.imgRes(this)
            val = this.mapSize ./ this.imgSize;
        end
    end
    
    methods
        function this = TracerBaseClass(vm)
            %CNN Construct an instance of this class
            if nargin > 0
                this.vm = vm;
            end
        end
        
        function img = GetSpikeImage(this, s)
            % 
            
            t = s.spikes.time;
            y = s.spikes.y;
            amp = s.spikes.amp;
            
            [tWin, yWin] = this.FindMapWindows(s.focus(1), s.focus(2));
            tEdges = linspace(tWin(1), tWin(2), this.imgSize(1)+1);
            yEdges = linspace(yWin(1), yWin(2), this.imgSize(2)+1);
            
            [~, ~, ~, yBin, tBin] = histcounts2(y, t, yEdges, tEdges);
            img = accumarray([yBin tBin], amp, this.imgSize, @sum);
        end
        
        function img = GetLFPImage(this, s)
            % 
            img = imresize(flip(s.LFP.v'), this.imgSize);
        end
        
        function img = GetTraceImage(this, s)
            % 
            
            t = cat(1, s.traces.time{:});
            y = cat(1, s.traces.y{:});
            y(t > s.focus(1)) = NaN;
            
            [tWin, yWin] = this.FindMapWindows(s.focus(1), s.focus(2));
            tEdges = linspace(tWin(1), tWin(2), this.imgSize(1)+1);
            yEdges = linspace(yWin(1), yWin(2), this.imgSize(2)+1);
            
            img = histcounts2(y, t, yEdges, tEdges);
        end
        
        function [tWin, yWin] = FindMapWindows(this, t, y)
            % 
            
            tHalf = this.mapSize(1)/2;
            tWin = [t-tHalf, t+tHalf];
            
            yHalf = this.mapSize(2)/2;
            yWin = [y-yHalf, y+yHalf];
        end
        
        function [dt, dy] = Frac2Img(this, dy)
            % Convert network prediction unit to pixel coordinates on image
            dt = repmat(this.dtPred, size(dy));
            dt = dt / this.mapSize(1) * this.imgSize(1);
            dy = dy * this.imgSize(2);
        end
        
        function [dt, dy] = Img2Map(this, dt, dy)
            % Convert pixel coordinates on image to relative map coordinates
            if isempty(dt)
                dt = repmat(this.dtPred, size(dy));
            else
                dt = dt / this.mapSize(1) * this.imgRes(1);
            end
            dy = dy * this.imgRes(2);
        end
        
        function [dt, dy] = Frac2Map(this, dy)
            % Convert network prediction unit to relative map coordinates
            [~, dy] = this.Frac2Img(dy);
            [dt, dy] = this.Img2Map([], dy);
        end
        
        function SaveAsTracer(this, filePath)
            % 
            tracer = MTracer.CNN();
            tracer.net = this.net;
            if nargin < 2
                filePath = 'tracer_alexnet0';
            end
            save(filePath, 'tracer');
        end
        
    end
end

