classdef MTracerTrace < handle
    % 
    
    properties
        vm MTracerVM;
        dataTb table;       % anchor t, y, and undo status
        hAnchor;            % handle of anchor points
        hInterp;            % handle of interpolated trace
        motionFig;          % handle of motion figure
        
        tEps = 0.2;         % temporal tolerance of selection in seconds
        yEps = 100;         % spatial tolerance of selection in microns
        tGap = Inf;         % interval greater than this value separates segments
        tSample = 0.02;     % sample time for interpolation
    end
    
    properties(Dependent)
        numPts;
        dispName;
        fileName;
    end
    
    methods
        function val = get.numPts(this)
            val = sum(this.dataTb.isDone);
        end
        function val = get.dispName(this)
            [t, y] = GetCoords(this);
            tStr = [num2str(round(t(1))) '-' num2str(round(t(end))) 's'];
            yStr = [num2str(round(mean(y,'omitnan'))) 'um'];
            val = ['t: ' tStr ', y: ' yStr];
        end
        function val = get.fileName(this)
            [t, y] = GetCoords(this);
            tStr = [num2str(round(t(1))) '-' num2str(round(t(end))) 's'];
            yStr = [num2str(round(mean(y,'omitnan'))) 'um'];
            val = ['MTracerTrace_' this.vm.recId '_' tStr '_' yStr];
        end
    end
    
    methods
        % Object construction
        function this = MTracerTrace(vm, tb)
            % Constructor
            this.vm = vm;
            if nargin < 2
                tb = table();
                tb.t = NaN;
                tb.y = NaN;
                tb.isDone = false;
            end
            this.dataTb = tb;
        end
        
        function obj = Duplicate(this, newVM)
            obj = MTracerTrace(newVM, this.dataTb);
        end
        
        function delete(this)
            % 
            delete(this.hAnchor);
            delete(this.hInterp);
            delete(this.motionFig);
        end
        
        % Interactions
        function SetPoint(this, tPt, yPt)
            % Add a new point to table
            
            T = this.dataTb;
            
            % Overwrite existing point if having duplicate time
            I = tPt == T.t;
            
            % Append if unique
            if ~any(I)
                T{end+1,:} = zeros(1, width(T));
                I = height(T);
            end
            
            T.t(I) = tPt;
            T.y(I) = yPt;
            T.isDone(I) = true;
            
            % Clear history
            T(~T.isDone,:) = [];
            
            this.dataTb = T;
            
            this.PlotTrace();
        end
        
        function RemovePoint(this, tPt, yPt)
            % Remove the point that's closest to [tPt, yPt] in time, within limits
            
            % Find the closest point within limits
            T = this.dataTb;
            T.t(~T.isDone) = Inf; % exclude undone points by making them Inf
            [tdMin, k] = min(abs(T.t - tPt));
            ydMin = abs(T.y(k) - yPt);
            
            if tdMin < this.tEps && ydMin < this.yEps
                % Undo the point and sort it to the first of all the undones
                T = this.dataTb;
                T.isDone(k) = false;
                T = sortrows(T, 'isDone', 'descend');
                this.dataTb = T;
                
                this.PlotTrace();
            end
        end
        
        function [tPt, yPt] = AutoSetPoint(this)
            % Add a new point to table
            tracer = this.vm.tracers.(this.vm.currentTracer);
            [tTr, yTr] = GetCoords(this);
            [tPt, yPt] = tracer.Query(tTr(end), yTr(end));
            this.SetPoint(tPt, yPt);
        end
        
        function Undo(this)
            % Undo the last point
            idx = find(this.dataTb.isDone, 1, 'last');
            if ~isempty(idx)
                this.dataTb.isDone(idx) = false;
                this.PlotTrace();
            end
        end
        
        function Redo(this)
            % Redo the last point
            idx = find(~this.dataTb.isDone, 1, 'first');
            if ~isempty(idx)
                this.dataTb.isDone(idx) = true;
                this.PlotTrace();
            end
        end
        
        % Data operations
        function [t, y, ti, yi] = GetCoords(this)
            % Get t and y coordinates
            
            tb = sortrows(this.dataTb, 't');
            t = tb.t(tb.isDone);
            y = tb.y(tb.isDone);
            
            if numel(t) > 1
                isGap = [0; diff(t)] > this.tGap;
                segInd = 1 + cumsum(isGap);
                for i = max(segInd) : -1 : 1
                    ind = segInd == i;
                    tSeg = t(ind);
                    ySeg = y(ind);
                    if sum(ind) > 1
                        tiSeg = (tSeg(1) : this.tSample : tSeg(end))';
                        yiSeg = interp1(tSeg, ySeg, tiSeg, 'makima');
                    else
                        tiSeg = tSeg;
                        yiSeg = ySeg;
                    end
                    if i ~= max(segInd)
                        ti{i} = [tiSeg; NaN];
                        yi{i} = [yiSeg; NaN];
                    else
                        ti{i} = tiSeg;
                        yi{i} = yiSeg;
                    end
                end
                ti = cat(1, ti{:});
                yi = cat(1, yi{:});
%                 ti = (t(1) : this.tSample : t(end))';
%                 yi = interp1(t, y, ti, 'makima');
            else
                ti = t;
                yi = y;
            end
        end
        
        function ydMin = Dist2Selection(this, tPt, yPt)
            % Find the Y distance of the cloest point to trace
            
            if this.numPts == 0
                ydMin = NaN;
                return
            end
            
            % Find the closest point
            T = this.dataTb;
            T.t(~T.isDone) = Inf; % exclude undone points by making them Inf
            [tdMin, k] = min(abs(T.t - tPt));
            ydMin = abs(T.y(k) - yPt);
            
            % Check if the closest point is out of limits
            if tdMin > this.tEps || ydMin > this.yEps
                ydMin = NaN;
            end
        end
        
        % Plotting
        function PlotTrace(this, isSelected)
            % Plot anchor points and the interpolated trace
            
            if nargin < 2
                isSelected = true;
            end
            
            % Get data
            if ~this.numPts
                return
            end
            [t, y, ti, yi] = this.GetCoords();
            
            % Select line width
            if isSelected
                w = 4;
            else
                w = 2;
            end
            
            % Plot new
            ax = this.vm.hAxes;
            if isempty(this.hAnchor) || isempty(this.hInterp) || ~isvalid(this.hAnchor) || ~isvalid(this.hInterp)
                cc = lines(7);
                plotArgs = {'Color', cc(5,:), 'HitTest', 'off'};
                this.hAnchor = plot(ax, t, y, 'o', plotArgs{:});
                this.hInterp = plot(ax, ti, yi, '-', 'LineWidth', w, plotArgs{:});
                this.vm.AddHandle2Layer(this.hAnchor, 'anchors');
                this.vm.AddHandle2Layer(this.hInterp, 'interp');
            else
                set(this.hAnchor, 'XData', t, 'YData', y);
                set(this.hInterp, 'XData', ti, 'YData', yi, 'LineWidth', w);
            end
        end
        
        function PlotMotions(this)
            % 
            
            if isempty(this.motionFig) || ~isvalid(this.motionFig)
                this.motionFig = MPlot.Figure( ...
                    'Name', ['MTracer: Inspection: trace ' this.dispName], ...
                    'NumberTitle', 'off', ...
                    'IntegerHandle', 'off');
            else
                clf(this.motionFig);
            end
            
            % Get coordinates
            [t, y, ti, yi] = this.GetCoords();
            ti(isnan(yi)) = [];
            yi(isnan(yi)) = [];
            
            % Highpass filter to get fluctuation
            hpFilt = designfilt('highpassiir', 'FilterOrder', 8, ...
                'PassbandFrequency', 0.1, 'PassbandRipple', 0.2, ...
                'SampleRate', 1/this.tSample);
            yHp = filtfilt(hpFilt, yi);
            
            % Compute instantaneous amplitude and frequency
            c = hilbert(yHp);
            amp = abs(c);
            
            % Smooth the instantanous amplitude
%             lpFilt = designfilt('lowpassiir', 'FilterOrder', 8, ...
%                 'PassbandFrequency', 0.1, 'PassbandRipple', 0.2, ...
%                 'SampleRate', 1/this.tSample);
%             amp = filtfilt(lpFilt, amp);
            amp = smooth(amp, 5/this.tSample);
            
            % Make plot
            figure(this.motionFig);
            
            ax = subplot(3,1,1); cla
            plot(ax, ti, yi, 'k'); hold on
            ax.XLabel.String = 'Time (s)';
            ax.YLabel.String = 'Distance from tip (um)';
            ax.Title.String = 'Trace (original)';
            ax.XLim = ti([1 end]);
            ax.XGrid = 'on';
            ax.YGrid = 'on';
            ax.XMinorGrid = 'on';
            ax.YMinorGrid = 'on';
            MPlot.Axes(ax);
            
            ax = subplot(3,1,2); cla
            MPlot.ErrorShade(ti, yi, amp, -amp, 'IsRelative', false, 'Alpha', 0.1); hold on
            plot(ax, ti, yHp);
            ax.XLabel.String = 'Time (s)';
            ax.YLabel.String = 'Amplitude (um)';
            ax.Title.String = 'Trace (highpassed at 0.1Hz) and amplitude (5s boxcar)';
            ax.XLim = ti([1 end]);
            ax.XGrid = 'on';
            ax.YGrid = 'on';
            ax.XMinorGrid = 'on';
            ax.YMinorGrid = 'on';
            MPlot.Axes(ax);
            
            ax = subplot(3,1,3); cla
            pspectrum(yHp, ti, 'FrequencyLimits', [0 5]);
            ax.YLim = [-20 60];
            ax.Title.String = 'Power spectrum of highpassed trace';
            MPlot.Axes(ax);
        end
    end
    
end
