classdef MotionPlot
    
    methods(Static)
        function AmpScaling(F)
            % Plot the fitted motion scaling function
            % 
            %   MotionPlot.AmpScaling(F)
            % 
            % Input
            %   F       A griddedInterpolant object that maps depth to motion magnitude.
            % 
            
            yq = (0 : 20 : 7660)';
            aq = F(yq);
            plot(aq, yq); hold on
            
            y = F.GridVectors{1};
            a = F.Values;
            plot(a, y, 'o');
            
            ax = MPlot.Axes(gca);
%             ax.XLim(1) = 0;
            ax.YLim = yq([1 end]);
            ax.XLabel.String = 'Amplitude (um)';
            ax.YLabel.String = 'Distance from tip (um)';
        end
        
        function [vidMat, Ye] = ExtrapolateTraces(Fscale, Y, t)
            % Make a movie of trace extrapolation
            % 
            %   [vidMat, Ye] = ExtrapolateTraces(Fscale, Y, t)
            % 
            
            % Initialize variables
            cc = lines(size(Y,2));
            vidMat = struct('cdata', [], 'colormap', []);
            frIdx = 1;
            
            % Prepare extrapolation
            [nTm, nTr] = size(Y);
            isVal = ~isnan(Y);
            len = sum(isVal)';
            Ye = Y;
            
            indTr = 1 : nTr;
            for i = indTr
                indTmp = setdiff(indTr, i);
                
                % Compute temporal distance (zero if overlap)
                dt = zeros(nTr, 1);
                for j = indTmp
                    bb = MMath.Logical2Bounds(isVal(:,i) | isVal(:,j));
                    if size(bb,1) == 1
                        dt(j) = 0;
                    else
                        dt(j) = bb(2,1) - bb(1,2);
                    end
                end
                
                % Sort templates first by temporal and then by trace length
                [~, order] = sortrows([dt len], {'ascend', 'descend'});
                
                % Extrapolate traces
                for j = order(:)'
                    y = MTracer.Motion.ExtrapolateTrace(Fscale, Ye(:,i), Y(:,j));
                    if sum(isnan(y)) == sum(isnan(Ye(:,i)))
                        continue
                    end
                    Ye(:,i) = y;
                    
                    cla;
                    plot(t, Ye); hold on
                    plot(t, Ye(:,i), 'Color', cc(i,:), 'LineWidth', 3);
                    ax = MPlot.Axes(gca);
                    ax.XLim = [250 1100];
                    ax.YLim = [500 5500];
                    ax.XLabel.String = 'Time (s)';
                    ax.YLabel.String = 'Distance from tip (um)';
                    ax.LooseInset = [0 0 0 0];
                    
                    drawnow;
                    vidMat(frIdx) = getframe(gcf);
                    frIdx  = frIdx + 1;
                    pause(0.5);
                end
            end
            
            vidMat = cat(4, vidMat.cdata);
        end
        
        function varargout = Traces(T, Y, varargin)
            % Plot traces on map
            % 
            %   PlotTraces(T, Y)
            %   PlotTraces(T, Y, cc)
            %   PlotTraces(T, Y, ..., 'ShowLabel', false)
            %   hh = PlotTraces(...)
            % 
            %   T       A vector or a cell array of such vectors for timestamps of traces.
            %   Y       A vector, a matrix of column vectors, or a cell array of them for the 
            %           cooresponding y coordinates.
            %   cc      A single character, or a 1-by-3 or #traces-by-3 array for color code
            %   hh      Handles of plotted traces
            % 
            
            % Make sure data are in cell array
            if ~iscell(T)
                T = {T};
            end
            if ~iscell(Y)
                Y = {Y};
            end
            
            % Parse inputs
            p = inputParser();
            p.KeepUnmatched = true;
            p.addOptional('cc', [], @(x) isnumeric(x) || (ischar(x) && isscalar(x)));
            p.addParameter('ShowLabel', false, @islogical);
            p.parse(varargin{:});
            cc = p.Results.cc;
            isLabel = p.Results.ShowLabel;
            lineArgs = p.Unmatched;
            
            % Set default line style
            if ~isfield(lineArgs, 'LineWidth')
                lineArgs.LineWidth = 2;
            end
            if ~isfield(lineArgs, 'HitTest')
                lineArgs.HitTest = 'off';
            end
            
            % Plot traces
            nTr = numel(T);
            hh = cell(nTr, 1);
            for i = 1 : nTr
                if isempty(Y{i})
                    continue
                end
                hh{i} = plot(T{i}, Y{i}, lineArgs); hold on
            end
            hh = cat(1, hh{:});
            
            % Apply style
            nTr = numel(hh);
            if isempty(cc)
                cc = lines(nTr);
            end
            if size(cc,1) == 1
                cc = repmat(cc, [nTr 1]);
            end
            for i = 1 : nTr
                hh(i).Color = cc(i,:);
                if isLabel
                    x = hh(i).XData;
                    y = hh(i).YData;
                    k = find(~isnan(y), 1, 'first');
                    text(x(k), y(k), num2str(i), 'Color', cc(i,:), ...
                        'FontSize', 16, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
                end
            end
            
            if nargout == 0
                ax = MPlot.Axes(gca);
                ax.XLabel.String = 'Time (sec)';
                ax.YLabel.String = 'Distance from tip (um)';
                ax.LooseInset = [0 0 0 0];
            else
                varargout{1} = hh;
            end
        end
        
        function varargout = Spikes(tb)
            % Plot spikes on map
            
            t = tb.time;
            y = tb.y;
            a = tb.amp;
            ax = gca;
            
            ampRange = 8 : 100;
            hh = cell(numel(ampRange), 1);
            for k = 1 : numel(ampRange)
                % for each amplitude bin, plot all the spikes of that size in the
                % same shade of gray
                ind = a == ampRange(k); % the amplitudes are rounded to integers
                if ~any(ind)
                    continue
                end
                hh{k} = plot(ax, t(ind), y(ind), ...
                    'LineStyle', 'none', ...
                    'Marker', '.', ...
                    'MarkerSize', 4, ...
                    'Color', [1 1 1] * max(0, 1-ampRange(k)/40), ... % the marker color here has been carefully tuned
                    'HitTest', 'off');
                hold on
            end
            
            if nargout == 0
                ax = MPlot.Axes(gca);
                ax.XLabel.String = 'Time (sec)';
                ax.YLabel.String = 'Distance from tip (um)';
                ax.LooseInset = [0 0 0 0];
            else
                varargout{1} = hh;
            end
        end
        
    end
end

