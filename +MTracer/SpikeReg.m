classdef SpikeReg < MTracer.TracerBaseClass
    
    properties
        
    end
    
    methods
        function [t2, y2] = Query(this, t1, y1)
            % query the next point
            
%             img = this.GetSpikeImage(t1, y1);
%             
%             dy2 = predict(this.net, img);
            
            [x, y, a] = this.GetSpikes(t1, y1);
            dx1 = x - t1;
            dy1 = y - y1;
            mdlr = fitlm(dx1, dy1, 'y ~ x1', 'RobustOpts', 'on');
            
            dy2 = predict(mdlr, this.dtPred);
            
            t2 = t1 + this.dtPred;
            y2 = y1 + dy2;
        end
        
        function [tSpk, ySpk, aSpk] = GetSpikes(this, t, y)
            % 
            
            if nargin < 3
                t = this.vm.focus(1);
                y = this.vm.focus(2);
            end
            
            if ~isscalar(t)
                [tSpk, ySpk, aSpk] = arrayfun(@(a,b) this.GetSpikes(a,b), t, y, 'Uni', false);
                tSpk = cat(1, tSpk{:});
                ySpk = cat(1, ySpk{:});
                aSpk = cat(1, aSpk{:});
                return
            end
            
            [~, yWin] = this.FindMapWindows(t, y);
            tWin = [t t+this.dtPred];
            s = this.vm.SliceMap(tWin, yWin, 'DataNames', {'spikes'});
            s.focus = [t y];
            tSpk = s.spikes.time;
            ySpk = s.spikes.y;
            aSpk = s.spikes.amp;
        end
        
        function ComputeValidation(this)
            % Validate the CNN with validation samples
            
            ind = this.cvp.test(1);
            XValid = this.X(:,:,:,ind);
%             YValid = this.Y(ind);
            this.Yhat(ind) = predict(this.net, XValid);
        end
        
        function f = PlotSample(this, k)
            % Plot input images and response position
            
            img = this.X(:,:,:,k);
            % t = linspace(tWin(1), tWin(2), this.imgRes(1));
            % y = linspace(yWin(1), yWin(2), this.imgRes(2));
            imgName = ["Spike cloud", "Lowpassed LPF", "Trace history"];
            
            pos0 = this.imgSize/2;

            [dtPx, dyPx] = this.Frac2Img(this.Y(k));
            pos = [dtPx dyPx] + this.imgSize/2;
            
            [dtPx, dyPx] = this.Frac2Img(this.Yhat(k));
            posHat = [dtPx dyPx] + this.imgSize/2;

            f = MPlot.Figure(456); clf
            for i = 1 : size(img,3)
                ax = subplot(1,3,i);
                % imagesc(ax, t, y, img(:,:,i));
                imagesc(ax,  img(:,:,i));
                colormap(ax, MPlot.PolarMap);
                hold on
                plot(pos0(1), pos0(2), 'k+', 'MarkerSize', 6, 'LineWidth', 2);
                plot(pos(1), pos(2), 'g*', 'MarkerSize', 6, 'LineWidth', 1);
%                 plot(posHat(1), posHat(2), 'g*');
                axis xy equal off
                ax.Title.String = imgName(i);
                MPlot.Axes(ax);
            end
        end
        
        function f = PlotExamplePred(this, k)
            % Plot input images and response position
            
            img = this.X(:,:,:,k);
            [h, w, ~, ~] = size(img);
            x = (1:w) - round(w/2);
            y = (1:h) - round(h/2);
            
            [dtPx, dyPx] = this.Frac2Img(this.Y(k));
            [dtPxH, dyPxH] = this.Frac2Img(this.Yhat(k));
            
            f = MPlot.Figure(789); clf
            N = numel(k);
            r = ceil(sqrt(N));
            c = ceil(N/r);
            r = 3;
            c = 6;
            for i = 1 : N
                ax = subplot(r, c, i);
                imagesc(ax, x, y, img(:,:,1,i));
                colormap(ax, 'gray');
                hold on
                
                plot(dtPx(i), dyPx(i), 'go', 'LineWidth', 1);
                plot(dtPxH(i), dyPxH(i), 'y*', 'LineWidth', 1);
                
                qArgs = {'AutoScale', 'off', 'LineWidth', 1, 'MaxHeadSize', 0.5, 'ShowArrowHead', 'off'};
                quiver(0, 0, dtPx(i), dyPx(i), 'Color', 'g', qArgs{:});
                quiver(0, 0, dtPxH(i), dyPxH(i), 'Color', 'y', qArgs{:});
                
                axis xy off
                ax.Title.String = "#" + k(i);
                MPlot.Axes(ax);
            end
            tightfig(f);
        end
        
        function f = PlotPerf(this)
            % Plot performance
            
            ind = this.cvp.test(1);
            [~, dy] = this.Frac2Map(this.Y(ind));
            [~, dyHat] = this.Frac2Map(this.Yhat(ind));
            
            err = dyHat - dy;
            v = dy ./ this.dtPred / 1e3; % to millisec
            
            T = MNeuro.Tuning(v, err);
            T = T{1};
            
            f = MPlot.Figure(136); clf
            plot(v, err, 'k.'); hold on
            MPlot.ErrorShade(T(:,1), T(:,2), T(:,3));
            plot(T(:,1), T(:,2));
            ax = gca;
            ax.XLabel.String = "Motion velocity (um/ms)";
            ax.YLabel.String = "Error (um)";
            MPlot.Axes(ax);
        end
        
        function SaveAsTracer(this, filePath)
            % 
            tracer = MTracer.SpikeReg();
            if nargin < 2
                filePath = 'tracer_spkreg';
            end
            save(filePath, 'tracer');
        end
    end
end

