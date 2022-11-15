classdef CNN < MTracer.TracerBaseClass
    %CNN Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        dataTb;                 % 
        X
        Y
        Yhat
        cvp;                    % crossvalidation partitioning object
        net;
    end
    
    methods
        function [t2, y2] = Query(this, t1, y1)
            % Use the CNN to query points
            img = this.GetImages(t1, y1);
            dy2 = predict(this.net, img);
            [dt2, dy2] = this.Frac2Map(dy2);
            t2 = t1 + dt2;
            y2 = y1 + dy2;
        end
        
        function img = GetImages(this, t, y)
            % 
            
            if nargin < 3
                t = this.vm.focus(1);
                y = this.vm.focus(2);
            end
            
            if ~isscalar(t)
                img = arrayfun(@(a,b) this.GetImages(a,b), t, y, 'Uni', false);
                img = cat(4, img{:});
                return
            end
            
            [tWin, yWin] = this.FindMapWindows(t, y);
            s = this.vm.SliceMap(tWin, yWin, 'DataNames', {'traces', 'spikes', 'LFP'}, 'TraceType', 'interp');
            s.focus = [t y];
            
            img = cat(3, ...
                this.GetSpikeImage(s), ...
                this.GetLFPImage(s), ...
                this.GetTraceImage(s));
        end
        
        function PrepareData(this)
            % 
            
            rng(61);
            
            tb = table;
            [~, ~, tb.t, tb.y] = arrayfun(@(x) x.GetCoords(), this.vm.traces, 'Uni', false);
            
            % Randomly sample from the data
            for i = 1 : height(tb)
                t = tb.t{i};
                y = tb.y{i};
                
                % Compute target number of samples
                dIndResp = round(this.dtPred / this.vm.traces(i).tSample);
                dur = t(end) - t(1);
                nSp = round(dur*5); % on avg 5 samples/sec
                nSrc = numel(t)-dIndResp;
                
                % Initialize weight
                dy = gradient(y(1:nSrc)) * 10;
                w = ones(size(dy));
                
                % Use the speed of motion as sampling weight
                w = w .* abs(dy).^2;
                
%                 % Use hardcoded function to compute weight
%                 w = w .* min((1+2*dy./0.5).^2, 10);
                
                % 
                ind = randsample(nSrc, nSp, true, w);
                ind = sort(ind);
                tb.tSamp{i} = t(ind);
                tb.ySamp{i} = y(ind);
                tb.img{i} = this.GetImages(tb.tSamp{i}, tb.ySamp{i});
                
                tb.tResp{i} = tb.t{i}(ind+dIndResp);
                tb.yResp{i} = tb.y{i}(ind+dIndResp);
                tb.dyResp{i} = (tb.yResp{i} - tb.ySamp{i}) / this.mapSize(2); % normalize
            end
            
            % Response normalization by augmenting samples
            oriX = cat(4, tb.img{:});
            oriY = cat(1, tb.dyResp{:});
            
%             edges = -0.5:0.1:0.5;
%             [N, ~, binInd] = histcounts(oriY, edges);
%             
%             isOOB = binInd == 0;
%             oriX = oriX(:,:,:,~isOOB);
%             oriY = oriY(~isOOB);
%             binInd = binInd(~isOOB);
%             
%             augN = min(round(max(N)./N), 10); % maximally augment 10 times
%             augX = cell(size(N));
%             augY = cell(size(N));
%             for i = 1 : numel(N)
%                 isBin = binInd == i;
%                 [augX{i}, augY{i}] = AugmentSample(augN(i), oriX(:,:,:,isBin), oriY(isBin));
%             end
%             
%             function [imgAug, dyAug] = AugmentSample(nAug, img, dy)
%                 % 
%                 imgAug = repmat(img, [1 1 1 nAug]);
%                 dyAug = repmat(dy, [nAug 1]);
%             end
            
            this.dataTb = tb;
            this.X = oriX;
            this.Y = oriY;
%             this.X = cat(4, augX{:});
%             this.Y = cat(1, augY{:});
            this.Yhat = NaN(size(this.Y));
            this.cvp = cvpartition(numel(this.Y), 'KFold', 10);
        end
        
        function TrainResnet(this)
            % Train the CNN with training samples
            
            ind = this.cvp.training(1);
            XTrain = this.X(:,:,:,ind);
            YTrain = this.Y(ind);
            
            ind = this.cvp.test(1);
            XValid = this.X(:,:,:,ind);
            YValid = this.Y(ind);
            
            miniBatchSize  = 64;
            options = trainingOptions('sgdm', ...
                'MiniBatchSize', miniBatchSize, ...
                'MaxEpochs', 60, ...
                'InitialLearnRate', 1e-4, ...
                'LearnRateSchedule', 'piecewise', ...
                'LearnRateDropFactor', 0.1, ...
                'LearnRateDropPeriod', 30, ...
                'Shuffle', 'every-epoch', ...
                'ValidationData', {XValid, YValid}, ...
                'ValidationFrequency', floor(numel(this.Y)/miniBatchSize), ...
                'Plots', 'training-progress', ...
                'CheckpointPath', 'C:\chang_lab\presentations\2022-11-15 SfN2022 minisymposium\checkpoints', ...
                'Verbose', false);
            
            if isempty(this.net)
                nn = resnet18;
                nn = layerGraph(nn);
%                 nn = removeLayers(nn, {'input_1', 'fc1000', 'fc1000_softmax', 'ClassificationLayer_fc1000'});
                nn = removeLayers(nn, {'data', 'fc1000', 'prob', 'ClassificationLayer_predictions'});
                
                layersInput = [
                    imageInputLayer([this.imgSize 3], 'Name', 'image')
                    resize2dLayer('OutputSize', [224 224], 'Name', 'resize224')];
                nn = addLayers(nn, layersInput);
                nn = connectLayers(nn, 'resize224', 'conv1');
                
                numResp = 1;
                layersReg = [
                    fullyConnectedLayer(numResp, 'WeightLearnRateFactor', 20, 'BiasLearnRateFactor', 20, 'Name', 'fc_reg')
                    regressionMAELayer('regMAE')];
                nn = addLayers(nn, layersReg);
                nn = connectLayers(nn, 'pool5', 'fc_reg');
                
            else
                nn = this.net;
            end
            
            this.net = trainNetwork(XTrain, YTrain, nn, options);
        end
        
        function TrainAlexnet(this)
            % Train the CNN with training samples
            
            ind = this.cvp.training(1);
            XTrain = this.X(:,:,:,ind);
            YTrain = this.Y(ind);
            
            ind = this.cvp.test(1);
            XValid = this.X(:,:,:,ind);
            YValid = this.Y(ind);
            
            miniBatchSize  = 128;
            options = trainingOptions('sgdm', ...
                'MiniBatchSize', miniBatchSize, ...
                'MaxEpochs', 60, ...
                'InitialLearnRate', 1e-4, ...
                'LearnRateSchedule', 'piecewise', ...
                'LearnRateDropFactor', 0.1, ...
                'LearnRateDropPeriod', 30, ...
                'Shuffle', 'every-epoch', ...
                'ValidationData', {XValid, YValid}, ...
                'ValidationFrequency', floor(numel(this.Y)/miniBatchSize), ...
                'Plots', 'training-progress', ...
                'CheckpointPath', 'C:\chang_lab\presentations\2022-11-15 SfN2022 minisymposium\checkpoints', ...
                'Verbose', false);
            
            if isempty(this.net)
                layers = alexnet('Weights', 'none');
                numResp = 1;
                layersInput = [
                    imageInputLayer([this.imgSize 3], 'Name', 'image')
                    resize2dLayer('OutputSize', [227 227], 'Name', 'resize227')];
                layersTransfer = layers(2:end-3);
                layersReg = [
%                     fullyConnectedLayer(numResp, 'WeightLearnRateFactor', 20, 'BiasLearnRateFactor', 20, 'Name', 'fc_reg')
                    fullyConnectedLayer(numResp, 'Name', 'fc_reg')
                    regressionMAELayer('regMAE')];
                layers = [
                    layersInput
                    layersTransfer
                    layersReg];
%                 layersTransfer = nn.Layers(1:end-3);
%                 numResp = 1;
%                 layersReg = [
%                     fullyConnectedLayer(numResp, 'WeightLearnRateFactor', 20, 'BiasLearnRateFactor', 20, 'Name', 'fc_reg')
%                     regressionMAELayer('regMAE')];
%                 layers = [layersTransfer; layersReg];
                
            else
                layers = this.net.Layers;
            end
            
            this.net = trainNetwork(XTrain, YTrain, layers, options);
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
            tracer = MTracer.CNN();
            tracer.net = this.net;
            if nargin < 2
                filePath = 'tracer_alexnet0';
            end
            save(filePath, 'tracer');
        end
    end
end

