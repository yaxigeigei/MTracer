classdef Motion
    
    properties(Constant)
        sd2ampFunc = @(x) x / 0.7171 * 2;  % this converts SD to pk2pk amplitude assuming sine wave
    end
    
    methods(Static)
        function tb = ResampleTraces(tb, t)
            % Resample traces over the time window
            
            % Get times to resample at
            tAll = cat(1, tb.time{:});
            if nargin < 2 || isempty(t)
                t = min(tAll) : 0.02 : max(tAll);
            end
            t = t(:);
            
            % Resample using MSessionExplorer
            se = MSessionExplorer();
            se.SetTable('traces', tb, 'timeSeries');
            dt = t(2) - t(1);
            tEdges = [t-dt/2; t(end)+dt/2];
            tb = se.ResampleTimeSeries('traces', tEdges);
            
            % Denest table by making all y's in a matrix
            t = tb.time{1};
            Y = cat(2, tb.y{:});
            tb = table;
            tb.time = t;
            tb.Y = Y;
        end
        
        function tbSeg = SegmentTraceCoverage(tb, allowGap)
            % Partition resampled traceTable by unique combinations of traces
            
            if nargin < 2
                allowGap = true;
            end
            
            % Use NaN in yy as bitcode for segmentation
            Y = tb.Y{1};
            bc = ~isnan(Y);
            id = sum(bc.*2.^(size(Y,2):-1:1), 2); % id = bit2int(bc', size(bc,1));
            if ~allowGap
                id = cumsum(abs([0; diff(id)]));
            end
            [~, ~, id] = unique(id, 'stable');
            
            tbSeg = table;
            for i = 1 : width(tb)
                colName = tb.Properties.VariableNames{i};
                tbSeg.(colName) = splitapply(@(x) {x}, tb.(colName){1}, id);
            end
            for i = 1 : height(tbSeg)
                ind = ~isnan(tbSeg.Y{i}(1,:));
                tbSeg.Y{i} = tbSeg.Y{i}(:,ind);
                tbSeg.ind{i} = ind;
            end
        end
        
        function F = FitInterpolant2(t, Y)
            % Compute 2d makima interpolant from traces
            % 
            %   F = FitInterpolant2(t, Y)
            % 
            % t         A vector of m timestamps.
            % Y         A m-by-n matrix of depth where n is the number of traces.
            % F         Computed interpolant that outputs the change of y.
            % 
            
            % If there's only one trace, make a shifted duplicate
            if size(Y,2) == 1
                Y = [Y Y+1];
            end
            
            % Ensure depths are monotonically incresing (as required by griddedInterpolant)
            y0 = Y(1,:);
            [y0, I] = sort(y0);
            Y = Y(:,I);
            
            % Fit forward interpolant: dy(t) = F(t, y(t0))
            [T, Y0] = ndgrid(t, y0);
            dY = Y - Y0;
            F = griddedInterpolant(T, Y0, dY, 'makima');
        end
        
        function [fFluc, fDrift] = FitAmpScaling(Y, fs, method)
            % Fit two models that scale motion amplitude as a function of depth
            % 
            %   [fFluc, fDrift] = FitAmpScaling(Y, fs, method)
            % 
            % Inputs
            %   Y           A m-by-n matrix of depth coordinates. m is the number of time points. n is the number of traces.
            %   fs          Sampling frequency in Hz.
            %   method      Method of interpolation. Use 'makima' for nonlinear models and 'linear' for linear models.
            % Outputs
            %   fFluc       A griddedInterpolant object that maps the depth to the amplitude of fluctuation.
            %   fDrift      A griddedInterpolant object that maps the depth to the size of drift.
            % 
            % See also griddedInterpolant
            
            if nargin < 3
                method = 'makima';
            end
            
            % Remove samples with incomplete set of Y values
            hasNaN = any(isnan(Y), 2);
            Y(hasNaN,:) = [];
            
            % Sort traces with ascending depth
            [~, I] = sort(Y(1,:));
            Y = Y(:,I);
            
            % Compute rapid fluctuation and the residual motion
            hpFilt = designfilt('highpassiir', 'FilterOrder', 8, ...
                'PassbandFrequency', 0.2, 'PassbandRipple', 0.2, ...
                'SampleRate', fs);
            hpY = filtfilt(hpFilt, Y);
            resY = Y - hpY;
            
            % Fit fluctuation scaling
            ym = median(Y)';
            sd = std(hpY)';
            amp = MTracer.Motion.sd2ampFunc(sd);
            if strcmp(method, 'linear')
                fitObj = fit(ym, amp, 'poly1');
                amp = fitObj(ym);
            end
            fFluc = griddedInterpolant(ym, amp, 'makima', 'linear');
            
            % Drift scaling
            sd = std(resY)';
            if strcmp(method, 'linear')
                fitObj = fit(ym, sd, 'poly1');
                sd = fitObj(ym);
            end
            fDrift = griddedInterpolant(ym, sd, 'makima', 'linear');
        end
        
        function Ye = ExtrapolateTraces(Y, varargin)
            % Extrapolate traces
            % 
            %   Ye = ExtrapolateTraces(Y, F1fluc)
            %   Ye = ExtrapolateTraces(Y, F1fluc, F1drift)
            % 
            % Inputs
            %   Y           A m-by-n matrix of depth coordinates. m is the number of time points. n is the number of traces.
            %   fFluc       A griddedInterpolant object that maps the depth to the amplitude of fluctuation.
            %   fDrift      A griddedInterpolant object that maps the depth to the size of drift.
            % Output
            %   Ye          A m-by-n matrix of extrapolated depth coordinates.
            % 
            
            [nTm, nTr] = size(Y);
            isVal = ~isnan(Y);
            len = sum(isVal)';
            depth = median(Y, 'omitnan');
            
            Ye = Y;
            indTr = 1 : nTr;
            for i = indTr
                indTmp = setdiff(indTr, i);
                
                % Find temporal distances between the trace to extend and other template traces (zero if overlap)
                dt = zeros(nTr, 1);
                for j = indTmp
                    bb = MMath.Logical2Bounds(isVal(:,i) | isVal(:,j));
                    if size(bb,1) == 1
                        dt(j) = 0;
                    else
                        dt(j) = bb(2,1) - bb(1,2);
                    end
                end
                
                % Find spatial distances between the trace to extend and other template traces
                dd = zeros(nTr, 1);
                for j = indTmp
                    dd(j) = abs(depth(j) - depth(i));
                end
                
                % Sort templates in the order of temporal distance, trace length, and spatial distance
                % Ignore trance length if overlapping
                tb = table(indTr', dt, len, dd);
                tb.len(dt == 0) = Inf;
                [tb, order] = sortrows(tb, {'dt', 'len', 'dd'}, {'ascend', 'descend', 'ascend'});
                
                % Extrapolate traces
                for j = order(:)'
                    Ye(:,i) = MTracer.Motion.ExtrapolateTrace(Ye(:,i), Y(:,j), varargin{:});
%                     plot(Ye);
%                     xlim([0 nTm+1]);
%                     ylim([0 7660]);
                end
            end
        end
        
        function yExt = ExtrapolateTrace(yExt, yBase, varargin)
            % Extrapolate a single trace
            % 
            %   yExt = ExtrapolateTrace(yExt, yBase)
            %   yExt = ExtrapolateTrace(yExt, yBase, fFluc)
            %   yExt = ExtrapolateTrace(yExt, yBase, fFluc, fDrift)
            % 
            % Inputs
            %   yExt        A vector of depth coordinates. Elements that are NaN in yExt but not NaN in yBase will be extrapolated.
            %   yBase       A vector of depth coordinates.
            %   fFluc       1) A griddedInterpolant object that maps the depth to the amplitude of fluctuation. The ratio 
            %                  between the fluctuation amplitudes of extrapolated and base trace will be used for scaling.
            %               2) A numeric scalar that will be directly applied to the scaling. We can use 1 to disable 
            %                  (or maintain the same) fluctuation scaling.
            %               Default is 1.
            %   fDrift      1) A griddedInterpolant object that maps the depth to the size of drift.The ratio between the 
            %                  drift magnitudes of extrapolated and base trace will be used for scaling.
            %               2) A numeric scalar that will be directly applied to the scaling. We can use 1 to disable 
            %                  (or maintain the same) drift scaling.
            %               Default is 1.
            % Output
            %   yExt        A vector of extrapolated depth coordinates.
            % 
            
            nExt = size(yExt, 2);
            if nExt > 1
                for i = 1 : nExt
                    yExt(:,i) = MTracer.Motion.ExtrapolateTrace(yExt(:,i), yBase, varargin{:});
                end
                return
            end
            
            switch numel(varargin)
                case 0
                    fFluc = 1;
                    fDrift = 1;
                case 1
                    fFluc = varargin{1};
                    fDrift = fFluc;
                case 2
                    [fFluc, fDrift] = varargin{:};
            end
            
            % Compute scaling factors
            if isnumeric(fFluc) && isscalar(fFluc)
                sFluc = fFluc;
            else
                sFluc = fFluc(median(yExt, 'omitnan')) / fFluc(median(yBase, 'omitnan'));
            end
            sFluc = max(sFluc, 0); % ensure non-negative
            
            if isnumeric(fDrift) && isscalar(fDrift)
                sDrift = fDrift;
            else
                sDrift = fDrift(median(yExt, 'omitnan')) / fDrift(median(yBase, 'omitnan'));
            end
            sDrift = max(sDrift, 0); % ensure non-negative
            
            % Find segments of base to copy
            isExt = ~isnan(yExt);
            isBase = ~isnan(yBase);
            isCopy = isBase & ~isExt;
            bbCopy = MMath.Logical2Bounds(isCopy);
            
            for i = 1 : size(bbCopy,1)
                % Get the copied segment
                b1 = bbCopy(i,1);
                b2 = bbCopy(i,2);
                yCopy = yBase(b1:b2);
                
                % Scale movement
                if numel(varargin) < 2
                    yCopy = (yCopy - median(yCopy)).*sFluc;
                else
                    hpFilt = designfilt('highpassiir', 'FilterOrder', 8, ...
                        'PassbandFrequency', 0.2, 'PassbandRipple', 0.2, ...
                        'SampleRate', 50); % hardcoding
                    if numel(yCopy) <= 24
                        yHP = zeros(size(yCopy));
                    else
                        yHP = filtfilt(hpFilt, yCopy);
                    end
                    yRes = yCopy - yHP;
                    yCopy = yHP.*sFluc + (yRes - median(yRes)).*sDrift;
                end
                
                % Match the joints
                iPrev = MMath.Bound(b1-1, [1 numel(isExt)]);
                iNext = MMath.Bound(b2+1, [1 numel(isExt)]);
                
                if isExt(iPrev) && ~isExt(iNext)
                    % Extending trace forward
                    yCopy = yCopy - yCopy(1) + yExt(iPrev);
                    
                elseif ~isExt(iPrev) && isExt(iNext)
                    % Extending trace backward
                    yCopy = yCopy - yCopy(end) + yExt(iNext);
                    
                elseif isExt(iPrev) && isExt(iNext)
                    % Filling a gap in trace
                    yCopy = interp1(yCopy([1 end])', yExt([iPrev iNext])', yCopy);
                end
                
                yExt(b1:b2) = yCopy;
            end
        end
        
    end
end

