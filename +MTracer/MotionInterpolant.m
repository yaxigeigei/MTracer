classdef MotionInterpolant
    
    properties
        t;
        Y;
        y0;
        fwF;        % forward interpolant, i.e. dY_observed = F(t, Y_stable)
        bwF;        % backward interpolant, i.e. dY_stable = F(t, Y_observed)
    end
    
    methods
        function this = MotionInterpolant(t, Y, t0)
            % 
            
            if ~exist('t0', 'var')
                t0 = 0;
            end
            
            % Remove NaNs
            hasNaN = any(isnan(Y), 2);
            nNaN = sum(hasNaN);
            if nNaN
                if all(hasNaN)
                    error('Cannot create interpolant because no sample (row) in Y has complete set of Y values');
                else
                    warning('%i samples (rwo) in Y do not have a complete set of values', nNaN);
                end
            end
            t = t(~hasNaN);
            Y = Y(~hasNaN,:);
            
            % If there's only one trace, make a shifted duplicate
            if size(Y,2) == 1
                Y = [Y Y+1];
                warning('Only one anchor trace is available. Make a shifted duplicate for rigid interpolation.');
            end
            
            % Ensure depths are monotonically incresing (as required by griddedInterpolant)
            [~, I] = min(abs(t-t0));
            y0 = Y(I,:);
            [y0, I] = sort(y0);
            Y = Y(:,I);
            
            % Fit forward interpolant
            [T, Y0] = ndgrid(t, y0);
            dY = Y - Y0;
            F = griddedInterpolant(T, Y0, dY, 'makima', 'linear');
            
            this.t = t;
            this.Y = Y;
            this.y0 = y0;
            this.fwF = F;
        end
        
        function [T, Y, dY] = ComputeDispField(this, t, y0)
            
            if nargin < 2 || isempty(t)
                t = this.t;
            end
            if nargin < 3
                y0 = (0 : 200 : 7660)';
            end
            [T, Y0] = ndgrid(t, y0);
            dY = this.fwF(T, Y0);
            Y = Y0 + dY;
        end
        
        function [Yk, dYk] = CorrectDepths(this, tq, Yq)
            % Apply backward interpolation to find the stable depths
            % 
            %   [Yk, dYk] = obj.CorrectDepths(tq, Yq)
            % 
            % Inputs
            %   tq              A m-element vector of sample times in second
            %   Yq              A m-by-n matrix or a n-element row vector of depths
            % Output
            %   Yk              The corrected depths
            %   dYk             The amount of correction
            
            tq = tq(:);
            if isrow(Yq)
                Yq = repmat(Yq, size(tq));
            end
            
            % Remove out of range samples
            ind = tq >= this.t(1) & tq <= this.t(end);
            if ~any(ind)
                Yk = Yq;
                dYk = zeros(size(Yq));
                return
            end
            tqFull = tq;
            YqFull = Yq;
            tq = tq(ind);
            Yq = Yq(ind,:);
            
            % Get anchor positions at query times
            for i = size(this.Y,2) : -1 : 1
                Ya(:,i) = interp1(this.t, this.Y(:,i), tq, 'makima');
            end
            
            % Predict queried y(t=1) given landmark Y(t=i), Y(t=1), and queried y(t=i)
            dYa = Ya - this.y0; % displacement of anchors
            dYk = zeros(size(Yq));
            for i = 1 : size(Yq,1)
                F = griddedInterpolant(Ya(i,:), dYa(i,:), 'makima', 'linear');
                dYk(i,:) = F(Yq(i,:));
%                 dYk(i,:) = interp1(Ya(i,:), dYa(i,:), Yq(i,:), 'makima', 'linear');
            end
            Yk = Yq - dYk;
            
            % Restore array size
            YqFull(ind,:) = Yk;
            Yk = YqFull;
            
            dYkFull = zeros(size(YqFull));
            dYkFull(ind,:) = dYk;
            dYk = dYkFull;
        end
        
        function Vk = CorrectVoltArray(this, t, y, V, varargin)
            % Apply reverse motion to an array of voltage data
            % 
            %   Vk = obj.CorrectVoltArray(t, y, V)
            %   Vk = obj.CorrectVoltArray(..., 'Highpass', false)
            % 
            % Inputs
            %   t               A m-element vector of sample times in second
            %   y               A n-element vector of channel depths for V
            %   V               A m-by-n matrix of voltage data
            %   'Highpass'      Whether or not to highpass filter the corrected voltage timeseries
            % Output
            %   Vk              A matrix of corrected voltage timeseries
            % 
            
            p = inputParser();
            p.addParameter('Highpass', false, @islogical);
            p.parse(varargin{:});
            isHp = p.Results.Highpass;
            
            % Sort channels by depth
            [y, indChan] = sort(y);
            V = V(:,indChan);
            
            % Find the correct depth
            t = t(:); % make t a column vector (not necessarily)
            y = y(:)'; % ensure y is a row vector
            Yk = this.CorrectDepths(t, y);
            
            % Interpolate voltage
            Vk = zeros(size(V));
            for i = 1 : size(V,1)
                Vk(i,:) = interp1(Yk(i,:), V(i,:), y, 'makima', 0); % somehow this mapping works
            end
            
            % Restore original channel order
            [~, yInd] = sort(indChan);
            Vk = Vk(:, yInd);
            
            if isHp
                % Prepare parameters for high-pass filtering
                %   Phy uses bandpass between 500Hz and 14.25kHz(.475*sample_rate) to visulize waveform;
                %   here uses 300Hz highpass to be consistent with Kilosort configuration
                D = designfilt('highpassiir', ...
                    'PassbandFrequency', 300, ...
                    'StopbandFrequency', 250, ...
                    'StopbandAttenuation', 60, ...
                    'PassbandRipple', 0.1, ...
                    'SampleRate', 1/diff(t([1 2])), ...
                    'DesignMethod', 'ellip'); % order of this filter is 8
                
                for i = 1 : size(Vk,2)
                    Vk(:,i) = filtfilt(D, Vk(:,i));
                end
            end
        end
        
        function outImec = CorrectBinary(this, imec, outDir, todoList)
            
            if ~exist(outDir, 'dir')
                mkdir(outDir);
            end
            
            todoList = lower(cellstr(todoList));
            isAP = ismember('ap', todoList);
            isLF = ismember('lf', todoList);
            fAP = {@(o,V,c,i) TransformRawData(o, V, c, i, 30e3)};
            fLF = {@(o,V,c,i) TransformRawData(o, V, c, i, 2500)};
            
            outImec = imec.saveTransformedDataset(outDir, ...
                'stem', imec.fileStem, ...
                'writeAP', isAP, ...
                'writeLF', isLF, ...
                'transformAP', fAP(isAP), ...
                'transformLF', fLF(isLF), ...
                'mappedChannelsOnly', false, ...
                'ChunkSize', 30e3*5);
            
            function Vk = TransformRawData(obj, V, chanInd, iSample, fs)
                t = double(iSample) / fs;
                y = obj.channelMap.ycoords; % for the 384 mapped channels
                V = double(V');
                Vk = V;
                isMapped = obj.channelMap.channelIdsMapped;
                Vk(:,isMapped) = this.CorrectVoltArray(t, y, V(:,isMapped));
                Vk = Vk';
            end
        end
        
        function traces = CorrectTraces(this, traces, outDir)
            % Apply backward interpolation to trace data
            
            if ischar(traces) || isstring(traces)
                traces = cellstr(traces);
            end
            
            if iscell(traces)
                % Apply correction to data tables
                for i = 1 : numel(traces)
                    % Check input type
                    tr = traces{i};
                    if ischar(tr)
                        load(tr, 'dataTb');
                    elseif istable(tr)
                        dataTb = tr;
                    else
                        error('The cell array must contain file paths or data table');
                    end
                    
                    % Apply correction
                    dataTb.y = this.CorrectDepths(dataTb.t, dataTb.y);
                    
                    % Save corrected data
                    if ischar(tr)
                        if ~exist(outDir, 'dir')
                            mkdir(outDir);
                        end
                        [~, fileName] = fileparts(tr);
                        savePath = fullfile(outDir, fileName);
                        save(savePath, 'dataTb');
                    end
                end
                
            elseif isa(traces, 'MTracerTrace')
                % Apply correction in each trace object
                for i = 1 : numel(traces)
                    dataTb = traces(i).dataTb;
                    dataTb.y = this.CorrectDepths(dataTb.t, dataTb.y);
                    traces(i).dataTb = dataTb;
                end
            end
        end
        
        function traces = RestoreTraces(this, traces, outDir)
            % Apply forward interpolation to trace data
            
            if ischar(traces) || isstring(traces)
                traces = cellstr(traces);
            end
            
            if iscell(traces)
                % Apply correction to data tables
                for i = 1 : numel(traces)
                    % Check input type
                    tr = traces{i};
                    if ischar(tr)
                        load(tr, 'dataTb');
                    elseif istable(tr)
                        dataTb = tr;
                    else
                        error('The cell array must contain file paths or data table');
                    end
                    
                    % Apply correction
                    dy = this.fwF(dataTb.t, dataTb.y);
                    dataTb.y = dataTb.y + dy;
                    
                    % Save corrected data
                    if ischar(tr)
                        if ~exist(outDir, 'dir')
                            mkdir(outDir);
                        end
                        [~, fileName] = fileparts(tr);
                        savePath = fullfile(outDir, fileName);
                        save(savePath, 'dataTb');
                    end
                end
                
            elseif isa(traces, 'MTracerTrace')
                % Apply correction in each trace object
                for i = 1 : numel(traces)
                    dataTb = traces(i).dataTb;
                    dataTb.y = this.CorrectDepths(dataTb.t, dataTb.y);
                    traces(i).dataTb = dataTb;
                end
            end
        end
    end
end
