classdef MMath
    % MMath is a collection of functions useful for doing math and manipulating data
    
    methods(Static)
        function edges = BinCenters2Edges(centers)
            % Convert bin centers to bin edges
            % 
            %   edges = BinCenters2Edges(centers)
            % 
            % Input
            %   centers             A numeric vector of n bin centers. n must be > 1.
            % Output
            %   edges               A numeric vector of n+1 bin edges
            % 
            assert(numel(centers) > 1, "Must provide more than one bin centers to derive the edges.");
            isCol = iscolumn(centers);
            if isCol
                centers = centers';
            end
            hw = diff(centers) / 2;
            edges = [centers(1)-hw(1), centers(2:end)-hw, centers(end)+hw(end)];
            if isCol
                edges = edges';
            end
        end
        
        function centers = BinEdges2Centers(edges)
            % Convert bin edges to bin centers
            % 
            %   centers = BinEdges2Centers(edges)
            % 
            % Input
            %   edges               A numeric vector of n bin edges. n must be > 1.
            % Output
            %   centers             A numeric vector of n-1 bin centers.
            % 
            assert(numel(edges) > 1, "Must provide more than one bin edges to derive the centers.");
            centers = edges(2:end) - diff(edges)/2;
        end
        
        function [ci, bootstat] = BootCI(nboot, bootfun, varargin)
            % Bootstraps confidence interval as the built-in bootci function but supports hierarchical bootstrap
            %
            %   [ci, bootstat] = BootCI(nboot, bootfun, D)
            %   [ci, bootstat] = BootCI(nboot, {bootfun, D})
            %   [ci, bootstat] = BootCI(..., 'alpha', 0.05)
            %   [ci, bootstat] = BootCI(..., 'Groups', [])
            %
            % Inputs:
            %   nboot           Number of resampling.
            %   bootfun         Function handle for computation.
            %   D               Data to be sent into bootfun.
            %   'alpha'         Significance level.
            %   'Groups'        A numeric array of grouping indices for hierarchical bootstrap. Each column 
            %                   contains indices for a level, from the highest (first column) to the lowest 
            %                   (last column). Do not (and no need to) include unique indices for individual 
            %                   samples as the last column. The default is empty and standard bootstrap 
            %                   will be used.
            % Output:
            %   ci              Confidence interval.
            %   bootstat        The distribution of computed values from bootstrap.
            
            % Process user inputs
            if ~ischar(varargin{1})
                D = varargin{1};
                varargin(1) = [];
            else
                D = bootfun{2};
                bootfun = bootfun{1};
            end
            
            p = inputParser();
            p.addParameter('alpha', 0.05, @(x) isscalar(x) && x<1 && x>0);
            p.addParameter('Groups', [], @isnumeric);
            p.addParameter('Options', statset('bootstrp'));
            p.parse(varargin{:});
            a = p.Results.alpha;
            G = p.Results.Groups;
            ops = p.Results.Options;
            
            nSample = size(D,1);
            if isempty(G)
                G = ones(nSample,1);
            end
            
            % Work out dimensions
            Dcolons = repmat({':'}, 1, ndims(D)-1);
            r = bootfun(D);
            bootstat = zeros(nboot, numel(r));
            
            % Compute bootstrap stats
            for n = 1 : nboot
                ind = MMath.HierBootSample(G);
                d = D(ind, Dcolons{:});
                r = bootfun(d);
                bootstat(n,:) = r(:);
            end
            
            % Compute CI
            bootstat = reshape(bootstat, [nboot size(r)]);
            ci = prctile(bootstat, 100*[a/2 1-a/2]', 1);
        end
        
        function [trainInd, testInd] = BootPartition(n, ptest, nboot)
            % Create crossvalidation partitions via bootstrap resampling
            %
            %   [trainInd, testInd] = MMath.BootPartition(n, ptest)
            %   [trainInd, testInd] = MMath.BootPartition(group, ptest)
            %   [trainInd, testInd] = MMath.BootPartition(..., nboot)
            % 
            % Inputs:
            %   n           a
            %   group       a
            %   ptest       a
            %   nboot       a
            % Outputs:
            %   trainInd    a
            %   testInd     a
            
            if nargin < 3
                nboot = 1;
            end
            if ~isscalar(n)
                n = numel(n);
            end
            
            indList = 1 : n;
            ntest = round(n * ptest);
            
            for i = nboot : -1 : 1
                testInd(:,i) = randsample(n, ntest);
                trainInd(:,i) = setdiff(indList, testInd(:,i));
            end
        end
        
        function [pval, sig] = EstimatePval(X, null, varargin)
            % Compute p-values from null distribution and the observed value(s)
            % 
            %   [pval, sig] = MMath.EstimatePval(X, null)
            %   [pval, sig] = MMath.EstimatePval(X, null, 'Tail', 'two')
            %   [pval, sig] = MMath.EstimatePval(X, null, 'AlphaList', [0.05 0.01 0.001])
            %   [pval, sig] = MMath.EstimatePval(X, null, 'Method', 'gumbel')
            % 
            % Inputs
            %   X               A scalar or vector of value(s) for test.
            %   null            A vector or a matrix of column vectors. Each vector contains values sampled from 
            %                   the null distribution for each value in X.
            %   'Tail'          'left', 'right', or 'two' (default) tailed test.
            %   'AlphaList'     A vector of significance levels. Default is [0.05 0.01 0.001].
            %   'Method'        'gumbel' (default): calculate p-values for maximum value of N samples drawn from Gaussian.
            %                   'empirical': calculate p-values using empirical quantiles.
            % Output
            %   pval            p-values in the same size as input x.
            %   sig             Level of significance given alpha value.
            % 
            % Credits
            %   Adapted from Montijn et al. eLife 2021, https://github.com/JorritMontijn/zetatest
            %   
            %   Version history:
            %   1.0 - March 11 2020
            %       Created by Jorrit Montijn
            %   
            %   Sources:
            %   Baglivo (2005)
            %   Elfving (1947), https://doi.org/10.1093/biomet/34.1-2.111
            %   Royston (1982), DOI: 10.2307/2347982
            %   https://stats.stackexchange.com/questions/394960/variance-of-normal-order-statistics
            %   https://stats.stackexchange.com/questions/9001/approximate-order-statistics-for-normal-random-variables
            %   https://en.wikipedia.org/wiki/Extreme_value_theory
            %   https://en.wikipedia.org/wiki/Gumbel_distribution
            %   
            
            p = inputParser;
            p.addParameter('Method', 'gumbel', @(x) any(strcmpi(x, {'empirical', 'gumbel'})));
            p.addParameter('AlphaList', [0.05 0.01 0.001], @(x) isnumeric(x));
            p.addParameter('Tail', 'two', @(x) any(strcmpi(x, {'left', 'right', 'two'})));
            p.parse(varargin{:});
            method = p.Results.Method;
            alphaList = p.Results.AlphaList;
            tailMode = p.Results.Tail;
            
            % Standardize inputs
            X = X(:)';
            if isvector(null)
                null = null(:);
                if numel(X) > 1
                    null = repmat(null, [1 numel(X)]);
                end
            end
            assert(numel(X)==size(null,2), "The number of elements in 'X' does not match the number of columns in 'null'.");
            
            if strcmpi(method, 'empirical')
                % Calculate statistical significance using empirical quantiles
                pval = nan(size(X));
                for i = 1 : numel(X)
                    x = X(i);
                    nullVec = unique(null(:,i));
                    if numel(nullVec) == 1 && x == nullVec
                        rk = NaN;
                    elseif x < min(nullVec) || isnan(x)
                        rk = 0;
                    elseif x > max(nullVec) || isinf(x)
                        rk = numel(nullVec);
                    else
                        rk = interp1(nullVec, 1:numel(nullVec), x);
                    end
                    pval(i) = rk / numel(nullVec);
                end
            else
                % Calculate statistical significance using Gumbel approximation
                nullMu = mean(null, 1);
                nullVar = var(null, 0, 1);
                
                % Define constants
                % define Euler-Mascheroni constant
                dblEulerMascheroni = 0.5772156649015328606065120900824; %vpa(eulergamma)
                
%                 % define apery's constant
%                 dblApery = 1.202056903159594285399738161511449990764986292;
                
                % Define Gumbel parameters from mean and variance
                % derive beta parameter from variance
                dblBeta = (sqrt(6).*sqrt(nullVar))./(pi);
                % dblSkewness = (12*sqrt(6)*dblApery)/(pi.^3);
                
                % Derive mode from mean, beta and E-M constant
                dblMode = nullMu - dblBeta.*dblEulerMascheroni;
                
                % Define Gumbel cdf
                fGumbelCDF = @(x) exp(-exp(-((x-dblMode)./dblBeta)));
                
                % Calculate output variables
                % calculate cum dens at X
                dblGumbelCDF = fGumbelCDF(X);
                
                % define p-value
                pval = dblGumbelCDF; % originally pval = (1-dblGumbelCDF);
%                 % transform to output z-score
%                 Z = -norminv(pval/2);
%                 
%                 % approximation for large X
%                 pval(isinf(Z)) = exp( (dblMode-X(isinf(Z)))./dblBeta );
            end
            
            % Get significance level
            alphaList = alphaList(:)'; % make row vector
            switch tailMode
                case 'two'
                    m = pval > 0.5;
                    pval(m) = 1 - pval(m);
                    pval = pval(:) * 2;
                    sig = pval < alphaList;
                case 'left'
                    sig = pval(:) < alphaList;
                case 'right'
                    pval = 1 - pval;
                    sig = pval(:) < alphaList;
                otherwise
                    error("'Tail' must be 'two', 'left', or 'right', but was '%s'.", tailMode);
            end
            sig = sum(sig, 2);
            if isrow(X)
                sig = sig';
            end
        end
        
        function [H, p] = BootTest2(s1, s2, varargin)
            % 
            
            p = inputParser;
            p.addParameter('Alpha', 0.05, @isscalar);
            p.addParameter('NBins', 100, @isscalar);
            p.parse(varargin{:});
            alpha = p.Results.Alpha;
            nBins = p.Results.NBins;
            
            if isempty(s1) || isempty(s2)
                H = NaN;
                p = NaN;
                return
            end
            
            % Upsampling
            s1 = s1(:);
            s2 = s2(:);
            % TBW
            
            % Fit GMM
            gm = fitgmdist([s1(:) s2(:)], 1);
            
            % Determine the range and resolution of quantization
            m = gm.mu;
            sd1 = std(s1);
            sd2 = std(s2);
            x = linspace(m(1)-5*sd1, m(1)+5*sd1, nBins);
            y = linspace(m(2)-5*sd2, m(2)+5*sd2, nBins);
            [X, Y] = meshgrid(x, y);
            
%             a = min([x y]);
%             b = max([x y]);
%             figure
%             plot(s1(:), s2(:), 'k.'); hold on
%             plot([a b]', [a b]', 'r', 'LineWidth', 2); 
%             xlabel('ctrl'); ylabel('opto'); axis xy equal tight
            
            % Evaluate probability density
            gmPDF = @(x,y) arrayfun(@(x0,y0) gm.pdf([x0 y0]), x, y);
            Z = gmPDF(X, Y);
            Z = Z ./ sum(Z(:));
            
            % Calculate propability beyond unity line
            isSide = X < Y;
            p = sum(Z(isSide(:)));
            if p > .5
                p = max(1-p, 0); % sometimes 1-p can produce 'negative zero'
            end
            p = p*2;
            H = p < alpha;
        end
        
        function val = Bound(val, valRange)
            %Bounds input values to specified range or allowed members
            %
            %   val = MMath.Bound(val, valRange)
            %
            % Inputs:
            %   val         Numeric scalar, vector or array
            %   valRange    Either (1) a tuple of boundaries (e.g. [ minVal, maxVal ])
            %               or (2) a numeric vector or array of all allowed values (number of elements must > 2) and
            %               each raw input value will be bound to the closest allowed value.
            % Output:
            %   val         Same size as the input val with bound values
            
            if numel(valRange) > 2
                for i = 1 : numel(val)
                    dists = abs(valRange(:) - val(i));
                    [ ~, idx ] = min(dists);
                    val(i) = valRange(idx);
                end
            else
                val = max(val, min(valRange));
                val = min(val, max(valRange));
            end
        end
        
        function mask = Bounds2Logical(bb, vecLength)
            % Convert a vector of indices to logical mask
            %
            %   mask = MMath.Bounds2Logical(bb, vecLength)
            %
            % Inputs:
            %   bb          A n-by-2 matrix, where each row bounds a range of indices (inclusive)
            %   vecLength   The length of the output vector
            % Output:
            %   mask        A logical vector where the bounded elements are true
            
            assert(size(bb,2)==2, 'The size of the second dimension of the bounds (bb) must equal two');
            
            mask = false(vecLength, 1);
            
            for i = 1 : size(bb,1)
                mask(bb(i,1) : bb(i,2)) = true;
            end
        end
        
        function B = CombineDims(A, dims)
            % Reshape array by combining specified dimensinons
            % 
            %   B = MMath.CombineDims(A, dims)
            %
            % Inputs
            %   A           Numeric array.
            %   dims        Dimensions to combine. Note that even [2 3] and [3 2] specify the same dimensions, 
            %               the order of elements in the combined dimension will be different. 
            % Output
            %   B           Resulting array where the first dimension is the combined. 
            
            % Move dimensions of interest to the front
            otherDims = setdiff(1:ndims(A), dims);
            A = permute(A, [dims otherDims]);
            
            % Collapse specified dimensions by reshaping array
            sz = [size(A) 1];
            dims = 1 : numel(dims);
            otherDims = numel(dims)+1 : numel(sz);
            B = reshape(A, [prod(sz(dims)) sz(otherDims)]);
        end
        
        function [condDist, margDist] = ConditionalDist(jointDist, givenWhich)
            %Compute conditional probability distribution from joint distribution
            %
            %   [condDist, margDist] = MMath.ConditionalDist(jointDist, givenWhich)
            %
            % Inputs
            %   jointDist       Multidimentional (>=2) array of joint probability
            %   givenWhich      Dimentions to condition on (currently has to be ndims(jointDist)-1 dimensions)
            % Outputs
            %   condDist        Array of conditional distribution
            %   margDist        Marginal distribution
            % 
            
            % Collapsing other dimensions
            numDims = ndims(jointDist);
            margWhich = setdiff(1:numDims, givenWhich);
            if isempty(margWhich)
                margDist = jointDist;
            else
                for i = 1 : length(margWhich)
                    margDist = sum(jointDist, margWhich(i), 'omitnan');
                end
            end
            
            % Normalization
            margSize = ones(1, numDims);
            margSize(1:ndims(margDist)) = size(margDist);
            condDist = jointDist ./ repmat(margDist, size(jointDist)./margSize);
            
            % Set NaN to zero
            condDist(isnan(condDist)) = 0;
        end
        
        function r = Crossval(func, X, Y, trainInd, testInd)
            % 
            
            for i = size(trainInd,2) : -1 : 1
                xTrain = X(trainInd(:,i), :);
                yTrain = Y(trainInd(:,i));
                xTest = X(testInd(:,i), :);
                yTest = Y(testInd(:,i));
                yTestHat = func(xTest, xTrain, yTrain);
                r(i) = mean(yTestHat == yTest);
            end
        end
        
        function [x, D] = Decimate(x, r, n, D)
            %Similar to the MATLAB decimate function but with the following modifications
            % 1) Accepts data types other than double
            % 2) If x is an matrix, treats each column as a time series
            % 3) To be consistent with MATLAB downsample function, it uses the (n+1)-th 
            % 	 value in each bin of r samples
            % 4) By default, it uses 10th-order elliptic (IIR) lowpass filter rather than 
            %    cheby1 or FIRs. 
            %    - fpass is fs/r/2 (i.e. Nyquist freq. of the output), fstop is ~1.2*fpass
            %    - passband has ~1/1000th uneveness (ripple) in amplitude (or 0.01dB)
            %    - stopband is attenuated >10000 times in amplitude (or 80dB)
            % 5) User can provide custom digitalFilter object
            %   
            %   [x, D] = MMath.Decimate(x, r)
            %   [x, D] = MMath.Decimate(x, r, n)
            %   [x, D] = MMath.Decimate(x, r, n, D)
            
            if nargin < 4
                D = designfilt('lowpassiir', ...
                    'FilterOrder', 10, ...
                    'PassbandFrequency', 1/r, ...
                    'PassbandRipple', .01, ...
                    'StopbandAttenuation', 80);
            end
            if nargin < 3
                n = 0;
            end
            
            isRow = isrow(x);
            if isRow
                x = x';
            end
            
            % Filter one column at a time to save memory and conserve data type
            for i = 1 : size(x,2)
                x(:,i) = filtfilt(D, double(x(:,i)));
            end
            
            x = downsample(x, r, n);
            
            if isRow
                x = x';
            end
        end
        
        function H = Entropy(probDist)
            %Calculates the Shannon's entropy of a given probability distribution
            %
            %   H = MMath.Entropy(probDist)
            %
            % Inputs:
            %   probDist    A numeric array of probability values
            % Output:
            %   H           Entropy in bits
            
            probDist = probDist(:);
            probDist(probDist < 0) = 0;
            probDist = probDist / sum(probDist, 'omitnan');
            H = -sum(probDist .* log2(probDist), 'omitnan');
        end
        
        function expect = Expectation(distMat, val, idxDim)
            %Calculate expected values in the dimension of interest given other variables
            %
            %   expect = MMath.Expectation(distMat, val, idxDim)
            %
            % Inputs:
            %   distMat     Joint probability distribution
            %   val         Values of each probabilities in the dimension of interest
            %   idxDim      Index of the dimension of interest
            % Output:
            %   expect      Expected values
            
            if nargin < 3
                idxDim = 1;
            end
            
            sizeDist = size(distMat);
            reshapeVect = ones(1, length(sizeDist));
            reshapeVect(idxDim) = length(val);
            val = reshape(val, reshapeVect);
            
            sizeDist(idxDim) = 1;
            expect = nansum(distMat .* repmat(val, sizeDist), idxDim);
        end
        
        function ind = HierBootSample(G)
            % Hierarchical bootstrap sampling
            %
            %   ind = HierBootSample(G)
            %
            % Input
            %   G       A numeric array of grouping indices for hierarchical bootstrap. Each column 
            %           contains indices for a level, from the highest (first column) to the lowest 
            %           (last column). Do not (and no need to) include unique indices for individual 
            %           samples as the last column.
            % Output
            %   ind     Bootstrap sample indices.
            
            % Add unique sample indices to the end of grouping indices
            [nSample, nLevel] = size(G);
            G = [G, (1:nSample)'];
            
            % Sampling
            GG = sampleHier(G);
            
            % Denest cell array
            for n = 1 : nLevel
                GG = cat(1, GG{:});
            end
            ind = GG;
            
            function GG = sampleHier(G)
                % Split the rest of grouping indices by the top level (1st column) indices
                G(:,1) = findgroups(G(:,1));
                GG = splitapply(@(x) {x}, G(:,2:end), G(:,1));
                
                % Random sampling with replacement
                GG = randsample(GG, numel(GG), true);
                
                if size(GG{1},2) > 1
                    % Recursively sample a lower level
                    for i = 1 : numel(GG)
                        GG{i} = sampleHier(GG{i});
                    end
                else
                    % Sample the last level
                    for i = 1 : numel(GG)
                        GG{i} = randsample(GG{i}, numel(GG{i}), true);
                    end
                end
            end
            
%             n1 = 4;
%             n2 = [3 1 4 4];
%             n3 = randi(10, sum(n2), 1);
%             
%             g1 = repelem((1:n1)', n2');
%             g1 = repelem(g1, n3');
%             g2 = repelem((1:sum(n2))', n3');
%             G = [g1 g2];
%             
%             ind = MMath.HierBootSample(G)
        end
        
        function roiInd = Ind2Roi(ind, winRoi, valRange)
            %Converts each index to indices for indexing a corresponding region of interest (ROI)
            %e.g. indices [ 23; 56 ] with a window of [ -1 0 1 2 ] => ROIs of [ 22 23 24 25; 55 56 57 58 ]
            %
            %   roiInd = MMath.Ind2Roi(ind, winRoi, valRange)
            %
            % Inputs:
            %   ind         An index or a vector of indices
            %   winRoi      A vector of relative indices (not a boundary tuple) indicating the region of interest
            %   valRange    Usually (1) a tuple of indexing boundaries (e.g. [ minIdx, maxIdx ])
            %               or (2) a numeric vector or array of all allowed index values (number of elements must > 2)
            %               and each index will be bound to the closest allowed value.
            % Output:
            %   roiInd      An array where each row is a set of ROI indices for one input index
            %               ROIs causing edge issue are removed from the reuslt.
            
            if islogical(ind) || ~isempty(find(ind == 0, 1))
                if nargin < 3
                    valRange = [ 1 length(ind) ];
                end
                ind = find(ind);
            end
            
            roiInd = zeros(length(ind), length(winRoi));
            for i = 1 : length(winRoi)
                roiInd(:,i) = ind + winRoi(i);
            end
            
            roiInd = MMath.Bound(roiInd, valRange);
            
            validRoiMask = true(size(roiInd,1), 1);
            roiLength = size(roiInd, 2);
            for i = 1 : size(roiInd, 1)
                if unique(roiInd(i,:)) < roiLength
                    validRoiMask(i) = false;
                end
            end
            roiInd = roiInd(validRoiMask, :);
        end
        
        function mask = Ind2Logical(ind, vecLength)
            % Convert a vector of indices to logical mask
            %
            %   mask = MMath.Ind2Logical(ind, vecLength)
            %
            % Inputs:
            %   ind         A vector of indices
            %   vecLength   The length of the output vector
            % Output:
            %   mask        A logical vector where elements at 'ind' are true
            
            if isrow(ind)
                mask = false(1, vecLength);
            else
                mask = false(vecLength, 1);
            end
            
            mask(ind) = true;
        end
        
        function [ jointDist, axisLabels, binCoor ] = JointDist(randVars, numBins, axisLabels)
            % Calculates the joint probability distribution of multiple random variables
            %
            %   [ jointDist, axisLabels, binInd ] = MMath.JointDist(randVars, numBins)
            %
            % Inputs:
            %   randVars	   	Random variables in an array of column vectors. Each row is an observation.
            %   numBins         The number(s) of bins for input random variables. Use NaN to indicate categorical
            %                   variable. You may assign different bin numbers to different variables in a
            %                   vector otherwise they all use the same bin number.
            %   axisLabels      For categorical variables, you can provide predefined categories (e.g. to account
            %                   for unobserved values). It expects a 1-by-n cell array where n is the number of
            %                   variables in randVars. Each element in this cell array is a vector of numeric
            %                   categories along that dimension.
            % Output:
            %   jointDist       N-dimensional matrix of joint probability. The order of dimensions is the same
            %                   as the order of input variables.
            %   axisLabels      Labels of each axis in cell array.
            %   binCoor         Each row is a coordinate in the joint distribution sapce for the corresponding
            %                   observation (row) in randVars.
            
            % Handles user inputs
            if nargin < 3
                axisLabels = cell(1, size(randVars,2));
            end
            
            if length(numBins) ~= size(randVars,2)
                numBins = repmat(numBins(1), 1, size(randVars,2));
            end
            
            % Bins variables
            binCoor = zeros(size(randVars));
            for i = size(randVars,2) : -1 : 1
                if isnan(numBins(i))
                    % Catagorical
                    axisLabels{i} = union(axisLabels{i}, unique(randVars(:,i)));
                    binCoor(:,i) = arrayfun(@(x) find(axisLabels{i} == x), randVars(:,i));
                    numBins(i) = length(axisLabels{i});
                else
                    % Continuous
                    edges = linspace(nanmin(randVars(:,i)), nanmax(randVars(:,i)), numBins(i)+1);
                    axisLabels{i} = mean([ edges(1:end-1); edges(2:end) ])';
                    edges(end) = edges(end) + 1;
                    [ ~, binCoor(:,i) ] = histc(randVars(:,i), edges);
                end
            end
            
            % Multidimensional tally
            dimInputStr = 'ndgrid(';
            for i = 1 : length(numBins)
                dimInputStr = [ dimInputStr, '1:numBins(', num2str(i), '), ' ];
            end
            dimInputStr = [ dimInputStr(1:end-2), ')' ];
            
            dimOutput = cell(1, length(numBins));
            [ dimOutput{:} ] = eval(dimInputStr);
            dimOutput = cell2mat(cellfun(@(x) x(:), dimOutput, 'UniformOutput', false));
            
            jointDist = zeros(numBins(:)');
            for i = size(dimOutput,1) : -1 : 1
                subBinCoor = binCoor;
                for j = 1 : length(numBins)
                    hitVect = subBinCoor(:,j) == dimOutput(i,j);
                    subBinCoor = subBinCoor(hitVect, :);
                end
                jointDist(i) = sum(hitVect);
            end
            
            jointDist = jointDist / size(randVars,1);
        end
        
        function [H, pValue, KSstatistic] = KStest2CDF(ecdf1, ecdf2, varargin)
            %KSTEST2 Two-sample Kolmogorov-Smirnov goodness-of-fit hypothesis test.
            
            if nargin > 2
                [varargin{:}] = convertStringsToChars(varargin{:});
            end
            
            if nargin < 2
                error(message('stats:kstest2:TooFewInputs'));
            end
            
            % Parse optional inputs
            alpha = []; tail = [];
            if nargin >=3
                if isnumeric(varargin{1})
                    % Old syntax
                    alpha = varargin{1};
                    if nargin == 4
                        tail = varargin{2};
                    end
                else
                    % New syntax
                    params = {'alpha', 'tail'};
                    dflts =  { []     , []};
                    [alpha, tail] =...
                        internal.stats.parseArgs(params, dflts, varargin{:});
                end
            end
            
            % Ensure each sample is a VECTOR.
            ecdf1  =  ecdf1(:);
            ecdf2  =  ecdf2(:);
            if isempty(ecdf1)
                error(message('stats:kstest2:NotEnoughData', 'X1'));
            end
            if isempty(ecdf2)
                error(message('stats:kstest2:NotEnoughData', 'X2'));
            end
            
            % Ensure the significance level, ALPHA, is a scalar
            % between 0 and 1 and set default if necessary.
            
            if ~isempty(alpha)
                if ~isscalar(alpha) || (alpha <= 0 || alpha >= 1)
                    error(message('stats:kstest2:BadAlpha'));
                end
            else
                alpha  =  0.05;
            end
            
            % Ensure the type-of-test indicator, TAIL, is a string or scalar integer
            % from the allowable set, and set default if necessary.
            if ~isempty(tail)
                if ischar(tail)
                    try
                        [~,tail] = internal.stats.getParamVal(tail, ...
                            {'smaller','unequal','larger'},'Tail');
                    catch
                        error(message('stats:kstest2:BadTail'));
                    end
                    tail = tail - 2;
                elseif ~isscalar(tail) || ~((tail==-1) || (tail==0) || (tail==1))
                    error(message('stats:kstest2:BadTail'));
                end
            else
                tail  =  0;
            end
            
            % Compute the test statistic of interest.
            switch tail
                case  0      %  2-sided test: T = max|F1(x) - F2(x)|.
                    deltaCDF  =  abs(ecdf1 - ecdf2);
                case -1      %  1-sided test: T = max[F2(x) - F1(x)].
                    deltaCDF  =  ecdf2 - ecdf1;
                case  1      %  1-sided test: T = max[F1(x) - F2(x)].
                    deltaCDF  =  ecdf1 - ecdf2;
            end
            KSstatistic   =  max(deltaCDF);
            
            % Compute the asymptotic P-value approximation and accept or
            % reject the null hypothesis on the basis of the P-value.
            n1     =  length(ecdf1);
            n2     =  length(ecdf2);
            n      =  n1 * n2 /(n1 + n2);
            lambda =  max((sqrt(n) + 0.12 + 0.11/sqrt(n)) * KSstatistic , 0);
            if tail ~= 0        % 1-sided test.
                pValue  =  exp(-2 * lambda * lambda);
            else                % 2-sided test (default).
                %  Use the asymptotic Q-function to approximate the 2-sided P-value.
                j       =  (1:101)';
                pValue  =  2 * sum((-1).^(j-1).*exp(-2*lambda*lambda*j.^2));
                pValue  =  min(max(pValue, 0), 1);
            end
            
            H  =  (alpha >= pValue);
        end
        
        function boundaryInd = Logical2Bounds(vect)
            %Converts logical(boolean/binary) vector to tuples indicating the bondaries of 1 regions
            %e.g. [ 0 0 0 1 1 1 0 0 1 ] => bondaries [ 4 6; 9 9 ]
            %
            %   boundaryInd = MMath.Logical2Bounds(vect)
            %
            % Inputs:
            %   vect            A logical(boolean/binary) vector
            % Output:
            %   boundaryInd     Tuples (each row) indicating the bondaries of 1 regions
            
            if isempty(vect)
                boundaryInd = zeros(0,2);
                return;
            end
            
            vect = logical(vect(:));
            vect = [ zeros(1,size(vect,2)); vect; zeros(1,size(vect,2)) ];
            dVect = diff(vect);
            boundaryInd = [ find(dVect == 1), find(dVect == -1) - 1 ];
        end
        
        function varargout = MeanStats(A, dim, varargin)
            % Compute means, SDs, SEMs and bootstrap CIs from samples
            %
            %   [m, sd, se, ci] = MMath.MeanStats(A)
            %   [m, sd, se, ci] = MMath.MeanStats(A, dim)
            %   [m, sd, se, ci] = MMath.MeanStats(A, dim, ..., 'IsOutlierArgs', {})
            %   [m, sd, se, ci] = MMath.MeanStats(A, dim, ..., 'NBoot', 1000)
            %   [m, sd, se, ci] = MMath.MeanStats(A, dim, ..., 'Alpha', 0.05)
            %   [m, sd, se, ci] = MMath.MeanStats(A, dim, ..., 'Options', statset)
            %
            % Inputs
            %   A                   Numeric array of samples.
            %   dim                 Dimension to operate along.
            %   'isoutlierArgs'     A cell array for one or more arguments of the isoutlier 
            %                       function to remove outliers in A.
            %   'NBoot'             The number of bootstrap sampling.
            %   'Alpha'             Significance level. Default is 0.05.
            % Outputs
            %   All otuputs has the same dimensionality as A with statistical values in the dim 
            %   dimension. 
            %   m                   Mean values.
            %   sd                  Standard deviations.
            %   se                  Standard error of the mean.
            %   ci                  95% bootstrap confidence intervals of the mean.
            
            p = inputParser;
            p.addParameter('IsOutlierArgs', {});
            p.addParameter('NBoot', 1000);
            p.addParameter('Alpha', 0.05);
            p.addParameter('Options', statset);
            p.parse(varargin{:});
            isoutlierArgs = p.Results.IsOutlierArgs;
            nboot = p.Results.NBoot;
            alphaVal = p.Results.Alpha;
            ops = p.Results.Options;
            
            if nargin < 2
                % Default using the first non-singleton dimension
                if numel(A) < 2
                    dim = 1;
                else
                    dim = find(size(A) > 1, 1);
                end
            end
            
            if isempty(A)
                % Create NaN array with matching size
                sz = size(A);
                sz(sz==0) = 1;
                A = NaN(sz);
            end
            
            if ~iscell(isoutlierArgs)
                isoutlierArgs = {isoutlierArgs};
            end
            if numel(isoutlierArgs) > 0
                A(isoutlier(A, isoutlierArgs{:}, dim)) = NaN;
            end
            
            m = mean(A, dim, 'omitnan');
            sd = std(A, 0, dim, 'omitnan');
            se = sd ./ sqrt(size(A,dim));
            varargout{1} = m;
            varargout{2} = sd;
            varargout{3} = se;
            
            if nargout > 3
                if size(A,dim) < 2
                    ci = cat(dim, m, m);
                else
                    % bootci can only sample along the first dimension, thus permuting A
                    dimOrder = [dim setdiff(1:ndims(A), dim)];
                    A = permute(A, dimOrder);
                    ci = bootci(nboot, {@(x) mean(x, 'omitnan'), A}, 'alpha', alphaVal, 'Options', ops);
                    
                    % Restore original dimension order
                    ci = permute(ci, [1 3:ndims(ci) 2]); % squeeze out the second (mean value) dimension
                    ci = ipermute(ci, dimOrder);
                end
                varargout{4} = ci;
            end
        end
        
        function [m, qt, ad] = MedianStats(A, dim)
            % Compute medians, 1st and 3rd quartiles and absolute deviations from samples
            %
            %   [m, qt, ad] = MMath.MedianStats(A)
            %   [m, qt, ad] = MMath.MedianStats(A, dim)
            %
            % Inputs
            %   A           Numeric array of samples.
            %   dim         Dimension to operate along. Default is the first non-singleton dimension. 
            % Outputs
            %   m           Median values.
            %   qt          1st and 3rd quartiles.
            %   ad          Median absolute deviation.
            
            if isempty(A)
                A = NaN;
            end
            
            if nargin < 2
                if isscalar(A)
                    dim = 1;
                else
                    dim = find(size(A) > 1, 1);
                end
            end
            
            m = median(A, dim, 'omitnan');
            qt = prctile(A, [25 75], dim);
            ad = mad(A, 1, dim);
        end
        
        function I = MutualInfo(jointDist)
            %Calculates the mutual information of two random variables from their joint probability distribution
            %
            %   I = MMath.MutualInfo(jointDist)
            %
            % Inputs:
            %   jointDist   The joint probability distribution matrix of two random variables
            % Output:
            %   I           Mutual information in bits
            
            marginDist{1} = nansum(jointDist, 2);
            marginDist{2} = nansum(jointDist, 1);
            [ marginDist{2}, marginDist{1} ] = meshgrid(marginDist{2}, marginDist{1});
            I = jointDist .* log2(jointDist ./ (marginDist{1}.*marginDist{2}));
            I = nansum(I(:));
        end
        
        function [N, c, k] = Normalize(A, normType, normParam)
            % Normalize data in the columns of an array. The 'soft' version of each method scales 
            % smaller signals less than larger signals.
            % 
            %   [N, c, k] = MMath.Normalize(A)
            %   [N, c, k] = MMath.Normalize(A, normType)
            %   [N, c, k] = MMath.Normalize(A, normType, normParam)
            % 
            % Inputs
            %   A           An m-by-n numeric array where each column will be normalized.
            %   normType    Standard normalization methods:
            %                 'zscore' (default): zscore each timeseries.
            %                 'max': normalize each timeseries by maximum or peak value.
            %                 'minmax': normalize each timeseries by maximum or peak value.
            %                 'none': no normalization.
            %               Soft normalization (e.g. Russo et al. 2018):
            %                 If appending a method name with 'soft', the normParam value will be 
            %                 added to the normalization target. For example, 'maxsoft' normalizes 
            %                 timeseries to max(A)+normParam rather than max(A). 'zscoresoft' 
            %                 treats normParam as a value added to the amplitude of a sinusoidal 
            %                 timeseries, and converts it to a value added to standard deviation.
            %   normParam   See the soft normalization above.
            % Outputs
            %   N           An m-by-n array where each columns is normalized by N = (A-c)./k
            %   c           A 1-by-n vector of offsets applied to A.
            %   k           A 1-by-n vector of scaling factors applied after offset.
            % 
            
            if nargin < 2
                normType = 'zscore';
            end
            normType = lower(char(normType));
            
            k = ones(1, size(A,2));
            c = zeros(1, size(A,2));
            switch normType
                case 'max'
                    k = max(A, [], 1, 'omitnan');
                case 'maxsoft'
                    k = max(A+normParam, [], 1, 'omitnan');
                case 'minmax'
                    c = min(A, [], 1, 'omitnan');
                    k = max(A-c, [], 1, 'omitnan');
                case 'minmaxsoft'
                    c = min(A, [], 1, 'omitnan');
                    k = max(A-c+normParam, [], 1, 'omitnan');
                case 'zscore'
                    c = mean(A, 1, 'omitnan');
                    k = std(A-c, 0, 1, 'omitnan');
                case 'zscoresoft'
                    c = mean(A, 1, 'omitnan');
                    k = std(A-c, 0, 1, 'omitnan') + normParam*0.707/2; % 0.707/2 converts max to equivalent std
                case 'none'
                    % do nothing
                otherwise
                    error("'%s' is not a valid normalization option.", normType);
            end
            k(k==0) = eps;
            N = (A-c)./k;
        end
        
        function [r2, r2adj] = RSquared(X, y, b, c)
            % Compute R-squared and adjusted R-squared
            %
            %   [r2, r2adj] = MMath.RSquared(X, y, b)
            %   [r2, r2adj] = MMath.RSquared(X, y, b, c)
            %
            % Inputs
            %   X       An array of samples where columns are predictors and rows are observations.
            %   y       A vector of responses.
            %   b       Regression coefficients.
            %   c       Regression intercept.
            % Output
            %   r2      R-Squared (coefficient of determination).
            %   r2adj   Adjusted R-Squared.
            if nargin < 4
                c = zeros(1,size(X,2));
                p = size(b,1);
            else
                p = size(b,1) + 1;
            end
            if iscolumn(c)
                c = c';
            end
            n = size(X,1);
            yhat = X*b + c;
            rhat = y - yhat;
            rtotal = y - mean(y, 1, 'omitnan');
            SSE = sum(rhat.^2, 1, 'omitnan'); % error sum of squares
            TSS = sum(rtotal.^2, 1, 'omitnan'); % total sum of squares
            r2 = 1 - SSE./TSS;
            r2adj = 1 - (n-1)/(n-p) * SSE./TSS;
        end
        
        function rsM = RepShiftMat(M, K, varargin)
            % Shift element positions in an array along the first dimension
            % 
            %   rsM = RepShiftMat(M, K)
            %   rsM = RepShiftMat(M, K, fillVal)
            %   rsM = RepShiftMat(..., 'CellOutput', false)
            % 
            % Inputs
            %   M               Matrix to replicate and shift.
            %   K               A vector of the indices to shift by. The size of this vector, size(K), 
            %                   determines how the replicated matrices will be arranged. 
            %   fillVal         What value to fill when elements are shifted away.
            %   'CellOutput'    Whether to output individual replicated and shifted matrices in cell 
            %                   array. Default is false where cell array will be denested into a single 
            %                   numeric array.
            % Output
            %   rsM             Replicated and shifted matrix.
            % 
            
            p = inputParser();
            p.addOptional('fillVal', NaN, @(x) isnumeric(x) && isscalar(x));
            p.addParameter('CellOutput', false, @islogical);
            p.parse(varargin{:});
            fillVal = p.Results.fillVal;
            isCell = p.Results.CellOutput;
            
            nRow = size(M, 1);
            ind = (1:nRow)' + K(:)';
            isOOB = ind < 1 | ind > nRow; % find out-of-bond indices
            ind(isOOB) = 1; % make indices legal
            
            rsM = cell(size(K));
            for k = 1 : numel(K)
                rsM{k} = M(ind(:,k), :);
                rsM{k}(isOOB(:,k), :) = fillVal;
            end
            
            if ~isCell
                rsM = cell2mat(rsM);
            end
        end
        
        function [aSorted, I] = SortLike(a, b, isVerbose)
            % Sort elements in 'a' like how they are ordered in 'b'
            %
            %   [aSorted, I] = MMath.SortLike(a, b)
            %   [aSorted, I] = MMath.SortLike(a, b, isVerbose)
            %
            % Inputs
            %   a           A vector to sort on.
            %   b           A vector with elements in target order to sort to.
            %   isVerbose   Whether or not to show warnings. Default is true.
            % Outputs
            %   aSorted     The sorted vector a.
            %   I           A vector of indices where aSorted = a(I).
            % 
            
            if nargin < 3
                isVerbose = true;
            end
            
            nNotUni = numel(a) - numel(unique(a));
            if nNotUni && isVerbose
                warning("%i elements in 'a' are not unique. Only the first element of each unique value will be used.", nNotUni);
            end

            I = cell(size(b));
            for i = 1 : numel(b)
                I{i} = find(a==b(i), 1);
            end
            
            isAbsent = cellfun(@isempty, I);
            if any(isAbsent) && isVerbose
                warning("%i members in 'b' cannot be found in 'a', and will be ignored.", sum(isAbsent));
            end
            
            I = cat(1, I{~isAbsent});
            aSorted = a(I);
        end
        
        function B = SqueezeDims(A, dims)
            % Squeeze only the specified dimensions
            %
            %   B = MMath.SqueezeDims(A, dims)
            %
            % Inputs
            %   A       Numeric array.
            %   dims    Dimensions to squeeze.
            % Output
            %   B       Resulting array.
            
            sz = size(A);
            assert(all(sz(dims) == 1), 'Specified dimensions must all have size of 1');
            dimKeep = setdiff(1:ndims(A), dims);
            B = permute(A, [dimKeep dims]);
        end
        
        function [se, sd] = StandardError(A, dim)
            % Calculates the standard error of the population mean
            %
            %   [se, sd] = MMath.StandardError(A)
            %   [se, sd] = MMath.StandardError(A, dim)
            %
            % Inputs
            %   A           Numeric array of samples.
            %   dim         Dimension to operate along.
            % Outputs
            %   se          Standard error of the population mean.
            %   sd          Standard deviation of the population mean.
            
            if isempty(A)
                A = NaN;
            end
            
            if nargin < 2
                if isscalar(A)
                    dim = 1;
                else
                    dim = find(size(A) > 1, 1);
                end
            end
            
            sd = std(A, 0, dim, 'omitnan');
            se = sd ./ sqrt(size(A,dim));
        end
        
        function ranges = ValueBounds(vect, varargin)
            % Find boundary indices of the same values in a vector
            %
            %   boundaryInd = MMath.ValueBounds(vect)
            %   boundaryInd = MMath.ValueBounds(vect, categories)
            %   boundaryInd = MMath.ValueBounds(vect, ..., 'UniformOutput', true)
            %
            % Inputs
            %   vect            A vector.
            %   categories      A set of values about which to return boundary indices. Useful when 
            %                   one needs to return boundaries for a subset or a superset of values 
            %                   present in vect. The default is the unique values in vect. 
            %   'UniformOutput' Whether to combine boundary indices of each category into a 
            %                   matrix or to store them in a cell array. Default true, to combine.
            % Output
            %   boundaryInd     Tuples (each row) indicating the bondaries of 1 regions
            
            p = inputParser;
            p.addOptional('categories', []);
            p.addParameter('UniformOutput', true, @islogical);
            p.parse(varargin{:});
            cats = p.Results.categories;
            isUni = p.Results.UniformOutput;
            
            if isempty(cats)
                cats = unique(vect, 'stable');
            end
            ranges = cell(numel(cats),1);
            for k = 1 : numel(cats)
                ranges{k} = MMath.Logical2Bounds(vect == cats(k));
            end
            if isUni
                assert(all(cellfun(@(x) size(x,1)==1, ranges)), ...
                    'For ''UniformOutput'' set to true, every category must produce one pair of bounds');
                ranges = cell2mat(ranges);
            end
        end
        
        function VE = VarExplained(X, B, varType)
            % Percent variance explained (VE)
            % 
            %   VE = MMath.VarExplained(X, B)
            %   VE = MMath.VarExplained(X, B, varType)
            % 
            % Inputs 
            %   X           A matrix whose columns are variables and rows are observations.
            %   B           A matrix of basis vectors in column.
            %   varType     'each'(default): compute VE by each basis in B independently.
            %               'span': compute VE by a set of orthonormal bases that span B.
            %               'unique': compute variance uniquely explained by each basis in B.
            % Output
            %   VE          A row vector of percent variance explained on each basis, or a
            %               scalar if varType is 'span'.
            % 
            
            if nargin < 3 || size(B,2) == 1
                varType = 'each';
            end
            varType = lower(char(varType));
            
            X = rmmissing(X);
            
            switch varType
                case 'each'
                    % Normalize bases
                    L = sqrt(sum(B.^2, 1));
                    L(L == 0) = eps;
                    B = B ./ L;
                    
                    % Compute variances
                    varSub = var(X*B);
                    
                case 'span'
                    % Compute total variance in the full span of B
                    varSub = sum(var(X*orth(B), 0, 1));
                    
                case 'unique'
                    % Compute total variance in the full span of B
                    varB = sum(var(X*orth(B), 0, 1));
                    
                    % Compute total variance in spans of B minus single axes
                    for i = size(B,2) : -1 : 1
                        B_ = B(:,[1:i-1 i+1:end]);
                        varB_(i) = sum(var(X*orth(B_), 0, 1));
                    end
                    
                    % Subtract to get the unique variances
                    varSub = varB - varB_;
                    
                otherwise
                    error("'%s' is not a valid varType.", varType);
            end
            
            varTotal = trace(cov(X));
            VE = varSub/varTotal*100;
        end
        
        function C = VecCosine(A, B)
            % Compute pairwise cosine of vectors
            %
            %   C = MMath.VecCosine(A)
            %   C = MMath.VecCosine(A, B)
            %
            % Inputs
            %   A,B     Matrices of column vectors. If only A is provided, the function computes
            %           pairwise cosine among columns of A, otherwise between columns of A and B.
            % Output
            %   C       A matrix of pairwise cosine values. If A and B are both provided, rows 
            %           correspond to vectors in A and columns correspond to vectors in B.
            
            L = sqrt(sum(A.^2));
            L(L == 0) = eps;
            A = A ./ L;
            
            if nargin < 2
                B = A;
            else
                L = sqrt(sum(B.^2));
                L(L == 0) = eps;
                B = B ./ L;
            end
            
            C = A'*B;
        end
        
    end
    
end

