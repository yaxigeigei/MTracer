function tbOut = SliceTimeSeries(this, tbIn, tWin, varargin)
% Return slices of time series data using time windows
%
%   tbOut = SliceTimeSeries(tbIn, tWin)
%   tbOut = SliceTimeSeries(tbIn, tWin, rowInd)
%   tbOut = SliceTimeSeries(tbIn, tWin, rowInd, colInd)
%   tbOut = SliceTimeSeries(..., 'Fill', 'none')
%   tbOut = SliceTimeSeries(..., 'ReferenceTime', [])
%
% Inputs
%   tbIn                A table of time series data or the name of a timeSeries table in the object.
%   tWin                1) An n-by-2 matrix where each row has the begin and end time of a window.
%                       When n equals 1, this window is applied to every rows. When n equals the height
%                       of the table or the number of selected rows, each window is applied to respective
%                       row. Inf values will be substituted by min or max time available, respectively.
%                       2) An n-element cell array where each element is a m-by-2 matrix of time windows.
%                       All m windows will be applied to a single corresponding row, resulting in m rows
%                       in tbOut. Options for the cell array length, n, are the same as described above.
%   rowInd              Integer or logical indices of rows to operate on and return. The default value is
%                       [] indicating all rows.
%   colInd              Integer or logical indices of columns to operate on and return. It can also be
%                       a cell array of column names of the input table. The default is [] indicating all
%                       columns.
%   'Fill'              When tWin exceeds data in an epoch, 'none' (default) does not fill anything in
%                       exceeded parts; 'bleed' will look for data from neighboring epochs up to the
%                       entire session.
%   'ReferenceTime'     Reference time to use with the 'bleed' option. The default value is [] and the
%                       method will look for reference times associated with the specified table.
% Output
%   tbOut               The output table with inquired time series data.

% Handle user inputs
p = inputParser();
p.addRequired('tbIn', @(x) ischar(x) || istable(x));
p.addRequired('tWin', @(x) isnumeric(x) || iscell(x));
p.addOptional('rowInd', [], @(x) isnumeric(x) || islogical(x));
p.addOptional('colInd', [], @(x) isnumeric(x) || islogical(x) || iscellstr(x));
p.addParameter('Fill', 'none', @(x) ismember(x, {'none', 'bleed', 'nan', 'nanbleed'}));
p.addParameter('ReferenceTime', [], @(x) isnumeric(x) && isvector(x));

p.parse(tbIn, tWin, varargin{:});
rowInd = p.Results.rowInd(:);
colInd = p.Results.colInd(:);
fillMethod = p.Results.Fill;
tRef = p.Results.ReferenceTime;

if ~istable(tbIn)
    assert(this.IValidateTableName(tbIn, true) == 3, '''%s'' is not a timeSeries table', tbIn);
    tRef = this.GetReferenceTime(tbIn);
    tbIn = this.GetTable(tbIn);
end

% Validate and standardize row and column indices
[rowInd, colInd] = this.IValidateTableIndexing(tbIn, rowInd, colInd);
if ~ismember(1, colInd)
    colInd = [1 colInd]; % make sure the time column is always included
end
tbIn = tbIn(:,colInd);

% Validate and standardize time windows
tWin = this.IValidateTimeWindows(tWin, height(tbIn), rowInd);

% Fill empty rows with NaN
for i = 1 : height(tbIn)
    if isempty(tbIn.time{i})
        tbIn.time{i} = NaN;
    end
end

% Replace Inf with ends of timestamp
tBound = cellfun(@(x) x([1 end])', tbIn.time, 'Uni', false);
for i = 1 : numel(tWin)
    bb = repmat(tBound{i}, size(tWin{i},1), 1);
    isWinInf = isinf(tWin{i});
    tWin{i}(isWinInf) = bb(isWinInf);
end

% Slicing
switch fillMethod
    case 'none'
        winCells = ISliceTsFillNone(tbIn, tWin, rowInd);
    case 'nan'
        warning('''Fill'', ''nan'' may not be supported in a future version');
        winCells = ISliceTsFillNaN(tbIn, tWin, rowInd);
    case 'bleed'
        assert(~isempty(tRef), 'Requires ''ReferenceTime'' to use ''%s''', fillMethod);
        winCells = ISliceTsFillBleed(tbIn, tWin, rowInd, tRef);
    case 'nanbleed'
        warning('''Fill'', ''nanbleed'' may not be supported in a future version');
        assert(~isempty(tRef), 'Requires ''ReferenceTime'' to use ''%s''', fillMethod);
        tPan = cellfun(@(x) [min(x(:)) max(x(:))], tWin, 'Uni', false);
        winCells = ISliceTsFillBleed(tbIn, tPan, rowInd, tRef);
        winCells = cat(1, winCells{:});
        tbIn = cell2table(winCells, 'VariableNames', tbIn.Properties.VariableNames);
        winCells = ISliceTsFillNaN(tbIn, tWin(rowInd), 1:numel(rowInd));
end
winCells = cat(1, winCells{:});

% Make table where all variable type is cell
tbOut = cell2table(cell(size(winCells)), 'VariableNames', tbIn.Properties.VariableNames);
tbOut{:,:} = winCells;
end


function winData = ISliceTsFillNone(tbIn, tWin, rowInd)
% Find time series for each epoch
winData = cell(size(tWin));
for i = rowInd
    t = tbIn.time{i};
    w = tWin{i};
    winData{i} = cell(size(w,1), width(tbIn));
    
    % Find time series for each window
    for j = 1 : size(w,1)
        isInWin = w(j,1) <= t & t < w(j,2);
        winData{i}(j,:) = cellfun(@(x) x(isInWin,:), tbIn{i,:}, 'Uni', false);
    end
end
end


function winData = ISliceTsFillNaN(tbIn, tWin, rowInd)
% Find time series for each epoch
winData = cell(size(tWin));
for i = rowInd
    t = tbIn.time{i};
    dtPre = diff(t(1:2));
    dtPost = diff(t(end-1:end));
    w = tWin{i};
    winData{i} = cell(size(w,1), width(tbIn));
    
    % Find time series for each window
    for j = 1 : size(w,1)
        isInWin = w(j,1) <= t & t < w(j,2);
        
        % Make time stamps in exceeded parts
        tPre = flip(t(1)-dtPre : -dtPre : w(j,1))';
        tPre(tPre >= w(j,2)) = [];
        
        tPost = (t(end)+dtPost : dtPost : w(j,2))';
        tPost(tPost < w(j,1)) = [];
        
        % Indexing and filling
        winData{i}{j,1} = [tPre; t(isInWin,:); tPost];
        winData{i}(2:end) = cellfun( ...
            @(x) [NaN(numel(tPre), size(x,2)); x(isInWin,:); NaN(numel(tPost), size(x,2))], ...
            tbIn{i,2:end}, 'Uni', false);
    end
end
end


function winData = ISliceTsFillBleed(tbIn, tWin, rowInd, tRef)
% Cache variables
isOldVer = verLessThan('matlab', '9.1');
tAbs = cellfun(@(x,r) x + r, tbIn.time, num2cell(tRef), 'Uni', false);
tAbsBegin = cellfun(@(x) x(1), tAbs);
tAbsEnd = cellfun(@(x) x(end), tAbs);

% Find time series for each epoch
winData = cell(size(tWin));
for i = rowInd
    wAbs = tWin{i} + tRef(i);
    winData{i} = cell(size(wAbs,1), width(tbIn));
    
    % Find time series for each window
    for j = 1 : size(wAbs,1)
        % Find relevant epochs
        if isOldVer
            dtWinBeforeEpoch = bsxfun(@minus, tAbsBegin, wAbs(j,:));
            dtWinAfterEpoch = bsxfun(@minus, wAbs(j,:), tAbsEnd);
        else
            dtWinBeforeEpoch = tAbsBegin - wAbs(j,:);
            dtWinAfterEpoch = wAbs(j,:) - tAbsEnd;
        end
        isWinBeforeEpoch = all(dtWinBeforeEpoch >= 0, 2);
        isWinAfterEpoch = all(dtWinAfterEpoch > 0, 2);
        isWinOverlapEpoch = ~(isWinBeforeEpoch | isWinAfterEpoch);
        srcRowInd = find(isWinOverlapEpoch);
        
        % Sort relevant epochs in temporal order
        [~, ord] = sort(tRef(srcRowInd));
        srcRowInd = srcRowInd(ord);
        
        % Get masks
        isInWin = cellfun(@(x) wAbs(j,1) <= x & x < wAbs(j,2), tAbs(srcRowInd), 'Uni', false);
        
        % Find times
        tCells = cellfun(@(x,y) x(y)-tRef(i), tAbs(srcRowInd), isInWin, 'Uni', false);
        winData{i}{j,1} = vertcat(tCells{:});
        
        % Find data
        for k = 2 : width(tbIn)
            dataCells = cellfun(@(x,y) x(y,:), tbIn{srcRowInd,k}, isInWin, 'Uni', false);
            winData{i}{j,k} = vertcat(dataCells{:});
        end
    end
end
end