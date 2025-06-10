function tbOut = SliceEventTimes(this, tbIn, tWin, varargin)
% Return slices of event time data using time windows
%
%   tbOut = SliceEventTimes(tbIn, tWin)
%   tbOut = SliceEventTimes(tbIn, tWin, rowInd)
%   tbOut = SliceEventTimes(tbIn, tWin, rowInd, colInd)
%   tbOut = SliceEventTimes(..., 'Fill', 'none')
%   tbOut = SliceEventTimes(..., 'ReferenceTime', [])
%
% Inputs
%   tbIn                A table of event times data or the name of a eventTimes table in the object.
%   tWin                1) An n-by-2 matrix where each row has the begin and end time of a window.
%                       When n equals 1, this window is applied to every rows. When n equals the height
%                       of the table or the number of selected rows, each window is applied to respective
%                       row.
%                       2) An n-element cell array where each element is a m-by-2 matrix of time windows.
%                       n must equal to the height of the table or the number of selected rows. All m
%                       windows in a m-by-2 matrix are applied to a corresponding row.
%   rowInd              Integer or logical indices of rows to operate on and return. The default value is
%                       [] indicating all rows.
%   colInd              Integer or logical indices of columns to operate on and return. It can also be
%                       a cell array of column names of the input table. The default is [] indicating all
%                       columns.
%   'Fill'              When tWin exceeds an epoch, 'none' (default) ignores exceeded parts whereas 'bleed'
%                       looks for events from neighboring epochs up to the entire session.
%   'ReferenceTime'     Reference time to use with the 'bleed' option. The default value is empty and
%                       the method will look for reference time associated with the specified table.
% Output
%   tbOut               The output table with inquired event times data.

% Handle user inputs
p = inputParser();
p.addRequired('tbIn', @(x) ischar(x) || isstring(x) || istable(x));
p.addRequired('tWin', @(x) isnumeric(x) || iscell(x));
p.addOptional('rowInd', [], @(x) isnumeric(x) || islogical(x));
p.addOptional('colInd', [], @(x) isnumeric(x) || islogical(x) || iscellstr(x) || isstring(x));
p.addParameter('Fill', 'none', @(x) any(strcmpi(x, {'none', 'bleed'})));
p.addParameter('ReferenceTime', [], @(x) isnumeric(x) && isvector(x));

p.parse(tbIn, tWin, varargin{:});
rowInd = p.Results.rowInd;
colInd = p.Results.colInd;
fillMethod = lower(p.Results.Fill);
tRef = p.Results.ReferenceTime;

if ~istable(tbIn)
    assert(this.IValidateTableName(tbIn, true) == 1, '''%s'' is not an eventTimes table', tbIn);
    tRef = this.GetReferenceTime(tbIn);
    tbIn = this.GetTable(tbIn);
end

% Validate and standardize row and column indices
[rowInd, colInd] = this.IValidateTableIndexing(tbIn, rowInd, colInd);
tbIn = tbIn(:,colInd);

% Validate and standardize time windows
tWin = this.IValidateTimeWindows(tWin, height(tbIn), rowInd);

% Slicing
switch fillMethod
    case 'none'
        winCells = ISliceEtFillNone(tbIn, tWin, rowInd);
    case 'bleed'
        assert(~isempty(tRef), 'Requires ''ReferenceTime'' to use ''%s''', fillMethod);
        winCells = ISliceEtFillBleed(tbIn, tWin, rowInd, tRef);
end
winCells = cat(1, winCells{:});

% Make table with inherited data types, if possible
tbOut = table();
for k = 1 : width(tbIn)
    varName = tbIn.Properties.VariableNames{k};
    if isa(tbIn.(k), 'cell')
        tbOut.(varName) = winCells(:,k);
    elseif all(cellfun(@isscalar, winCells(:,k)))
        try % this is a hack to handle different datatype
            tbOut.(varName) = cat(1, winCells{:,k});
        catch
            tbOut.(varName) = winCells(:,k);
        end
    else
        tbOut.(varName) = winCells(:,k);
    end
end

end


function winData = ISliceEtFillNone(tbIn, tWin, rowInd)
% Find event times for each epoch
winData = cell(size(tWin));
for i = rowInd
    w = tWin{i};
    winData{i} = cell(size(w,1), width(tbIn));
    
    % Find event times for each variable
    for j = 1 : width(tbIn)
        % Get event times
        t = tbIn{i,j};
        if iscell(t)
            t = t{1};
        end
        
        % Find event times in each window
        for k = 1 : size(w,1)
            tHit = t(t >= w(k,1) & t < w(k,2));
            if isempty(tHit)
                tHit = NaN;
            end
            winData{i}{k,j} = tHit;
        end
    end
end
end

function winData = ISliceEtFillBleed(tbIn, tWin, rowInd, tRef)
% Cache variables
isOldVer = verLessThan('matlab', '9.1');
tAbs = cell(size(tbIn));
for i = 1 : width(tbIn)
    if iscell(tbIn{:,i})
        tAbs(:,i) = cellfun(@(x,r) double(x)+r, tbIn{:,i}, num2cell(tRef), 'Uni', false);
    else
        tAbs(:,i) = num2cell(double(tbIn{:,i}) + tRef);
    end
end
for i = 1 : numel(tAbs)
    if isempty(tAbs{i})
        tAbs{i} = NaN;
    end
end
tAbsBegin = cellfun(@(x) x(1), tAbs);
tAbsEnd = cellfun(@(x) x(end), tAbs);

% Find event times for each epoch
winData = cell(size(tWin));
for i = rowInd
    wAbs = tWin{i} + tRef(i);
    winData{i} = cell(size(wAbs,1), width(tbIn));
    
    % Find event times for each variable
    for j = 1 : width(tbIn)
        % Find event times in each window
        for k = 1 : size(wAbs,1)
            % Collect source epochs
            if isOldVer
                dtWinBeforeEpoch = bsxfun(@minus, tAbsBegin(:,j), wAbs(k,:));
                dtWinAfterEpoch = bsxfun(@minus, wAbs(k,:), tAbsEnd(:,j));
            else
                dtWinBeforeEpoch = tAbsBegin(:,j) - wAbs(k,:);
                dtWinAfterEpoch = wAbs(k,:) - tAbsEnd(:,j);
            end
            isWinBeforeEpoch = all(dtWinBeforeEpoch >= 0, 2);
            isWinAfterEpoch = all(dtWinAfterEpoch > 0, 2);
            isWinOverlapEpoch = ~(isWinBeforeEpoch | isWinAfterEpoch);
            srcRowInd = find(isWinOverlapEpoch);
            
            % Sort source epochs in temporal order
            [~, ord] = sort(tRef(srcRowInd));
            srcRowInd = srcRowInd(ord);
            
            % Get masks
            isInWin = cellfun(@(x) wAbs(k,1) <= x & x < wAbs(k,2), tAbs(srcRowInd,j), 'Uni', false);
            
            % Find times
            tCells = cellfun(@(x,y) x(y)-tRef(i), tAbs(srcRowInd,j), isInWin, 'Uni', false);
            t = vertcat(tCells{:});
            if isempty(t)
                t = NaN;
            end
            winData{i}{k,j} = t;
        end
    end
end
end