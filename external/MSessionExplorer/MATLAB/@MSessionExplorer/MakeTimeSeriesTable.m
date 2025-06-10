function [tb, preTb] = MakeTimeSeriesTable(t, s, varargin)
% Organize time series data into a 'timeSeries' table. 
% Specifically, each row of this table contains data from a unique period of time, or 
% an epoch (e.g. a trial). The first column named 'time' is reserved to contain vectors 
% of timestamp. Other columns contain the corresponding samples. Epochs that do not have 
% any data are left empty. 
% 
%   [tb, preTb] = MSessionExplorer.MakeTimeSeriesTable(t, s)
%   [tb, preTb] = MSessionExplorer.MakeTimeSeriesTable(..., 'DelimiterTimes', [])
%   [tb, preTb] = MSessionExplorer.MakeTimeSeriesTable(..., 'VariableNames', [])
%   [tb, preTb] = MSessionExplorer.MakeTimeSeriesTable(..., 'Verbose', true)
% 
% Inputs
%   t       1) A vector of timestamps. 
%           2) A cell array of 1), each for one epoch. 'DelimeterTimes' is not supported 
%              for this input. 
%   s       1) An numeric array of samples. The number of rows should match the number 
%              of timestamps. The number of columns can be 1 or more. 
%           2) A 1-D cell array of 1). Each cell is treated as a different signal. 
%           3) A 2-D cell array of 1) where columns are different signals and rows 
%              are different epochs. 'DelimeterTimes' is not supported for this input. 
%   'DelimiterTimes'
%           A vector of time values indicating when to cut data into different epochs. 
%           Data before the first delimiter time is separately stored in preTb. 
%           Timestamps will be converted to relative times wrt the preceeding delimiter 
%           time. The default value is [] which performs no delimiting. 
%   'VariableNames'
%           A cell array of column names for the output tables. 'time' is reserved for 
%           the timestamp column and is optional to specify. The default names are 
%           {'time', 'series1', 'series2', ...}. 
%   'Verbose'
%           A logical value that controls the display of progress. Default is true. 
% Outputs
%   tb      A table that works in MSessionExplorer objects as a 'timeSeries' table. 
%   preTb   Similar to tb but includes data before the first delimiter time, if any. 

% Parse inputs
p = inputParser();
p.addRequired('t', @(x) isnumeric(x) || iscell(x));
p.addRequired('s', @(x) isnumeric(x) || iscell(x));
p.addParameter('DelimiterTimes', [], @isnumeric);
p.addParameter('VariableNames', []);
p.addParameter('Verbose', true, @islogical);
p.parse(t, s, varargin{:});
delimiterTimes = p.Results.DelimiterTimes;
varNames = p.Results.VariableNames;
isVerbose = p.Results.Verbose;

% Unify to cell array
if ~iscell(t)
    t = {t};
end
if ~iscell(s)
    s = {s};
end

% Verify timestamp vectors
for i = 1 : numel(t)
    if isempty(t{i})
        continue;
    end
    assert(isvector(t{i}), 'Timestamps must be in vectors');
    if isVerbose && ~all(diff(t{i}) > 0)
        warning('Values in the #%d timestamp vector are not monotonically increasing', i);
    end
    if isrow(t{i})
        t{i} = t{i}';
    end
end

% Make variable names
if isempty(varNames)
    varNames = [{'time'}, arrayfun(@(x) ['series' num2str(x)], 1:size(s,2), 'Uni', false)];
elseif ~strcmp('time', varNames{1})
    varNames = [{'time'}; varNames(:)];
end

if isempty(delimiterTimes)
    % Each cell is an epoch
    if isVerbose
        disp('Delimiter is not specified. Each cell is treated as an epoch.');
    end
    assert(numel(t) == size(s,1), ...
        'There are %d vectors of timestamp and but %d epochs of data', numel(t), size(s,1));
    
    tEpoch = [cell(1); t(:)];
    sEpoch = [cell(1,size(s,2)); s];
else
    % Delimit signals by delimiter times
    if isVerbose
        fprintf('Delimit %d time series using %d delimiter times\n', numel(s), numel(delimiterTimes));
    end
    assert(all(diff(delimiterTimes) > 0), 'Delimiter times must be monotonically increasing');
    assert(isvector(s), 'Cannot delimit time series data in cell array with more than one dimension');
    assert(numel(t) == 1, 'Cannot delimit more than one vector of timestamps');
    
    % Group data
    [tEpoch(:,i), L] = MSessionExplorer.IDelimitTimestamps(t{1}, delimiterTimes);
    for i = numel(s) : -1 : 1
        sEpoch(:,i) = mat2cell(s{i}, L);
    end
end

% Put data into table
preTb = cell2table([tEpoch(1), sEpoch(1,:)], 'VariableNames', varNames);
tb = cell2table([tEpoch(2:end), sEpoch(2:end,:)], 'VariableNames', varNames);

end




