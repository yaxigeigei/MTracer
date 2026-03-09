function [tb, preTb] = MakeEventTimesTable(et, varargin)
% Organize event times data into an 'eventTimes' table. 
% Specifically, each row of this table contains data from a unique period of time, or 
% an epoch (e.g. a trial). Each column contains timestamps of a given event. Epochs 
% that do not have any data are filled with NaN (or an empty object). 
% 
%   [tb, preTb] = MSessionExplorer.MakeEventTimesTable(et)
%   [tb, preTb] = MSessionExplorer.MakeEventTimesTable(..., 'DelimiterTimes', [])
%   [tb, preTb] = MSessionExplorer.MakeEventTimesTable(..., 'VariableNames', [])
%   [tb, preTb] = MSessionExplorer.MakeEventTimesTable(..., 'UniformOutput', true)
%   [tb, preTb] = MSessionExplorer.MakeEventTimesTable(..., 'Verbose', true)
% 
% Inputs
%   et      1) A numeric vector of timestamps of an event. 
%           2) A vector of objects of MSessionExplorer.Event class or superclass. 
%           3) A 1-D cell array of 1) or 2). Each cell is treated as a different event. 
%           4) A 2-D cell array of 1) or 2) where columns are different events and rows 
%              are different epochs. 'DelimeterTimes' is not supported for this input. 
%   'DelimiterTimes'
%           A vector of time values indicating when to cut data into different epochs. 
%           Data before the first delimiter time is separately stored in preTb. 
%           Timestamps will be converted to relative times wrt the preceeding delimiter 
%           time. The default value is [] which performs no delimiting. 
%   'VariableNames'
%           A cell array of column names for the output table. The default names are 
%           {'event1', 'event2', 'event3', ...}. 
%   'UniformOutput'
%           A logical value that indicates whether to allow a column to be a numeric vector
%           if all elements are scalar. If set to false, all columns will be cell arrays.
%   'Verbose'
%           A logical value that controls the display of progress. Default is true. 
% Outputs
%   tb      A table that works in MSessionExplorer objects as an 'eventTimes' table. 
%   preTb   Similar to tb but includes data before the first delimiter time, if any. 

% Parse inputs
p = inputParser();
p.addRequired('et', @(x) isnumeric(x) || iscell(x));
p.addParameter('DelimiterTimes', [], @isnumeric);
p.addParameter('VariableNames', []);
p.addParameter('UniformOutput', true, @islogical);
p.addParameter('Verbose', true, @islogical);
p.parse(et, varargin{:});
delimiterTimes = p.Results.DelimiterTimes;
varNames = p.Results.VariableNames;
isUni = p.Results.UniformOutput;
isVerbose = p.Results.Verbose;

% Unify to cell array
if ~iscell(et)
    et = {et};
end

% Verify timestamp vectors
for i = 1 : numel(et)
    if isempty(et{i})
        continue;
    end
    assert(isvector(et{i}), 'Event times must be in vectors');
    if isVerbose
        etReal = et{i}(~isnan(et{i}));
        if ~all(diff(etReal) > 0)
            warning('Values in the #%d event time vector are not monotonically increasing', i);
        end
    end
    if isrow(et{i})
        et{i} = et{i}';
    end
end

% Make variable names
if isempty(varNames)
    varNames = arrayfun(@(x) ['event' num2str(x)], 1:size(et,2), 'Uni', false);
end

if isempty(delimiterTimes)
    % Each cell is an epoch
    if isVerbose
        disp('Delimiter is not specified. Each cell is treated as an epoch.');
    end
    etEpoch = [cell(1, size(et,2)); et];
    
else
    % Delimit event times by delimiter times
    if isVerbose
        fprintf('Delimit %d events using %d delimiter times\n', numel(et), numel(delimiterTimes));
    end
    assert(all(diff(delimiterTimes) > 0), 'Delimiter times must be monotonically increasing');
    assert(isvector(et), 'Cannot delimit event times in cell array with more than one dimension');
    
    for i = numel(et) : -1 : 1
        etEpoch(:,i) = MSessionExplorer.IDelimitTimestamps(et{i}, delimiterTimes);
    end
end

% Fill empty cells
etEpoch = FillEmpty(etEpoch);

% Put data into tables
preTb = cell2table(etEpoch(1,:), 'VariableNames', varNames);
tb = cell2table(etEpoch(2:end,:), 'VariableNames', varNames);

% Convert any numeric arrays to cell arrays
if ~isUni
    for i = 1 : width(tb)
        if isnumeric(tb.(i))
            preTb.(i) = num2cell(preTb.(i));
            tb.(i) = num2cell(tb.(i));
        end
    end
end

end

function C = FillEmpty(C)
% Fill empty cells
% If all values in a cloumn are numeric, fill empty cells with NaN
% Otherwise fill with an object from the non-numeric class (constructed with no input)

for j = 1 : size(C,2)
    % Cache values
    indEpt = find(cellfun(@isempty, C(:,j)));
    if isempty(indEpt)
        continue;
    end
    isNum = cellfun(@isnumeric, C(:,j));
    objIdx = find(~isNum, 1);
    
    for i = indEpt'
        if isempty(objIdx)
            % Fill NaN
            C{i,j} = NaN;
        else
            % Fill an empty object
            className = class(C{objIdx,j});
            C{i,j} = eval(className);
        end
    end
end

end

