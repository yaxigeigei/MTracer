classdef MSessionExplorer < handle
    % MSessionExplorer is a data container that makes data munging easy
    % 
    %   For a tutorial, please check out examples in the following order (just copy a line and run).
    %   Methods without examples should be well explained by help documents alone.
    %   
    %   Object Construction
    %     MSessionExplorer.Examples.MakeEventTimesTable
    %     MSessionExplorer.Examples.MakeTimeSeriesTable
    %     MSessionExplorer.Examples.Duplicate
    %     MSessionExplorer.Examples.Merge
    % 
    %   Content Management
    %     MSessionExplorer.Examples.Peek
    %     MSessionExplorer.Examples.SetTable
    %     MSessionExplorer.Examples.SetReferenceTime
    %     MSessionExplorer.Examples.SetColumn
    % 
    %   Data Operations
    %     MSessionExplorer.Examples.AlignTime
    %     MSessionExplorer.Examples.SliceSession
    %     MSessionExplorer.Examples.SliceEventTimes
    %     MSessionExplorer.Examples.SliceTimeSeries
    %     MSessionExplorer.Examples.ResampleEventTimes
    %     MSessionExplorer.Examples.ResampleTimeSeries
    % 
    % See also MSessionExplorer.Event, MPlot, MPlotter, MMath
    
    properties(Constant)
        supportedTableTypes = ...   % 'eventTimes', 'eventValues', 'timeSeries' (read-only)
            {'eventTimes', 'eventValues', 'timeSeries'};
    end
    properties(GetAccess = public, SetAccess = private)
        tot                         % A table of data tables, reference times and related metadata (read-only)
    end
    properties
        userData                    % A struct where user stores arbitrary data
    end
    properties(Dependent)
        tableNames                  % A cell array of table names (read-only)
        isEventTimesTable           % Whether each table is eventTimes table (read-only)
        isEventValuesTable          % Whether each table is eventValues table (read-only)
        isTimesSeriesTable          % Whether each table is timeSeries table (read-only)
        numEpochs                   % The number of epochs (read-only)
    end
    properties
        epochInd                    % Indices of epoch
        isVerbose logical = true    % Whether or not to display progress and some warnings
    end
    properties(Access = private)
        originalTrialInd            % for backward compatibility
    end
    properties(Dependent, Hidden)
        numTrials                   % for backward compatibility
    end
    
    % These methods are exposed to users
    methods
        % Constructor
        function this = MSessionExplorer()
            % Constructor of MSessionExplorer
            % 
            %   se = MSessionExplorer()
            
            % Initialize table of tables
            totHeaders = {'tableType', 'tableData', 'referenceTime'};
            this.tot = cell2table(cell(0, numel(totHeaders)), 'VariableNames', totHeaders);
        end
        
        % Property IO
        function val = get.tableNames(this)
            val = this.tot.Properties.RowNames;
        end
        function val = get.isEventTimesTable(this)
            val = strcmp(this.tot.tableType, this.supportedTableTypes{1})';
        end
        function val = get.isEventValuesTable(this)
            val = strcmp(this.tot.tableType, this.supportedTableTypes{2})';
        end
        function val = get.isTimesSeriesTable(this)
            val = strcmp(this.tot.tableType, this.supportedTableTypes{3})';
        end
        function val = get.numEpochs(this)
            if ~isempty(this.tot)
                val = height(this.tot.tableData{1});
            else
                val = 0;
            end
        end
        function set.epochInd(this, val)
            assert(all(val > 0 & mod(val,1) == 0), 'Epoch indices must be positive integers');
            assert(numel(val) == numel(unique(val)), 'Epoch indices must be unique numbers');
            assert(numel(val) == this.numEpochs, 'The number of indices must match the number of epochs');
            this.epochInd = val(:);
        end
        function val = get.epochInd(this)
            if isempty(this.epochInd) && this.numEpochs > 0 && ~isempty(this.originalTrialInd)
                warning(['originalTrialInd property will be removed in a future version. ' ...
                    'Use epochInd instead and consider saving the updated object.']);
                this.epochInd = this.originalTrialInd;
            end
            val = this.epochInd;
        end
        function set.isVerbose(this, val)
            assert(islogical(val) && isscalar(val), 'isVerbose must be a logical scalar');
            this.isVerbose = val;
        end
        function val = get.numTrials(this)
            warning('numTrials property will be removed in a future version. Use numEpochs instead.');
            val = this.numEpochs;
        end
        
        % Content Management
        function Peek(this, varargin)
            % Print a summary of the main content of this object or the begining of specific table(s)
            %   
            %   Peek()
            %   Peek(tbName1, tbName2, tbName3, ...)
            % 
            % Inputs
            %   tbNameN         Name of a data table to display. 
            %                   If not specified, this method will display tot and userData. 
            
            if nargin < 2
                disp('tot');
                disp(this.tot);
                disp('userData');
                disp(this.userData);
            else
                this.IValidateTableNames(varargin, true);
                for i = 1 : numel(varargin)
                    tb = this.tot{varargin{i}, 'tableData'}{1};
                    indRow = 1 : min(6, height(tb));
                    disp(tb(indRow,:));
                end
            end
        end
        
        function se = Duplicate(this, varargin)
            % Make a hard copy of the current object
            % 
            %   se = Duplicate()
            %   se = Duplicate(tbNames)
            %   se = Duplicate(tbNames, isUserData)
            % 
            % Inputs
            %   tbNames             A cell array of table names specifying which tables to include. An empty 
            %                       array (default) indicates all tables. 
            %   isUserData          A logical variable indicating whether or not to copy userData. Default true. 
            % Output
            %   se                  The new MSessionExplorer object
            
            if numel(this) > 1
                % Duplicate se array recursively
                se = arrayfun(@(x) x.Duplicate(varargin{:}), this);
                return
            end
            
            % Handle user input
            p = inputParser();
            p.addOptional('tbNames', []);
            p.addOptional('isUserData', true, @islogical);
            p.parse(varargin{:});
            tbNames = p.Results.tbNames;
            isUserData = p.Results.isUserData;
            
            if isempty(tbNames)
                tbNames = this.tableNames;
            end
            this.IValidateTableNames(tbNames, true);
            
            % Find table indices
            tbNames = cellstr(tbNames);
            tbInd = cellfun(@(x) find(strcmp(x, this.tableNames)), tbNames);
            
            % Copying
            se = eval(class(this));
            se.isVerbose = this.isVerbose;
            for i = tbInd(:)'
                se.SetTable(this.tableNames{i}, this.tot.tableData{i}, this.tot.tableType{i}, this.tot.referenceTime{i});
            end
            if isUserData
                se.userData = this.userData;
            end
            se.epochInd = this.epochInd;
        end
        
        function se = Merge(this, varargin)
            % Combine multiple MSessionExplorer objects into one by vertically concatenating data
            % 
            %   se = Merge(se1, se2, se3, ...)
            % 
            % Inputs
            %   se1, se2, se3, ...      Arbitrary number of MSessionExplorer objects or object arrays
            % Output
            %   se                      The merged MSessionExplorer object
            
            % SEs to merge
            for i = 1 : numel(varargin)
                varargin{i} = varargin{i}(:);
            end
            seArray = [this; cat(1, varargin{:})];
            this = seArray(1); % handle the case where the input "this" is a vector
            
            % Output SE
            se = eval(class(this));
            se.isVerbose = this.isVerbose;
            
            % Concatenate and set each table
            for i = 1 : numel(this.tableNames)
                % Get all tables of the same name
                tbName = this.tableNames{i};
                tbs = arrayfun(@(x) x.GetTable(tbName), seArray, 'Uni', false);
                
                % Concatenation
                tbCat = this.ICatTables(tbs, this.isTimesSeriesTable(i));
                
                % Set table to merged se
                refTime = arrayfun(@(x) x.GetReferenceTime(tbName), seArray, 'Uni', false);
                refTime = cat(1, refTime{:});
                se.SetTable(tbName, tbCat, this.tot.tableType{i}, refTime);
            end
            
            % Add incremented epoch indices
            cumNumEp = cumsum(arrayfun(@(x) max(x.epochInd), seArray));
            epInd = arrayfun(@(x) x.epochInd, seArray, 'Uni', false);
            epInd(2:end) = cellfun(@(x,y) x + y, num2cell(cumNumEp(1:end-1)), epInd(2:end), 'Uni', false);
            se.epochInd = cat(1, epInd{:});
            
            % Add all userdata in a struct array, filling any missing fields with []
            ud = arrayfun(@(x) x.userData, seArray, 'Uni', false);
            isEptUd = cellfun(@isempty, ud);
            if all(isEptUd)
                return
            else
                ud(isEptUd) = repmat({struct}, sum(isEptUd));
            end
            fdNames = cellfun(@fieldnames, ud, 'Uni', false);
            uniFdNames = unique(cat(1, fdNames{:}), 'stable');
            for i = 1 : numel(ud)
                for j = 1 : numel(uniFdNames)
                    fn = uniFdNames{j};
                    if ~isfield(ud{i}, fn)
                        ud{i}.(fn) = [];
                    end
                end
            end
            se.userData = cat(1, ud{:});
        end
        
        function seArray = Split(this, epochArg, isSplitUserData)
            % Split epochs of a MSessionExplorer object to individual objects
            % 
            %   seArray = Split(epochDist)
            %   seArray = Split(epochInd)
            %   seArray = Split(..., isSplitUserData)
            % 
            % Inputs
            %   epochDist       An n-element numeric vector for the numbers of epochs to split. 
            %                   This is similar to rowDist in MATLAB's mat2cell function.
            %   epochInd        A cell array where each element is a vector of epoch indices.
            % Output
            %   seArray         A vector of MSessionExplorer objects.
            %
            % See also mat2cell
            
            if nargin < 3
                isSplitUserData = true;
            end
            
            % Convert epochDist to epochInd
            if isnumeric(epochArg)
                epDist = epochArg(:);
                epB = cumsum(epDist);
                epA = [0; epB(1:end-1)] + 1;
                epInd = arrayfun(@(a,b) (a:b)', epA, epB, 'Uni', false);
            else
                epInd = epochArg;
            end
            
            % Initialize se objects with userData
            for k = numel(epInd) : -1 : 1
                seArray(k) = eval(class(this)); % inheritance compatible construction
                seArray(k).isVerbose = this.isVerbose;
                if isSplitUserData && numel(this.userData) == numel(epInd)
                    seArray(k).userData = this.userData(k);
                else
                    seArray(k).userData = this.userData;
                end
            end
            
            % Add tables
            for i = 1 : numel(this.tableNames)
                tbName = this.tableNames{i};
                tb = this.GetTable(tbName);
                rt = this.GetReferenceTime(tbName);
                for k = 1 : numel(seArray)
                    m = epInd{k};
                    if isempty(rt)
                        seArray(k).SetTable(tbName, tb(m,:), this.tot.tableType{i});
                    else
                        seArray(k).SetTable(tbName, tb(m,:), this.tot.tableType{i}, rt(m));
                    end
                end
            end
            
            seArray = reshape(seArray, size(epochArg));
        end
        
        function seTb = SplitConditions(this, condVars, sourceTbName)
            % Split an SE into a table of SEs by conditions in the behavValue table
            % Trials with NaN condition will be excluded except for opto
            % 
            %   seTb = SplitConditions(condVars)
            %   seTb = SplitConditions(condVars, sourceTbName)
            % 
            % Inputs
            %   condVars        1) One or more table column names. Conditions are defined by every unique conbination across 
            %                      the values of these variables.
            %                   2) A table with condition variables as columns.
            %                   3) An empty [] variable. A condition variable called dummyCond will be used to group all epochs 
            %                      in one condition.
            %   sourceTbName    The name of the table to find columns specified by condVars.
            %   
            
            % Construct grouping table
            if isempty(condVars)
                % Initialize table with a dummy grouping variable that includes all epochs in one group
                dummyCond = ones(se.numEpochs, 1);
                T = table(dummyCond);
            elseif istable(condVars)
                % Condition variables are provided as a table
                T = condVars;
            else
                % Get and modify variables from the specified table
                tb = this.GetTable(sourceTbName);
                T = tb(:, condVars);
            end
            assert(height(T)==this.numEpochs, "The number of rows of the condVars table must match the number of epochs.");
            
            % Find epoch indices for each group
            [condId, seTb] = findgroups(T);
            for i = 1 : max(condId)
                m = condId==i;
                seTb.epochInd{i} = find(m);
                seTb.numEpochs(i) = sum(m);
            end
            
            % Split se
            seTb.se = this.Split(seTb.epochInd);
        end
        
        function s = ToStruct(this)
            % Convert all the contents of this object to structures
            % 
            %   s = ToStruct()
            % 
            % Output
            %   s       The output structure with the following fields
            %           tableName       a cell array of table names
            %           tableType       a cell array of table types
            %           referenceTime   a cell array of referece time
            %           tableData       a cell array of structs whose fields are each table's variables 
            %                           (output from MATLAB table2struct function)
            %           epochInd        a vector of epoch indices
            %           userData        same as userData property
            
            s.tableName = this.tableNames;
            s.tableType = this.tot.tableType;
            s.referenceTime = this.tot.referenceTime;
            s.tableData = cellfun(@(x) table2struct(x), this.tot.tableData, 'Uni', false);
            s.epochInd = this.epochInd;
            s.userData = this.userData;
        end
        
        function SetTable(this, tbName, tb, varargin)
            % Add or update a data table
            % 
            %   SetTable(tbName, tb)
            %   SetTable(tbName, tb, tableType)
            %   SetTable(tbName, tb, tableType, referenceTime)
            %
            % Inputs:
            %   tbName              A string of the name of table to add or update. 
            %   tb                  The table variable.
            %   tableType           'eventTimes', 'eventValues', or 'timeSeries' indicating the type of 
            %                       data table. It is required when adding a new table but is ignored when 
            %                       updating an existing table. 
            %   referenceTime       A numeric vector where each element indicate the corresponding session 
            %                       (absolute) time of the time zero in each epoch. Default is empty. 
            
            % Handle user input
            p = inputParser();
            p.addRequired('tbName', @(x) ischar(x) || (isscalar(x) && isstring(x)));
            p.addRequired('tb', @istable);
            p.addOptional('tableType', [], @(x) any(strcmp(x, this.supportedTableTypes)));
            p.addOptional('referenceTime', [], @isnumeric);
            p.parse(tbName, tb, varargin{:});
            tbName = char(tbName);
            tbType = char(p.Results.tableType);
            refTimes = p.Results.referenceTime;
            
            % Check epoch number conflict
            if ~isempty(this.tableNames)
                assert(size(tb,1) == this.numEpochs, ...
                    'The new table cannot be added since it has %d rows whereas existing table has %d.', ...
                    size(tb,1), this.numEpochs);
            end
            
            % Set table
            if ismember(tbName, this.tableNames)
                % Replace an existing table
                this.tot{tbName, 'tableData'}{1} = tb;
                if ~isempty(tbType)
                    this.tot{tbName, 'tableType'}{1} = tbType;
                end
            else
                isFirst = isempty(this.tot);
                
                % Add a new table
                assert(~isempty(tbType), 'Table type is required for adding a new table');
                totHeaders = this.tot.Properties.VariableNames;
                totRow = cell2table(cell(1,3), 'VariableNames', totHeaders, 'RowNames', {tbName});
                totRow.tableType{1} = tbType;
                totRow.tableData{1} = tb;
                totRow.referenceTime{1} = [];
                this.tot = [this.tot; totRow];
                
                % Initialize epoch indices
                if isFirst
                    this.epochInd = (1 : height(tb))';
                end
            end
            
            % Set reference time
            if ~isempty(refTimes)
                this.SetReferenceTime(refTimes, tbName);
            end
        end
        
        function varargout = GetTable(this, varargin)
            % Return specific data table
            % 
            %   [tb1, tb2, tb3, ...] = GetTable(tbName1, tbName2, tbName3, ...)
            %
            % Input
            %   tbNameN         The name of a table to return. 
            % Output
            %   tbN             The table data. 
            
            this.IValidateTableNames(varargin, true);
            for i = numel(varargin) : -1 : 1
                varargout{i} = this.tot{varargin{i}, 'tableData'}{1};
            end
        end
        
        function RemoveTable(this, varargin)
            % Remove specific data tables
            % 
            %   RemoveTable(tbName1, tbName2, tbName3, ...)
            %
            % Input
            %   tbNameN         A string indicating the name of a table to remove.
            
            varargin = cellstr(varargin);
            notTb = setdiff(varargin, this.tableNames);
            if this.isVerbose
                warning backtrace off
                cellfun(@(x) warning('There is no table named %s', x), notTb);
            end
            varargin = intersect(varargin, this.tableNames);
            this.tot(varargin,:) = [];
            if isempty(this.tot)
                this.epochInd = [];
            end
        end
        
        function SetReferenceTime(this, rt, tbNames)
            % Set reference times to table(s)
            % 
            %   SetReferenceTime(rt)
            %   SetReferenceTime(rt, tbNames)
            %
            % Inputs
            %   rt              A numeric vector where each element stores the absolute time of the 
            %                   zero time in each epoch. Or use [] to clear existing reference times. 
            %   tbNames         A string or cell array of the name(s) of table which reference time is 
            %                   set to. The default is empty and rt is set to all eligible tables.
            
            if nargin < 3
                tbNames = this.tableNames(~this.isEventValuesTable);
            end
            this.IValidateTableNames(tbNames, true);
            
            rt = rt(:);
            if ~all(diff(rt) > 0)
                disp('Reference times are not monotonically increasing');
            end
            
            tbNames = cellstr(tbNames);
            for i = 1 : numel(tbNames)
                assert(this.numEpochs == numel(rt) || isempty(rt), ...
                    'The reference time has %d elements which does not match the %d epochs', ...
                    numel(rt), this.numEpochs);
                
                if ~this.isEventValuesTable(strcmp(tbNames{i}, this.tableNames))
                    this.tot{tbNames{i}, 'referenceTime'}{1} = rt;
                else
                    warning('''%s'' is an eventValues table and reference time is not applicable', tbNames{i});
                end
            end
        end
        
        function varargout = GetReferenceTime(this, varargin)
            % Return reference times
            % 
            %   rt = GetReferenceTime()
            %   [rt1, rt2, rt3, ...] = GetReferenceTime(tbName1, tbName2, tbName3, ...)
            %
            % Input
            %   tbNameN         The name of a table which the reference time is associated with. If not 
            %                   specified, the first availble reference time will be returned. 
            % Output
            %   rtN             A vector of reference times. 
            
            if nargin > 1
                this.IValidateTableNames(varargin, true);
            else
                isRt = ~cellfun(@isempty, this.tot.referenceTime);
                assert(any(isRt), 'No table has reference time');
                varargin = this.tableNames(find(isRt,1));
            end
            
            for i = numel(varargin) : -1 : 1
                varargout{i} = this.tot{varargin{i}, 'referenceTime'}{1};
            end
        end
        
        function RemoveEpochs(this, row2rm)
            % Remove specific rows across all tables
            % 
            %   RemoveEpochs(row2rm)
            %
            % Input
            %   row2rm          Integer or logical indices of rows to remove. 
            
            for k = 1 : height(this.tot)
                this.tot.tableData{k}(row2rm,:) = [];
                if ~isempty(this.tot.referenceTime{k})
                    this.tot.referenceTime{k}(row2rm) = [];
                end
            end
            this.epochInd(row2rm) = [];
        end
        
        function colData = GetColumn(this, tbName, colIDs)
            % Return specific columns from a data table
            %
            %   colData = GetColumn(tbName, colIDs)
            %
            % Inputs
            %   tbName          The name of a table that contains the columns of interest. 
            %   colIDs          1) Column name(s) as a string or a cell array of strings. 
            %                   2) Integer or logical indices of column. 
            % Output
            %   colData         Requsted data. 
            
            this.IValidateTableName(tbName, true);
            tb = this.tot{tbName, 'tableData'}{1};
            colData = tb{:,colIDs};
        end
        
        function SetColumn(this, tbName, colIDs, colData, rowScope)
            % Add, update or delete column(s) in a data table
            %
            %   SetColumn(tbName, colIDs, colData)
            %   SetColumn(tbName, colIDs, func, rowScope)
            %
            % Inputs
            %   tbName          The name of a table that contains the columns of interest. 
            %   colIDs          1) Column name(s) as a string or a cell array of strings. 
            %                   2) Integer or logical indices of column. 
            %   colData         Data to add or update. Use [] to delete. 
            %   func            A function handle that receives one input and gives one output. 
            %   rowScope        'each' or 'all'. 'each' applies func to each row separately, whereas 
            %                   'all' applies to all rows (concatenated) as one continuous series. 
            
            this.IValidateTableName(tbName, true);
            tb = this.tot{tbName, 'tableData'}{1};
            
            if isempty(colData)
                % Remove variable
                colInd = this.IValidateColumnIndexing(tb, colIDs);
                colInd = unique(colInd);
                for i = numel(colInd) : -1 : 1
                    tb.(colInd(i)) = [];
                end
            elseif isa(colData, 'function_handle') && nargin > 4
                % Apply user function
                func = colData;
                colInd = this.IValidateColumnIndexing(tb, colIDs);
                for k = 1 : numel(colInd)
                    C = tb.(colInd(k));
                    if strcmp(rowScope, 'each')
                        if iscell(C)
                            C = cellfun(func, C, 'Uni', false);
                        else
                            C = arrayfun(func, C);
                        end
                    elseif strcmp(rowScope, 'all')
                        if iscell(C)
                            L = cellfun(@(x) size(x,1), C);
                            C = func(cell2mat(C));
                            C = mat2cell(C, L);
                        else
                            C = func(C);
                        end
                    else
                        error('scope must be ''each'' or ''all'' but was ''%s''', rowScope);
                    end
                    tb.(colInd(k)) = C;
                end
            else
                % Set values
                tb{:,colIDs} = colData;
            end
            
            this.tot{tbName, 'tableData'}{1} = tb;
        end
        
        function Column2Cell(this, tbName, colIDs)
            % Force columns to be cell arrays
            %
            %   Column2Cell(tbName)
            %   Column2Cell(tbName, colIDs)
            %
            % Inputs
            %   tbName          The name of a table that contains the columns of interest.
            %   colIDs          1) Column name(s) as a string or a cell array of strings.
            %                   2) Integer or logical indices of column.
            %                   3) If this is not provided or empty, all columns will be included.
            
            this.IValidateTableName(tbName, true);
            tb = this.tot{tbName, 'tableData'}{1};
            
            if nargin < 3 || isempty(colIDs)
                colIDs = 1 : width(tb);
            end
            
            colInd = this.IValidateColumnIndexing(tb, colIDs);
            for k = 1 : numel(colInd)
                C = tb.(colInd(k));
                if ~iscell(C)
                    C = num2cell(C);
                end
                tb.(colInd(k)) = C;
            end
            
            this.tot{tbName, 'tableData'}{1} = tb;
        end
        
        % Data Operations
        function AlignTime(this, et, sourceTbName)
            % Align epochs by changing the origin of time in each epoch to the specified event times
            % 
            %   AlignTime(et)
            %   AlignTime(et, sourceTbName)
            %
            % Inputs
            %   et              An event name in an eventTimes table or a numeric vector of times to align 
            %                   to. Each time in the vector is relative wrt each row (epoch), which becomes 
            %                   the new zero time after alignment. 
            %   sourceTbName    If et is a string, then you must specify the name of an eventTimes table 
            %                   where this event name should be found. This avoids ambiguity of the same 
            %                   variable name found in multiple tables. The default is empty. 
            
            % Check eligible tables
            indAlignable = find(this.isEventTimesTable | this.isTimesSeriesTable);
            assert(~isempty(indAlignable), 'Requires at least one eventTimes or timeSeries table to operate on');
            
            % Get reference event times
            if ischar(et) || isstring(et)
                assert(nargin == 3, 'Requires refSourceTableName to indicate where ''%s'' is in', et);
                this.IValidateTableName(sourceTbName, true);
                et = this.tot{sourceTbName, 'tableData'}{1}.(et);
            end
            if isscalar(et)
                et = repmat(et, [this.numEpochs 1]);
            end
            
            % Validate reference event times
            assert(isnumeric(et) && isvector(et), 'Reference times must be a numeric vector');
            assert(~any(isnan(et)), 'Reference times cannot have NaN');
            assert(numel(et) == this.numEpochs, ...
                'The number of reference times (%d) does not match the number of rows (%d) in data table.', ...
                numel(et), this.numEpochs);
            et = et(:);
            
            % Align times
            for k = indAlignable
                if this.isEventTimesTable(k)
                    % For eventTimes table
                    tb = this.tot.tableData{k};
                    for i = 1 : size(tb, 2)
                        if iscell(tb{:,i})
                            % Cell vector of numeric vectors
                            tb{:,i} = cellfun(@(x,r) x-r, tb{:,i}, num2cell(et), 'Uni', false);
                        else
                            % Numeric vector
                            tb{:,i} = tb{:,i} - et;
                        end
                    end
                    this.tot.tableData{k} = tb;
                    
                elseif this.isTimesSeriesTable(k)
                    % For timeSeries table
                    this.tot.tableData{k}.time = cellfun(@(x,r) x-r, ...
                        this.tot.tableData{k}.time, num2cell(et), 'Uni', false);
                end
                
                % Change referenceTime
                if ~isempty(this.tot.referenceTime{k})
                    this.tot.referenceTime{k} = this.tot.referenceTime{k} + et;
                end
            end
        end
        
        function SortEpochs(this, ind)
            % Sort epochs across all tables
            % 
            %   SortEpochs()
            %   SortEpochs(ind)
            %
            % Input
            %   ind             Indices of the new order. If not specified, epochs are sorted back to the 
            %                   original order based on the epochInd property. 
            
            % Handle user inputs
            if nargin < 2
                [~, ind] = sort(this.epochInd);
            end
            ind = unique(ind(:), 'stable');
            assert(numel(ind) == this.numEpochs, ...
                'The number of unique indices (%d) does not equal to the number of epochs (%d)', ...
                numel(ind), this.numEpochs);
            
            % Process original indices
            this.epochInd = this.epochInd(ind);
            
            % Sort epochs
            for k = 1 : height(this.tot)
                this.tot.tableData{k} = this.tot.tableData{k}(ind,:);
                if ~isempty(this.tot.referenceTime{k})
                    this.tot.referenceTime{k} = this.tot.referenceTime{k}(ind);
                end
            end
        end
        
        function SliceSession(this, tSlice, tType)
            % Slice session into a different set of epochs
            % 
            %   SliceSession(tSlice, tType)
            %
            % Inputs
            %   tSlice      A vector of times at which slicing occurs. Note that data before the first slicing 
            %               time will be irreversibly discarded. 
            %               If refType is 'relative', tSlice can be a scalar which is used by all epochs, or 
            %               a vector whose length matches the number epoch. 
            %               If refType is 'absolute', tSlice can be a vector of arbitrary length. 
            %   tType       'absolute' or 'relative', indicating whether tSlice represents absolute times in 
            %               the session or is relative to epoch reference times. If using 'absolute', existing 
            %               epoch sorting, time alignment and all eventValues table will be lost due to an 
            %               unknown relationship between current and new epochs. New data will sort epochs in 
            %               temporal order and align epoch time to respective tSlice. 
            
            % Validate inputs
            indTimeTable = find(this.isEventTimesTable | this.isTimesSeriesTable);
            assert(~isempty(indTimeTable), 'Requires at least one eventTimes or timeSeries table to operate on');
            assert(isnumeric(tSlice) && isvector(tSlice), 'Slicing times must be a numeric vector');
            assert(~any(isnan(tSlice)), 'Slicing times cannot have NaN');
            assert(ismember(tType, {'absolute', 'relative'}), ...
                'refType must be ''absolute'' or ''relative'' but instead was ''%s''', tType);
            
            % Restore epoch order
            [~, indBack] = sort(this.epochInd);
            this.SortEpochs();
            
            % Loop through tables
            tSlice = tSlice(:);
            for k = indTimeTable
                % Show table name
                if this.isVerbose
                    fprintf('%s\t', this.tableNames{k});
                end
                
                % Get reference times
                tRef = this.tot.referenceTime{k};
                assert(~isempty(tRef), 'To reslice, a table must have associated reference times');
                
                % Convert slicing times
                if strcmp(tType, 'absolute')
                    tDelim = tSlice;
                else
                    if isscalar(tSlice)
                        tSlice = repmat(tSlice, [this.numEpochs 1]);
                    end
                    tDelim = tSlice + tRef;
                end
                
                % Slice table
                tb = this.tot.tableData{k};
                vect = cell(1,width(tb));
                if this.isEventTimesTable(k)
                    % For eventTimes table
                    for i = 1 : width(tb)
                        vect{i} = this.ICatColumn(tb.(i), tRef);
                    end
                    this.tot.tableData{k} = this.MakeEventTimesTable(vect, ...
                        'DelimiterTimes', tDelim, ...
                        'VariableNames', tb.Properties.VariableNames, ...
                        'Verbose', this.isVerbose);
                    
                elseif this.isTimesSeriesTable(k)
                    % For timeSeries table
                    vect{1} = this.ICatColumn(tb.time, tRef);
                    for i = 2 : width(tb)
                        vect{i} = this.ICatColumn(tb.(i));
                    end
                    this.tot.tableData{k} = this.MakeTimeSeriesTable(vect(1), vect(2:end), ...
                        'DelimiterTimes', tDelim, ...
                        'VariableNames', tb.Properties.VariableNames, ...
                        'Verbose', this.isVerbose);
                end
                this.tot.referenceTime{k} = tDelim;
            end
            
            if strcmp(tType, 'absolute')
                % Remove any eventValues table
                if any(this.isEventValuesTable)
                    warning('All eventValues table were removed when using ''%s'' times', tType);
                    this.tot(this.isEventValuesTable,:) = [];
                end
                % Reset epochInd
                this.epochInd = (1 : this.numEpochs)';
            else
                % Sort epoch order back
                this.AlignTime(-tSlice);
                this.SortEpochs(indBack);
            end
        end
        
        function tbOut = ResampleTimeSeries(this, tbIn, tEdges, varargin)
            % Resample timeSeries table by interpolation
            % 
            %   tbOut = ResampleTimeSeries(tbIn, tEdges)
            %   tbOut = ResampleTimeSeries(tbIn, tEdges, rowInd)
            %   tbOut = ResampleTimeSeries(tbIn, tEdges, rowInd, colInd)
            %   tbOut = ResampleTimeSeries(..., 'Method', 'linear')
            %   tbOut = ResampleTimeSeries(..., 'Method', 'linear', 'Extrapolation', 'none')
            %   tbOut = ResampleTimeSeries(..., 'Antialiasing', false)
            % 
            % Inputs
            %   tbIn            A table of time series data or a name of a timeSeries table in the current object.
            %   tEdges          Edges of time bins. This can be one numeric vector that defines edges for every 
            %                   rows, or a cell array of vectors where each applies to a corresponding row. The 
            %                   number of element in cell array must equal the height of the table or the number 
            %                   of selected rows. 
            %   rowInd          Integer or logical indices of rows to operate on and return. The default value is 
            %                   empty indicating all rows. 
            %   colInd          Integer or logical indices of columns to operate on and return. It can also be 
            %                   a cell array of column names of the input table. The default is empty indicating 
            %                   all columns. 
            %   'Method' and 'Extrapolation'
            %                   Use these Name-Value pairs to customize the behavior of interpolation. Please see 
            %                   options of the MATLAB griddedInterpolant function for details. 
            %   'Antialiasing'  Whether or not to resample with antialiasing filters (default is false). 
            %                   If tEdges do not have uniform bin sizes, the antialiasing targets the frequency 
            %                   determined by the smallest bin size. 
            %                   Antialiasing is always ignored when upsampling, or (in the case of non-uniform 
            %                   sampling) when the largest bin size of the query is smaller than the smallest 
            %                   sample duration of the timeseries to be resampled.
            % Output
            %   tbOut           The output table of time series data where each value is the number of occurance. 
            %
            % See also griddedInterpolant
            
            % Parse inputs
            p = inputParser();
            p.addRequired('tbIn', @(x) ischar(x) || isstring(x) || istable(x));
            p.addRequired('tEdges', @(x) iscell(x) || isnumeric(x));
            p.addOptional('rowInd', [], @(x) isnumeric(x) || islogical(x));
            p.addOptional('colInd', [], @(x) isnumeric(x) || islogical(x) || ischar(x) || iscellstr(x) || isstring(x));
            p.addParameter('Method', 'linear');
            p.addParameter('Extrapolation', 'none');
            p.addParameter('Antialiasing', false, @islogical);
            p.parse(tbIn, tEdges, varargin{:});
            rowInd = p.Results.rowInd;
            colInd = p.Results.colInd;
            interpMethod = p.Results.Method;
            extrapMethod = p.Results.Extrapolation;
            isAA = p.Results.Antialiasing;
            
            % Get table
            if ~istable(tbIn)
                assert(this.IValidateTableName(tbIn, true) == 3, '''%s'' is not a timeSeries table', tbIn);
                tbIn = this.GetTable(tbIn);
            end
            
            % Validate and standardize row and column indices
            [rowInd, colInd] = this.IValidateTableIndexing(tbIn, rowInd, colInd);
            if ~ismember(1, colInd)
                colInd = [1 colInd]; % make sure the time column is always included
            end
            
            % Validate and standardize time edges
            tEdges = this.IValidateTimeEdges(tEdges, height(tbIn), rowInd);
            
            % Select rows and columns
            tbIn = tbIn(rowInd, colInd);
            tEdges = tEdges(rowInd);
            
            % Interpolate time series
            tbOut = tbIn;
            for i = 1 : height(tbIn)
                % Get timestamps
                t = tbIn.time{i};
                tq = tEdges{i}(1:end-1) + diff(tEdges{i})/2;
                tbOut.time{i} = tq;
                
                if isAA && numel(t) > 1
                    % Find input and query sampling frequency for antialiasing
                    maxFs = 1 / min(diff(t));
                    minFsq = 1 / max(diff(tEdges{i}));
                end
                
                for j = 2 : width(tbIn)
                    t = tbIn.time{i};
                    v = tbIn.(j){i};
                    dtype = class(v);
                    
                    % Antialiasing
                    if isAA && sum(~isnan(v)) > 1 && minFsq < maxFs
                        [v, t] = resample(double(v), t, minFsq, 'Dimension', 1);
                    end
                    
                    % Interpolation
                    if numel(t) < 2
                        v = NaN(size(tq));
                    else
                        F = griddedInterpolant(t, double(v), interpMethod, extrapMethod);
                        v = F(tq);
                    end
                    
                    tbOut.(j){i} = cast(v, dtype);
                end
            end
        end
        
        function tbOut = ResampleEventTimes(this, tbIn, tEdges, varargin)
            % Resample an eventTimes table to a timeSeries table by counting events in each bin
            % 
            %   tbOut = ResampleEventTimes(tbIn, tEdges)
            %   tbOut = ResampleEventTimes(tbIn, tEdges, rowInd)
            %   tbOut = ResampleEventTimes(tbIn, tEdges, rowInd, colInd)
            %   tbOut = ResampleEventTimes(..., 'Normalization', 'count')
            %   
            % Inputs
            %   tbIn            A table of event times data or a name of an eventTimes table in the current object.
            %   tEdges          Edges of time bins. This can be one numeric vector that defines edges for every 
            %                   rows, or a cell array of vectors where each applies to a corresponding row. The 
            %                   number of element in cell array must equal the height of the table or the number 
            %                   of selected rows. 
            %   rowInd          Integer or logical indices of rows to operate on and return. The default value is 
            %                   empty indicating all rows. 
            %   colInd          Integer or logical indices of columns to operate on and return. It can also be 
            %                   a cell array of column names of the input table. The default is empty indicating 
            %                   all columns. 
            %   'Normalization' See MATLAB histcounts function for detail.
            % Output
            %   tbOut           The output table of time series data where each value is the number of occurance. 
            %
            % See also histcounts
            
            % Parse inputs
            p = inputParser();
            p.addRequired('tbIn', @(x) ischar(x) || istable(x));
            p.addRequired('tEdges', @(x) iscell(x) || isnumeric(x));
            p.addOptional('rowInd', [], @(x) isnumeric(x) || islogical(x));
            p.addOptional('colInd', [], @(x) isnumeric(x) || islogical(x) || ischar(x) || iscellstr(x) || isstring(x));
            p.addParameter('Normalization', 'count');
            p.parse(tbIn, tEdges, varargin{:});
            rowInd = p.Results.rowInd;
            colInd = p.Results.colInd;
            normMethod = p.Results.Normalization;
            
            % Get table
            if ~istable(tbIn)
                assert(this.IValidateTableName(tbIn, true) == 1, '''%s'' is not an eventTimes table', tbIn);
                tbIn = this.GetTable(tbIn);
            end
            assert(~ismember('time', tbIn.Properties.VariableNames), ...
                'The input table cannot contain variable named ''time'' since the output table requires it');
            
            % Validate and standardize row and column indices
            [rowInd, colInd] = this.IValidateTableIndexing(tbIn, rowInd, colInd);
            
            % Validate and standardize time edges
            tEdges = this.IValidateTimeEdges(tEdges, height(tbIn), rowInd);
            
            % Bin event times
            tbOut = tbIn(rowInd, colInd);
            tEdges = tEdges(rowInd);
            for k = 1 : width(tbOut)
                et = tbOut{:,k};
                if isnumeric(et)
                    et = num2cell(et);
                end
                ts = cell(size(et));
                for i = 1 : numel(et)
                    ts{i} = histcounts(et{i}, tEdges{i}, 'Normalization', normMethod)';
                end
                tbOut.(k) = ts;
            end
            
            % Add time column
            tbOut.time = cellfun(@(x) x(1:end-1)+diff(x)/2, tEdges, 'Uni', false);
            tbOut = tbOut(:,[end, 1:end-1]);
        end
        
        tbOut = SliceTimeSeries(this, tbIn, tWin, varargin)
        tbOut = SliceEventTimes(this, tbIn, tWin, varargin)
        items = Plot(this, varargin)
    end
    methods(Static)
        [tb, preTb] = MakeTimeSeriesTable(t, s, varargin)
        [tb, preTb] = MakeEventTimesTable(et, varargin)
        se = UpdateOldObject(se)
    end
    
    % These methods are used internally
    methods(Hidden, Access = protected)
        function val = IValidateTableName(this, tbName, isAssert)
            % The table name must be a string
            val = ischar(tbName) || (isstring(tbName) && isscalar(tbName));
            assert(~isAssert || val, 'A table name must be char or a single string rather than %s.', class(tbName));
            if ~val
                return;
            end
            % Return 1, 2, 3 for the three table types, 0 if the table does not exist
            tbTypes = this.isEventTimesTable + this.isEventValuesTable*2 + this.isTimesSeriesTable*3;
            val = strcmp(tbName, this.tableNames) .* tbTypes';
            val = sum(val);
            assert(~isAssert || val, 'A table named ''%s'' does not exist', tbName);
        end
        
        function val = IValidateTableNames(this, tbNames, isAssert)
            if iscell(tbNames)
                % A cell array of table names
                val = cellfun(@(x) this.IValidateTableName(x, isAssert), tbNames);
            else
                % One single table name
                val = this.IValidateTableName(tbNames, isAssert);
            end
        end
        
        function colInd = IValidateColumnIndexing(~, tb, colIds)
            % Validate and convert string and logical indexing to integer indices
            if ischar(colIds) || iscellstr(colIds)
                colIds = cellstr(colIds);
                for i = numel(colIds) : -1 : 1
                    mask = strcmp(colIds{i}, tb.Properties.VariableNames);
                    assert(any(mask), 'Column ''%s'' cannot be found in the table', colIds{i});
                    colInd(i) = find(mask);
                end
            elseif islogical(colIds)
                assert(numel(colIds) == width(tb), ...
                    'The number of logical indices (%d) do not match the number of table columns (%d)', ...
                    numel(colIds), width(tb));
                colInd = find(colIds);
            else
                colInd = colIds;
            end
            colInd = colInd(:)';
        end
        
        function [rowInd, colInd] = IValidateTableIndexing(~, tb, rowInd, colInd)
            % Validate and convert row and column indices to integers 
            
            if isempty(rowInd)
                % include all rows
                rowInd = 1 : height(tb);
            elseif islogical(rowInd)
                % convert to numerical indices
                rowInd = find(rowInd);
            end
            assert(isnumeric(rowInd), 'Cannot interpret row indexing');
            
            if isempty(colInd)
                % include all columns
                colInd = 1 : width(tb);
            elseif islogical(colInd)
                % convert to numerical indices
                colInd = find(colInd);
            elseif ischar(colInd) || iscellstr(colInd) || isstring(colInd)
                % find column ind based on variable names
                colInd = cellstr(colInd);
                varNames = tb.Properties.VariableNames;
                colNames = colInd;
                for i = 1 : numel(colInd)
                    colInd{i} = find(strcmp(colNames{i}, varNames));
                    assert(~isempty(colInd{i}), '%s is not a valid column name in the table.', colNames{i});
                end
                colInd = cell2mat(colInd);
            end
            assert(isnumeric(colInd), 'Cannot interpret column indexing');
            
            % Reshape as row vectors
            rowInd = rowInd(:)';
            colInd = colInd(:)';
        end
        
        function winOut = IValidateTimeWindows(~, winIn, tbHeight, rowInd)
            % Validate and convert time windows to a cell array of window matrices
            
            if isnumeric(winIn)
                winIn = num2cell(winIn, 2);
            end
            winIn = winIn(:);
            
            winOut = num2cell(NaN(tbHeight,2), 2);
            if numel(winIn) == 1
                % Propagate value to all epochs
                winOut(rowInd) = repmat(winIn, [numel(rowInd) 1]);
            elseif numel(winIn) == numel(rowInd)
                % Add windows to selected rows
                winOut(rowInd) = winIn;
            elseif numel(winIn) == tbHeight
                % Full size windows
                winOut = winIn;
            else
                error('Incorrect number of time windows');
            end
            
            for i = 1 : numel(winOut)
                % Verify array size
                assert(~isempty(winOut{i}), 'Time window cannot be empty. Consider using [NaN NaN].');
                assert(size(winOut{i},2) == 2, 'Time window array must have 2 elements in each row');
                % Verify time increment (ignore NaN windows)
                dt = diff(winOut{i}, 1, 2);
                assert(all(isnan(dt) | dt > 0), 'The window end time must be greater than the begin time');
            end
        end
        
        function edgeOut = IValidateTimeEdges(~, edgeIn, tbHeight, rowInd)
            % Verify time edges and turn it into a standard format
            
            edgeOut = cell(tbHeight, 1);
            if isnumeric(edgeIn)
                % Propagate edges to all rows
                edgeOut(rowInd) = repmat({edgeIn}, [numel(rowInd) 1]);
            elseif numel(edgeIn) == numel(rowInd)
                % Convert edges array for full table
                edgeOut(rowInd) = edgeIn;
            elseif numel(edgeOut) == tbHeight
                % Use full edges
                edgeOut = edgeIn;
            else
                error('Incorrect size of time edges array');
            end
            edgeOut = cellfun(@(x) x(:), edgeOut, 'Uni', false);
            
            for i = rowInd
                x = edgeOut{i};
                assert(isnumeric(x) && isvector(x), 'Each time edges must be a numeric vector');
                assert(all(diff(x) > 0), 'Time edges must be increasing numbers');
            end
        end
        
        function C = ICatColumn(~, C, tRef)
            if nargin < 3
                tRef = zeros(size(C));
            end
            if ~iscell(C)
                % Numeric vector
                C = C + tRef;
            else
                % Cell vector of numeric array
                for i = 1 : numel(C)
                    C{i} = C{i} + tRef(i);
                end
                C = cat(1, C{:});
            end
        end
        
        function tbCat = ICatTables(~, tbs, isTimeseries)
            % Vertically concatenates tables and handles inconsistent column variable names.
            % Each unique variable across tables will be a column in the concatenated table. 
            % Missing parts will be filled with matching NaN arrays for timeseries tables 
            % (i.e. isTimeseries == true), otherwise leaving automatic placeholder (e.g. []).
            
            if nargin < 3
                isTimeseries = false;
            end
            
            warning('off', 'MATLAB:table:RowsAddedExistingVars');
            tbCat = table;
            r = 0; % will store the end index of the last concatenated table
            for i = 1 : numel(tbs)
                vn = tbs{i}.Properties.VariableNames;
                h = height(tbs{i});
                for j = 1 : numel(vn)
                    tbCat.(vn{j})(r+1:r+h) = tbs{i}.(vn{j}); % currently cannot assign inconsistent data type
                end
                r = r + h;
            end
            warning('on', 'MATLAB:table:RowsAddedExistingVars');
            
            if ~isTimeseries
                return
            end
            nSp = cellfun(@numel, tbCat.time);
            for c = 2 : width(tbCat)
                nSig = max(cellfun(@(x) size(x,2), tbCat.(c)));
                for r = 1 : height(tbCat)
                    if isempty(tbCat.(c){r})
                        tbCat.(c){r} = NaN(nSp(r), nSig);
                    end
                end
            end
        end
        
        function IWarnBeta(this)
            if this.isVerbose
                warning('This method or option is in beta version and should only be used in exploratory analysis');
            end
        end
    end
    methods(Static, Access = private)
        function [T, L] = IDelimitTimestamps(t, d)
            % Find number of samples for each epoch based on timestamps and delimiter times
            d = [d(:); Inf];        % delimiter times appended by Inf
            L = zeros(size(d));     % epoch lengths
            T = cell(size(d));      % timestamps in epochs
            isIn = false(size(t));
            for i = 1 : numel(d)
                isBefore = t < d(i);
                T{i} = t(isBefore & ~isIn);
                L(i) = numel(T{i});
                isIn = isBefore;
            end
            for i = 2 : numel(T)
                T{i} = T{i} - d(i-1); % convert to relative time
            end
        end
    end
    
    methods(Hidden)
        % These methods are under development
        % none
        
        % These methods are for backward compatibility
        function SortTrials(this, ind)
            warning('SortTrials method will be removed in a future version. Use SortEpochs instead.');
            this.SortEpochs(ind);
        end
        function RemoveTrials(this, ind2rm)
            warning('RemoveTrials method will be removed in a future version. Use RemoveEpochs instead.');
            this.RemoveEpochs(ind2rm);
        end
        
        % Hide methods from handle superclass
        function listener(~)
        end
        function addlistener(~)
        end
        function notify(~)
        end
        function findobj(~)
        end
        function findprop(~)
        end
    end
end


