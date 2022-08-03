classdef Event
    % MSessionExplorer.Event binds times and values of an event together and works like a numeric event 
    % time in eventTimes table of MSessionExplorer objects
    
    properties(Dependent)
        t;  % Primary timestamp of this object
        T;  % A struct whose fields contain other timestamps
    end
    properties
        V;  % Where user stores arbitrary data
    end
    properties(Access = private)
        t_ = NaN;
        T_ = struct();
    end
    
    methods
        function this = Event(t, T, V)
            % MSessionExplorer.Event constructs an instance of this class
            %
            %   obj = MSessionExplorer.Event()
            %   obj = MSessionExplorer.Event(t)
            %   obj = MSessionExplorer.Event(t, T)
            %   obj = MSessionExplorer.Event(t, T, V)
            %
            if nargin > 0
                this.t = t;
            end
            if nargin > 1
                this.T = T;
            end
            if nargin > 2
                this.V = V;
            end
        end
        
        % Setters and Getters
        function this = set.t(this, val)
            assert(isscalar(val) && isnumeric(val), 't must be a numeric scalar');
            this.t_ = val;
        end
        function val = get.t(this)
            val = this.t_;
        end
        
        function this = set.T(this, val)
            if isempty(val)
                val = struct();
            end
            assert(isstruct(val), 'T must be a struct or empty array but instead was %s', class(val));
            isNum = structfun(@isnumeric, val);
            assert(all(isNum), 'Fields of T must be numeric array');
            this.T_ = val;
        end
        function val = get.T(this)
            val = this.T_;
        end
        
        function this = SetTfield(this, fieldName, val)
            % Assign value to a field in each T property of an object array
            % 
            %   objs = objs.SetTfield(fieldName, val)
            % 
            % Inputs
            %   fieldName       Name of the field to assign.
            %   val             Values to assign. The number of element must match that of the objects. 
            % Output
            %   objs            Modified objects
            assert(numel(this) == numel(val), 'Numbers of element do not match');
            if iscell(val)
                for i = 1 : numel(this)
                    this(i).T_.(fieldName) = val{i};
                end
            else
                for i = 1 : numel(this)
                    this(i).T_.(fieldName) = val(i);
                end
            end
        end
        function val = GetTfield(this, fieldName, isCat)
            % Return value of a field in each T property of an object array
            % 
            %   val = objs.GetTfield(fieldName, isCat)
            % 
            % Inputs
            %   fieldName       Name of the field to return.
            %   isCat           By default, field values are returned in a cell array unless all values
            %                   are scalar or there is only one object. Setting isCat to true force 
            %                   concatenating the cell array, whereas false retains it. 
            % Output
            %   val             Returned values. 
            val = cell(size(this));
            isScalar = false(size(this));
            for i = 1 : numel(this)
                val{i} = this(i).T_.(fieldName);
                isScalar(i) = isscalar(val{i});
            end
            if nargin < 3
                isCat = all(isScalar(:)) || numel(val) == 1;
            end
            if isCat
                val = this.IDenest(val);
            end
        end
        
        function this = SetVfield(this, fieldName, val)
            % Assign value to a field in each V property of an object array
            % 
            %   objs = objs.SetTfield(fieldName, val)
            % 
            % Inputs
            %   fieldName       Name of the field to assign.
            %   val             Values to assign. The number of element must match that of the objects. 
            % Output
            %   objs            Modified objects
            assert(numel(this) == numel(val), 'Numbers of element do not match');
            if iscell(val)
                for i = 1 : numel(this)
                    this(i).V.(fieldName) = val{i};
                end
            else
                for i = 1 : numel(this)
                    this(i).V.(fieldName) = val(i);
                end
            end
        end
        function val = GetVfield(this, fieldName, isCat)
            % Return value of a field in each V property of an object array
            % 
            %   val = objs.GetTfield(fieldName, isCat)
            % 
            % Inputs
            %   fieldName       Name of the field to return.
            %   isCat           By default, field values are returned in a cell array unless all values
            %                   are scalar or there is only one object. Setting isCat to true force 
            %                   concatenating the cell array, whereas false retains it. 
            % Output
            %   val             Returned values. 
            val = cell(size(this));
            isScalar = false(size(this));
            for i = 1 : numel(this)
                val{i} = this(i).V.(fieldName);
                isScalar(i) = isscalar(val{i});
            end
            if nargin < 3
                isCat = all(isScalar(:));
            end
            if isCat
                val = this.IDenest(val);
            end
        end
        
        % Function Overloading
        function val = double(this)
            % Return t properties of an object array in double
            val = double(this.IGet_t());
        end
        function val = single(this)
            % Return t properties of an object array in single
            val = single(this.IGet_t());
        end
        function val = isnan(this)
            % Return whether or not t properties of an object array are NaN
            val = isnan(this.IGet_t());
        end
        function val = diff(this, varargin)
            % Differences and Approximate Derivatives
            % It uses the same inputs and outputs as MATLAB diff function
            val = diff(double(this), varargin{:});
        end
        function [this, ind] = sort(this, varargin)
            % Sorting elements in object array based on t property
            % It uses the same inputs and outputs as MATLAB sort function. 
            [~, ind] = sort(this.IGet_t(), varargin{:});
            this = this(ind);
        end
        function varargout = histcounts(this, varargin)
            % Compute histogram based on values in t property
            % It uses the same inputs and outputs as MATLAB histcounts function
            val = double(this);
            switch nargout
                case 1
                    varargout{1} = histcounts(val, varargin{:});
                case 2
                    [varargout{1}, varargout{2}] = histcounts(val, varargin{:});
                case 3
                    [varargout{1}, varargout{2}, varargout{3}] = histcounts(val, varargin{:});
                otherwise
                    error('Unexpected number of output arguments');
            end
        end
        
        % Operator Overloading
        function obj = plus(obj1, obj2)
            % + operator adds objects with numbers
            % It applies to t and all fields in T, and returns the objects
            if isa(obj2, 'double')
                obj = IOffsetTime(obj1, obj2);
            else
                obj = IOffsetTime(obj2, obj1);
            end
        end
        function obj = minus(obj1, obj2)
            % - operator subtracts objects with numbers
            % It applies to t and all fields in T, and returns the objects
            if isa(obj2, 'double')
                obj = IOffsetTime(obj1, -obj2);
            else
                obj = IOffsetTime(obj2, -obj1);
            end
        end
        function obj = times(obj1, obj2)
            % .* operator multiplies objects with numbers element-wise
            % It applies to t and all fields in T, and returns the objects
            if isa(obj2, 'double')
                obj = IScaleTime(obj1, obj2);
            else
                obj = IScaleTime(obj2, obj1);
            end
        end
        function obj = mtimes(obj1, obj2)
            % * operator multiplies objects with numbers element-wise
            % It applies to t and all fields in T, and returns the objects
            obj = times(obj1, obj2);
        end
        function val = lt(obj1, obj2)
            % * operator multiplies objects with numbers element-wise
            % It applies to t and all fields in T, and returns the objects
            val = double(obj1) < double(obj2);
        end
        function val = gt(obj1, obj2)
            % > operator compares t property of the objects with numbers
            val = double(obj1) > double(obj2);
        end
        function val = le(obj1, obj2)
            % <= operator compares t property of the objects with numbers
            val = double(obj1) <= double(obj2);
        end
        function val = ge(obj1, obj2)
            % >= operator compares t property of the objects with numbers
            val = double(obj1) >= double(obj2);
        end
        function val = eq(obj1, obj2)
            % == operator compares t property of the objects with numbers
            val = double(obj1) == double(obj2);
        end
        function val = ne(obj1, obj2)
            % ~= operator compares t property of the objects with numbers
            val = double(obj1) ~= double(obj2);
        end
    end
    
    methods(Access = private)
        % Utilities
        function val = IGet_t(this)
            val = zeros(size(this));
            for i = 1 : numel(this)
                val(i) = this(i).t_;
            end
        end
        function this = IOffsetTime(this, val)
            if isscalar(val)
                val = repmat(val, size(this));
            end
            assert(all(size(this) == size(val)), 'Size of the object and numeric array does not agree');
            for i = 1 : numel(this)
                this(i).t_ = this(i).t_ + val(i);
                this(i).T_ = structfun(@(x) x + val(i), this(i).T_, 'Uni', false);
            end
        end
        function this = IScaleTime(this, val)
            if isscalar(val)
                val = repmat(val, size(this));
            end
            assert(all(size(this) == size(val)), 'Size of the object and numeric array does not agree');
            for i = 1 : numel(this)
                this(i).t_ = this(i).t .* val(i);
                this(i).T_ = structfun(@(x) x .* val(i), this(i).T_, 'Uni', false);
            end
        end
        function val = IDenest(~, val)
            if isrow(val)
                val = cat(2, val{:});
            elseif iscolumn(val)
                val = cat(1, val{:});
            else
                val = cell2mat(val);
            end
        end
    end
end

