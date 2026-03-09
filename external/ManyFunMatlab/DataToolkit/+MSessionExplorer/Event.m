classdef Event
    % MSessionExplorer.Event objects work just like numeric values, thus can undergo
    % numeric operations and be used in eventTimes table of MSessionExplorer objects. 
    % Unlike simple numbers, however, these objects can carry additional timestamps 
    % and attributes of the events, and keep all times consistent under operations.
    %
    %   obj = MSessionExplorer.Event()
    %   obj = MSessionExplorer.Event(t)
    %   obj = MSessionExplorer.Event(t, T)
    %   obj = MSessionExplorer.Event(t, T, V)
    %
    % Inputs
    %   t       One or a vector of event timestamps.
    %   T       One or a vector of structs containing associated event times.
    %           T must have the same number of elements as t.
    %   V       One or a vector of structs containing associated event values.
    %           V must have the same number of elements as t.
    % Output
    %   obj     One or a vector of MSessionExplorer.Event objects.
    % 
    
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
        function obj = Event(t, T, V)
            % See class documentation above for the use of the constructor
            
            % This allow for constructing object array
            if nargin == 0
                return
            end
            
            % Add values to object
            for i = numel(t) : -1 : 1
                obj(i,1).t = t(i);
                if nargin > 1 && numel(T) == numel(t)
                    obj(i,1).T = T(i);
                end
                if nargin > 2 && numel(V) == numel(t)
                    obj(i,1).V = V(i);
                end
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
            %   objs = objs.SetVfield(fieldName, val)
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
            %                   are scalar or there is only one object. Setting isCat to true forces 
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
        
        % Special Modifiers
        function obj = Morph(obj, mdl)
            % Allows a model to transform timestamps in t and T
            % 
            %   obj = Morph(obj, mdl)
            % 
            % Input
            %   mdl         Any object or function supporting the syntax "tNew = mdl(tOld)".
            % Output
            %   obj         Objects with transformed time.
            % 
            for i = 1 : numel(obj)
                obj(i).t = mdl(obj(i).t);
                obj(i).T = structfun(@(x) mdl(x), obj(i).T, 'Uni', false);
            end
        end
        
        % Utilities
        function mask = MaskTimestamps(obj, t, tPad)
            % Get a mask for a vector of timestamps where samples spanned by utterances are set to true
            % obj.T must contain tOn and tOff fields
            % 
            %   mask = MaskTimestamps(obj, t)
            %   mask = MaskTimestamps(obj, t, tPad)
            % 
            % Inputs
            %   t           A numeric vector of timestamps.
            %   tPad        The amount of time to pad at the onset and offset of each object. Positive 
            %               padding extends the boundaries, negative padding shortens them. Default tPad
            %               is zero and no padding is performed.
            %               1) If tPad is a scalar, the same amount is padded to both onsets and offsets
            %                  for all objects.
            %               2) If tPad is a 1-by-2 vector, the first element is used to pad onsets and the 
            %                  second element is used to pad offsets for all objects.
            %               3) If tPad is an n-by-1 or n-by-2 vector where n is the number of objects, 
            %                  each object will use the value(s) in the corresponding row of tPad.
            % Output
            %   mask        A logical vector with the same size as t.
            % 
            
            if nargin < 3
                tPad = 0;
            end
            if size(tPad,2) == 1
                tPad = repmat(tPad, [1 2]);
            end
            if size(tPad,1) == 1
                tPad = repmat(tPad, [numel(obj) 1]);
            end
            
            mask = false(size(t));
            for i = 1 : numel(obj)
                if isnan(obj(i))
                    continue
                end
                
                t1 = obj(i).T.tOn - tPad(i,1);
                t2 = obj(i).T.tOff + tPad(i,2);
                
                ind = t >= t1 & t < t2;
                mask(ind) = true;
            end
        end
        
        % Function Overloading
        function obj = round(obj, varargin)
            % Round timestamp variables including t and fields in T
            % It follows the same input and output convention as MATLAB's round function. 
            for i = 1 : numel(obj)
                obj(i).t = round(obj(i).t, varargin{:});
                obj(i).T = structfun(@(x) round(x, varargin{:}), obj(i).T, 'Uni', false);
            end
        end
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
            % It follows the same input and output convention as MATLAB's sort function. 
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
        % Private utilities
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
            assert(all(size(this) == size(val)), 'Size of the object and the numeric array does not agree');
            for i = 1 : numel(this)
                this(i).t_ = this(i).t_ + val(i);
                this(i).T_ = structfun(@(x) x + val(i), this(i).T_, 'Uni', false);
            end
        end
        function this = IScaleTime(this, val)
            if isscalar(val)
                val = repmat(val, size(this));
            end
            assert(all(size(this) == size(val)), 'Size of the object and the numeric array does not agree');
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

