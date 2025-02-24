classdef MUtil
    %MUtil Summary of this class goes here
    %   Detailed explanation goes here
    
    methods(Static)
        function x = ConvertNpyTypes(x, varargin)
            % Convert numpy to matlab datatypes
            % 
            %   x = ConvertNpyTypes(x)
            %   x = ConvertNpyTypes(x, ..., 'KeepUnsupported', true)
            %   x = ConvertNpyTypes(x, ..., 'Recursive', true)
            % 
            % Inputs
            %   'KeepUnsupported'       If set to false, unsupported datatypes are changed to [].
            %   'Recursive'             If set to true, any fields of the struct or dict x will also be converted.
            %                           The conversion will be performed recursively until variables are not struct 
            %                           or dict.
            % Output
            %   x                       Converted variable
            % 
            
            p = inputParser;
            p.addParameter('KeepUnsupported', true, @islogical);
            p.addParameter('Recursive', true, @islogical);
            p.parse(varargin{:});
            keepUnsupported = p.Results.KeepUnsupported;
            isRecursive = p.Results.Recursive;
            
            if isa(x, 'py.dict') || isstruct(x)
                x = struct(x);
                if ~isRecursive
                    return
                end
                n = fieldnames(x);
                for i = 1 : numel(n)
                    x.(n{i}) = MUtil.ConvertNpyTypes(x.(n{i}), 'Recursive', true);
                end
            elseif isa(x, 'py.int') || isa(x, 'py.numpy.float32') || isa(x, 'py.numpy.ndarray')
                x = double(x);
            elseif isa(x, 'py.str') || isa(x, 'py.list')
                x = string(x);
            elseif isa(x, 'py.NoneType')
                x = [];
            elseif startsWith(class(x), 'py,') && ~keepUnsupported
                x = [];
            end
        end
        
        function s = CombineStructs(varargin)
            % Combine multiple structs into one. Field names must be unique across structures. 
            %
            %   s = MUtil.CombineStructs(s1, s2, s3, ...)
            % 
            % Inputs
            %   s1, s2, s3...   Each input should be a scalar struct with mutually exclusive field names.
            % Output
            %   s               Combines struct.
            % 
            tbs = cellfun(@(x) struct2table(x, 'AsArray', true), varargin, 'Uni', false);
            s = table2struct(cat(2, tbs{:}));
        end
        
        function s = EvalFields(s)
            % Evaluate and convert interpretable fields (e.g. numeric arrays) in the struct
            % 
            %   s = MUtil.EvalFields(s)
            % 
            % Input
            %   s       A struct or a struct array.
            % Output
            %   s       A struct or a struct array with evaluated field values.
            % 
            if numel(s) > 1
                s = arrayfun(@MUtil.EvalFields, s);
                return
            end
            fn = fieldnames(s);
            for i = 1 : numel(fn)
                n = fn{i};
                try
                    s.(n) = eval(s.(n));
                catch
                end
            end
        end
        
        function hit = MatchAnyRegExp(str, expressions)
            % Match any of the regular expressions
            % 
            %   hit = MUtil.MatchAnyRegExp(str, expressions)
            %
            
            expressions = cellstr(expressions);
            expressions = expressions(:);
            hit = cellfun(@(x) any(regexp(str, x, 'once')), expressions);
            hit = any(hit);
        end
        
        function xlsTb = ReadXls(xlsPath, sheetId, varargin)
            % General-purposed function to read a spreadsheet from Excel file
            % 
            %   xlsTb = Tongue.Util.ReadXls(xlsPath, sheetId)
            % 
            % Inputs
            %   xlsPath         The path of an Excel file. If left empty, a file selection window 
            %                   will show up.
            %   sheetId         The index or name of a sheet for multi-sheet file. The default is
            %                   the first sheet. 
            % Output
            %   xlsTb           The output table. 
            % 
            
            if nargin < 2
                sheetId = 1;
            end
            
            % Browse for an Excel file
            if ~exist(xlsPath, 'file')
                xlsPath = MBrowse.File('', 'Choose an Excel spreadsheet', {'*.xlsx', '*.xls'});
            end
            
            % Find spreadsheet index by name
            if ischar(sheetId)
                [~, sheetNames] = xlsfinfo(xlsPath);
                sheetId = find(strcmpi(sheetId, sheetNames), 1);
                if isempty(sheetId)
                    error('The specified spreadsheet does not exist');
                end
            end
            
            % Load spreadsheet
            xlsTb = readtable(xlsPath, 'Sheet', sheetId, varargin{:});
        end
        
        function [txtStr, txtLines] = ReadTxt(svPath)
            % Read plain text of SatellitesViewer log to string(s). 
            % 
            %   [txtStr, txtLines] = MUtil.ReadTxt()
            %   [txtStr, txtLines] = MUtil.ReadTxt(svPath)
            % 
            % Input
            %   svPath              The path of a SatellitesViewer log file. In fact, you can use this method
            %                       to read any text file where lines are delimited by newline return \n. If 
            %                       svPath is not specified, a file selection window will be prompted. 
            % Outputs
            %   txtStr              The entire text file in a single string (character array). 
            %   txtLines            A cell array where each element is a line of the text. 
            % 
            
            % Handles user input
            if nargin < 1 || isempty(svPath)
                [svName, svDir] = uigetfile('*.txt', 'Please select a SatellitesViewer log file');
                if ~svName
                    txtStr = '';
                    txtLines = {};
                    return;
                end
                svPath = fullfile(svDir, svName);
            end
            
            % Read SatellitesViewer log file to text
            try
                fid = fopen(svPath);
                txtStr = fread(fid, '*char')';
                fclose(fid);
            catch
                fclose(fid);
                error(['Error occured when reading ' svPath]);
            end
            
            % Split text string into lines
            txtLines = MUtil.StringToLines(txtStr);
        end
        
        function txtLines = StringToLines(txtStr)
            % Split a string into lines based on return characters. Empty lines are removed. 
            %
            %   txtLines = MUtil.StringToLines(txtStr)
            %
            
            txtLines = strsplit(txtStr, {'\n', '\r'})';
            isEmpty = cellfun(@isempty, txtLines);
            txtLines(isEmpty) = [];
        end
    end
end

