classdef MUtil
    %MUtil Summary of this class goes here
    %   Detailed explanation goes here
    
    methods(Static)
        function s = CombineStructs(varargin)
            % Combine multiple structs into one. Field names must be unique across structures. 
            %
            %   s = MUtil.CombineStructs(varargin)
            % 
            
            tbs = cellfun(@(x) struct2table(x, 'AsArray', true), varargin, 'Uni', false);
            s = table2struct(cat(2, tbs{:}));
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

