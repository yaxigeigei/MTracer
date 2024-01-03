classdef MBrowse
    % Handy functions for searching and accessing files and folders
    
    methods(Static)
        function [filePath, folderPath, bareName, ext] = File(defaultFolderPath, dialogTitle, filterSpec)
            % Browse to select a file and get path parts
            %
            %   [filePath, folderPath, bareName, ext] = MBrowse.File()
            %   [filePath, folderPath, bareName, ext] = MBrowse.File(defaultFolderPath)
            %   [filePath, folderPath, bareName, ext] = MBrowse.File(defaultFolderPath, dialogTitle)
            %   [filePath, folderPath, bareName, ext] = MBrowse.File(defaultFolderPath, dialogTitle, filterSpec)
            %
            
            % Handle user inputs
            if nargin < 3
                filterSpec = {'*.*'};
            end
            
            if nargin < 2
                dialogTitle = [];
            end
            
            if nargin < 1
                defaultFolderPath = [];
            end
            
            if isempty(dialogTitle)
                dialogTitle = 'Please select a file';
            end
            
            if isempty(defaultFolderPath)
                try
                    load('lastimeDir', 'defaultFolderPath');
                catch
                    defaultFolderPath = pwd;
                end
            end
            
            % Get file path info
            [fileName, folderPath] = uigetfile(filterSpec, dialogTitle, defaultFolderPath);
            
            % Check selection
            if folderPath
                % Decompose paths
                filePath = fullfile(folderPath, fileName);
                [~, bareName, ext] = fileparts(fileName);
                
                % Remember new folder path
                defaultFolderPath = folderPath;
                save('lastimeDir.mat', 'defaultFolderPath');
            else
                filePath = '';
                bareName = '';
                ext = '';
            end
        end
        
        function [filePath, folderPath, bareNames, exts] = Files(defaultFolderPath, dialogTitle, filterSpec)
            % Browse to select multiple files and get path parts
            %
            %   [filePath, folderPath, bareNames, exts] = MBrowse.Files()
            %   [filePath, folderPath, bareNames, exts] = MBrowse.Files(defaultFolderPath)
            %   [filePath, folderPath, bareNames, exts] = MBrowse.Files(defaultFolderPath, dialogTitle)
            %   [filePath, folderPath, bareNames, exts] = MBrowse.Files(defaultFolderPath, dialogTitle, filterSpec)
            %
            
            % Handle user inputs
            if nargin < 3
                filterSpec = {'*.*'};
            end
            
            if nargin < 2
                dialogTitle = [];
            end
            
            if nargin < 1
                defaultFolderPath = [];
            end
            
            if isempty(dialogTitle)
                dialogTitle = 'Please select files';
            end
            
            if isempty(defaultFolderPath)
                try
                    load('lastimeDir', 'defaultFolderPath');
                catch
                    defaultFolderPath = pwd;
                end
            end
            
            % Get file path info
            [fileName, folderPath] = uigetfile(filterSpec, dialogTitle, defaultFolderPath, 'MultiSelect', 'on');
            
            if folderPath
                % Decompose paths
                fileName = cellstr(fileName)';
                filePath = cellfun(@(x) fullfile(folderPath, x), fileName, 'UniformOutput', false);
                [~, bareNames, exts] = cellfun(@fileparts, fileName, 'UniformOutput', false);
                
                % Remember new folder path
                defaultFolderPath = folderPath;
                save('lastimeDir.mat', 'defaultFolderPath');
            else
                filePath = {};
                bareNames = {};
                exts = {};
            end
        end
        
        function folderPath = Folder(defaultFolderPath, dialogTitle)
            % Browse to select a folder and get path parts
            %
            %   folderPath = MBrowse.Folder()
            %   folderPath = MBrowse.Folder(defaultFolderPath)
            %   folderPath = MBrowse.Folder(defaultFolderPath, dialogTitle)
            %
            
            % Handle user inputs
            if nargin < 2
                dialogTitle = [];
            end
            
            if nargin < 1
                defaultFolderPath = [];
            end
            
            if isempty(dialogTitle)
                dialogTitle = 'Please select a folder';
            end
            
            if isempty(defaultFolderPath)
                try
                    load('lastimeDir', 'defaultFolderPath');
                catch
                    defaultFolderPath = pwd;
                end
            end
            
            % Get file path info
            folderPath = uigetdir(defaultFolderPath, dialogTitle);
            
            if folderPath
                % Remember new folder path
                defaultFolderPath = folderPath;
                save('lastimeDir.mat', 'defaultFolderPath');
            end
        end
        
        function dirTable = Dir2Table(varargin)
            % List folder contents. It is exactly the same function as dir(...) but output a table. 
            %
            %   dirTable = MBrowse.Dir2Table(fileDir, fileName)
            %
            
            dirStruct = dir(varargin{:});
            dirTable = struct2table(dirStruct);
            
            for i = 1 : width(dirTable)
                if ischar(dirTable.(i))
                    dirTable.(i) = {dirTable.(i)};
                end
            end
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
        
    end
    
end




