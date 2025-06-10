%% 

isKeepName = true;
isMakeCopy = false;

[filePaths, folderPath, fileNames, fileExts] = MBrowse.Files([], 'Choose media files', ...
    {'*.jpg; *.jpeg; *.png; *.bmp; *.tif; *.tiff'});

fileExts = lower(fileExts);

%% 

tb = table(filePaths, fileNames, fileExts);

for i = 1 : numel(filePaths)
    
    % Find datetime
    try
        % For image file
        I = imfinfo(filePaths{i});
        tb.imgInfo{i} = I;
        
        if isfield(I, 'DigitalCamera') && isfield(I.DigitalCamera, 'DateTimeOriginal')
            % Best to use camera time
            C = I.DigitalCamera;
            tb.camInfo{i} = C;
            dtField = C.DateTimeOriginal;
            tb.datetimeSource{i} = 'DateTimeOriginal';
            
        elseif isfield(I, 'DateTime')
            % Second best to use DateTime field
            dtField = I.DateTime;
            fprintf('%s does not have camera info. Use DateTime field.\n', [fileNames{i} fileExts{i}]);
            tb.datetimeSource{i} = 'DateTime';
            
        else
            % Otherwise use date modified
            dtField = I.FileModDate;
            warning('%s did not record date taken. Use last modified date.', [fileNames{i} fileExts{i}]);
            tb.datetimeSource{i} = 'FileModDate';
        end
    catch
        % For non-image file
        I = dir(filePaths{i});
        tb.fileInfo{i} = I;
        
        dtField = I.date;
        warning('%s is not an image file. Use last modified date.', [fileNames{i} fileExts{i}]);
        tb.datetimeSource{i} = 'date';
    end
    
    % Figure out the new file name
    formats = {'yyyy:MM:dd HH:mm:ss', 'yyyy-MM-dd HH:mm:ss', 'dd-MMM-yyyy HH:mm:ss'};
    for k = 1 : numel(formats)
        try
            dt = datetime(dtField, 'InputFormat', formats{k});
            break;
        catch
        end
    end
    dtStr = datestr(dt, 'yyyy-mm-dd HH.MM.SS');
    
    oldName = [fileNames{i} fileExts{i}];
    if strncmp(oldName, dtStr, length(dtStr))
        fprintf('%s has the correct name\n', oldName)
        tb.isRenamed(i) = false;
        continue
    end
    if isKeepName
        newName = [dtStr ' ' oldName];
    else
        newName = [dtStr fileExts{i}];
    end
    tb.newNames{i} = newName;
    
    % Rename file
    newPath = fullfile(folderPath, newName);
    if isMakeCopy
        copyfile(filePaths{i}, newPath);
    else
        movefile(filePaths{i}, newPath);
    end
    tb.isRenamed(i) = true;
    
end
