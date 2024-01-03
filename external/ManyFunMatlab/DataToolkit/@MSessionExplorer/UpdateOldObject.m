function se = UpdateOldObject(se)
% Update objects generated and/or saved before Feb, 2019
%
%   se = UpdateOldObject(se)
%
% Input
%   se      An array of MSessionExplorer objects or file paths
% Output
%   se      An array of updated MSessionExplorer objects or file paths saved to

if ischar(se) || iscellstr(se) || isstring(se)
    se = cellstr(se);
end

for i = 1 : numel(se)
    if iscellstr(se)
        % Load, update and save object
        fprintf('Load file\n  %s', se{i});
        s = load(se{i});
        structfun(@Update, s);
        fprintf('Save file\n');
        save(se{i}, '-struct', 's');
    else
        % Update objet
        Update(se(i));
    end
end

end


function Update(se)

if ~isa(se, 'MSessionExplorer')
    return;
end

for i = 1 : numel(se)
    if isempty(se(i).epochInd) && se(i).numEpochs > 0
        disp('This object was generated before 2/26/2019');
        disp('Please refactor originalTrialInd to epochInd and numTrials to numEpochs in your code');
        se(i).epochInd = se(i).originalTrialInd;
    else
        disp('This object is of the latest version');
    end
end

end