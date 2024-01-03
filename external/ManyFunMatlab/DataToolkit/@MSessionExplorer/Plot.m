function items = Plot(this, varargin)
% Plot any table variables for a quick view.
% Timeseries will be plotted as traces, event times as dots. Event values will be displayed in Command Window.
% 
%   items = Plot()
%   items = Plot(epoch)
%   items = Plot(epoch, tWin)
%   items = Plot(..., 'Items', items)
% 
% Inputs
%   se              A MSessionExplorer obejct.
%   epoch           Index of epoch to plot from.
%   tWin            A two-element vector indicating the time window to plot. The window can exceed the epoch.
%   'Items'         A struct where field names are table names, and the value of each field is a cell array of
%                   variable names for that table. It specifies the items to plot. For example,
%                   
%                   struct with fields:
%                           ni: {'mic'  'speaker1'}
%                     taskTime: {'trialOn'  'cue1'  'stim'  'cue3'  'prod'}
%                   
%                   If this struct is not provided, users will be prompted with list dialogs to select variables
%                   from each table interactively.
% Output
%   items           The struct that specifies the plots, which can be used as input for 'Items' to avoid prompts.
% 

p = inputParser();
p.addOptional('epoch', 1, @(x) isscalar(x) && x>0 && x<=this.numEpochs);
p.addOptional('tWin', [0 Inf], @(x) numel(x)==2 && x(2)>x(1));
p.addParameter('Items', [], @isstruct);
p.parse(varargin{:});
k = p.Results.epoch;
tWin = p.Results.tWin;
items = p.Results.Items;

% Select variables to plot
if isempty(items)
    items = struct;
    for i = 1 : numel(this.tableNames)
        tn = this.tableNames{i};
        tb = this.GetTable(tn);
        vn = tb.Properties.VariableNames;
        ind = listdlg('ListString', vn, ...
            'ListSize', [300 400], ...
            'PromptString', ['Select variables from the ''' tn ''' table to show'], ...
            'SelectionMode', 'multiple');
        items.(tn) = vn(ind);
    end
end

tRef = this.GetReferenceTime;
if isinf(tWin(2)) || isnan(tWin(2))
    if isempty(tRef)
        error('Reference times are not available. You must provide a complete tWin.');
    end
    tWin(2) = tRef(k+1)-tRef(k);
end

fn = fieldnames(items);
nPlots = 0;
for i = 1 : numel(fn)
    % Remove empty fields
    if isempty(fn{i})
        items = rmfield(items, fn{i});
    end
    
    % Tally the total number of axes needed
    m = strcmp(fn{i}, this.tableNames);
    if this.isEventTimesTable(m)
        nPlots = nPlots + 1;
    elseif this.isTimesSeriesTable(m)
        nPlots = nPlots + numel(items.(fn{i}));
    end
end

% Plot for each table
fn = fieldnames(items);
iPlot = 0;
for i = 1 : numel(fn)
    tn = fn{i};
    tbIdx = strcmp(this.tableNames, tn);
    vn = items.(tn);
    if isempty(vn)
        continue
    end
    
    if this.isTimesSeriesTable(tbIdx)
        subTb = this.SliceTimeSeries(tn, tWin, k, vn, 'Fill', 'bleed');
        for j = 1 : numel(vn)
            iPlot = iPlot+1;
            if iPlot > nPlots
                warning('Plot items exceed the total number of subplots');
                break
            end
            ax = subplot(nPlots, 1, iPlot); cla
            plot(subTb.time{1}, subTb.(vn{j}){1});
            xlim(tWin)
            ylabel(vn{j});
            MPlot.Axes(ax);
        end
        
    elseif this.isEventTimesTable(tbIdx)
        iPlot = iPlot+1;
        if iPlot > nPlots
            warning('Plot items exceed the total number of subplots');
            break
        end
        subTb = this.SliceEventTimes(tn, tWin, k, vn, 'Fill', 'bleed');
        ax = subplot(nPlots, 1, iPlot); cla
        for j = 1 : numel(vn)
            if iscell(subTb.(vn{j}))
                t = subTb.(vn{j}){1};
            else
                t = subTb.(vn{j});
            end
            if isnumeric(t) && isvector(t)
                y = repmat(j, size(t));
            else
                continue
            end
            plot(t, y, '.', 'MarkerSize', 6); hold on
        end
        axis ij
        xlim(tWin); ylim([0 j+1]);
        ax.YTick = 1 : j; ax.YTickLabel = vn;
        MPlot.Axes(ax);
        
    elseif this.isEventValuesTable(tbIdx)
        for j = 1 : numel(vn)
            fprintf("\nTable: '%s' | Column: '%s'\n", tn, vn{j});
            col = this.GetColumn(tn, vn{j});
            if iscell(col)
                disp(col{k});
            else
                disp(col(k));
            end
        end
        
    end
end

end