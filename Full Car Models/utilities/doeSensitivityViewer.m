function f = doeSensitivityViewer(analysis,opts)
%DOESENSITIVITYVIEWER Interactive response and input sensitivity viewer.
%   Opens one UI for DOE output selection, Sobol input ranking, a predicted
%   main effect, and a two-input surrogate interaction surface.

if nargin < 2, opts = struct(); end
validateAnalysis(analysis);
visible = getOption(opts,'visible','on');
gridSize = getOption(opts,'gridSize',45);
if ~isnumeric(gridSize) || ~isscalar(gridSize) || gridSize < 2 || ...
        gridSize ~= floor(gridSize)
    error('doeSensitivityViewer:badGridSize', ...
        'opts.gridSize must be an integer of at least 2.')
end

responses = string(fieldnames(analysis.preferredModel));
predictors = string(analysis.settings.predictors);
f = uifigure('Name','DOE sensitivity viewer','Color','w', ...
    'Position',[100 100 1450 850],'Visible',visible);
root = uigridlayout(f,[2 1],'RowHeight',{38,'1x'}, ...
    'Padding',[8 8 8 8],'RowSpacing',6);
controls = uigridlayout(root,[1 7],'ColumnWidth',{60,180,48,180,48,180,'1x'}, ...
    'Padding',[0 0 0 0],'ColumnSpacing',6);
uilabel(controls,'Text','output','HorizontalAlignment','right');
responseDrop = uidropdown(controls,'Items',cellstr(responses));
uilabel(controls,'Text','x input','HorizontalAlignment','right');
xDrop = uidropdown(controls,'Items',cellstr(predictors));
uilabel(controls,'Text','y input','HorizontalAlignment','right');
yDrop = uidropdown(controls,'Items',cellstr(predictors));
qualityLabel = uilabel(controls,'Text','','HorizontalAlignment','left');

plots = uigridlayout(root,[2 2],'ColumnWidth',{'0.72x','1.28x'}, ...
    'RowHeight',{'1x','1x'},'Padding',[0 0 0 0],'RowSpacing',8,'ColumnSpacing',8);
rankAxes = uiaxes(plots);
rankAxes.Layout.Row = [1 2]; rankAxes.Layout.Column = 1;
mainAxes = uiaxes(plots);
mainAxes.Layout.Row = 1; mainAxes.Layout.Column = 2;
surfaceAxes = uiaxes(plots);
surfaceAxes.Layout.Row = 2; surfaceAxes.Layout.Column = 2;

S = struct('analysis',analysis,'responses',responses,'predictors',predictors, ...
    'gridSize',gridSize,'responseDrop',responseDrop,'xDrop',xDrop,'yDrop',yDrop, ...
    'qualityLabel',qualityLabel,'rankAxes',rankAxes,'mainAxes',mainAxes, ...
    'surfaceAxes',surfaceAxes);
S.refresh = @() refresh(f);
f.UserData = S;
responseDrop.ValueChangedFcn = @(~,~) resetInputs(f);
xDrop.ValueChangedFcn = @(~,~) refresh(f);
yDrop.ValueChangedFcn = @(~,~) refresh(f);
resetInputs(f);
setPlotFont(f);
end

function resetInputs(f)
S = f.UserData;
[names,~] = ranking(S.analysis,string(S.responseDrop.Value),S.predictors);
S.xDrop.Value = char(names(1));
if numel(names) >= 2
    S.yDrop.Value = char(names(2));
else
    S.yDrop.Value = char(otherPredictor(S.predictors,names(1)));
end
refresh(f);
end

function refresh(f)
S = f.UserData;
response = string(S.responseDrop.Value);
xName = string(S.xDrop.Value);
yName = string(S.yDrop.Value);
if xName == yName
    yName = otherPredictor(S.predictors,xName);
    S.yDrop.Value = char(yName);
end
[names,values] = ranking(S.analysis,response,S.predictors);
drawRanking(S.rankAxes,names,values,response);
drawMainEffect(S.mainAxes,S.analysis,response,xName,S.gridSize);
drawSurface(S.surfaceAxes,S.analysis,response,xName,yName,S.gridSize);
S.qualityLabel.Text = qualityText(S.analysis,response);
end

function drawRanking(ax,names,values,response)
cla(ax);
barh(ax,values,'FaceColor',[0.16 0.42 0.72]);
ax.YTick = 1:numel(names);
ax.YTickLabel = cellstr(names);
ax.YDir = 'reverse';
ax.TickLabelInterpreter = 'none';
xlabel(ax,'Total-order Sobol sensitivity');
title(ax,sprintf('%s input ranking',response),'Interpreter','none');
grid(ax,'on');
end

function drawMainEffect(ax,analysis,response,xName,gridSize)
T = analysis.cleanTable;
x = linspace(min(T.(char(xName)),[],'omitnan'), ...
    max(T.(char(xName)),[],'omitnan'),gridSize)';
Q = meanPredictorTable(T,string(analysis.settings.predictors),numel(x));
Q.(char(xName)) = x;
y = predictDOEModel(analysis,response,Q);
cla(ax);
plot(ax,x,y,'LineWidth',2,'Color',[0.13 0.40 0.70]);
hold(ax,'on');
if ismember(char(response),T.Properties.VariableNames)
    scatter(ax,T.(char(xName)),T.(char(response)),18, ...
        [0.2 0.2 0.2],'filled','MarkerFaceAlpha',0.35);
    legend(ax,{'surrogate, others at mean','simulated samples'}, ...
        'Location','best');
end
xlabel(ax,xName,'Interpreter','none');
ylabel(ax,response,'Interpreter','none');
title(ax,sprintf('%s main effect',response),'Interpreter','none');
grid(ax,'on');
end

function drawSurface(ax,analysis,response,xName,yName,gridSize)
T = analysis.cleanTable;
x = linspace(min(T.(char(xName)),[],'omitnan'), ...
    max(T.(char(xName)),[],'omitnan'),gridSize);
y = linspace(min(T.(char(yName)),[],'omitnan'), ...
    max(T.(char(yName)),[],'omitnan'),gridSize);
[X,Y] = meshgrid(x,y);
Q = meanPredictorTable(T,string(analysis.settings.predictors),numel(X));
Q.(char(xName)) = X(:);
Q.(char(yName)) = Y(:);
Z = reshape(predictDOEModel(analysis,response,Q),size(X));
cla(ax);
surf(ax,X,Y,Z,'EdgeColor','none','FaceAlpha',0.92);
xlabel(ax,xName,'Interpreter','none');
ylabel(ax,yName,'Interpreter','none');
zlabel(ax,response,'Interpreter','none');
title(ax,sprintf('%s interaction',response),'Interpreter','none');
colorbar(ax); grid(ax,'on'); view(ax,-45,30);
end

function [names,values] = ranking(analysis,response,predictors)
names = predictors(:);
values = zeros(numel(names),1);
if isfield(analysis,'sobol') && isfield(analysis.sobol,char(response))
    sobol = analysis.sobol.(char(response));
    if istable(sobol) && all(ismember({'parameter','totalOrder'},sobol.Properties.VariableNames))
        for i = 1:numel(names)
            row = string(sobol.parameter) == names(i);
            if any(row), values(i) = sobol.totalOrder(find(row,1)); end
        end
    end
end
[values,order] = sort(values,'descend');
names = names(order);
end

function value = otherPredictor(predictors,current)
value = predictors(find(predictors ~= current,1));
if isempty(value), error('doeSensitivityViewer:needTwoPredictors', ...
        'The sensitivity viewer needs at least two DOE predictors.')
end
end

function Q = meanPredictorTable(T,predictors,n)
Q = table();
for i = 1:numel(predictors)
    name = char(predictors(i));
    Q.(name) = repmat(mean(T.(name),'omitnan'),n,1);
end
end

function text = qualityText(analysis,response)
preferred = string(analysis.preferredModel.(char(response)));
text = sprintf('%s | %s surrogate',response,preferred);
if isfield(analysis,'validation') && isfield(analysis.validation,char(response))
    V = analysis.validation.(char(response));
    if isfield(V,char(preferred)) && isfield(V.(char(preferred)),'normalizedRMSE')
        nrmse = V.(char(preferred)).normalizedRMSE;
        if isfinite(nrmse), text = sprintf('%s | normalized RMSE %.3f',text,nrmse); end
    end
end
end

function validateAnalysis(analysis)
if ~isstruct(analysis) || ~isfield(analysis,'cleanTable') || ...
        ~istable(analysis.cleanTable) || ~isfield(analysis,'settings') || ...
        ~isfield(analysis.settings,'predictors') || ...
        ~isfield(analysis,'preferredModel') || isempty(fieldnames(analysis.preferredModel))
    error('doeSensitivityViewer:badAnalysis', ...
        'analysis must be the structure returned by doeAnalyze.')
end
if numel(analysis.settings.predictors) < 2
    error('doeSensitivityViewer:needTwoPredictors', ...
        'The sensitivity viewer needs at least two DOE predictors.')
end
end

function value = getOption(opts,name,default)
if ~isstruct(opts) || ~isscalar(opts)
    error('doeSensitivityViewer:badOptions','opts must be a scalar struct.')
end
if isfield(opts,name), value = opts.(name); else, value = default; end
end
