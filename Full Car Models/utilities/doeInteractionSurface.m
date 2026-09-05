function f = doeInteractionSurface(analysis,response,xParameter,yParameter,opts)
%DOEINTERACTIONSURFACE Rotatable two-parameter DOE surrogate surface.
%   f = doeInteractionSurface(analysis,response,xParameter,yParameter)
%   uses the selected surrogate in ANALYSIS, holds all other DOE inputs at
%   their valid-design mean, and returns a normal MATLAB figure. Rotate it
%   with the mouse to inspect the predicted interaction.

if nargin < 5, opts = struct(); end
response = scalarText(response,'response','doeInteractionSurface:badResponse');
xParameter = scalarText(xParameter,'xParameter', ...
    'doeInteractionSurface:badPredictor');
yParameter = scalarText(yParameter,'yParameter', ...
    'doeInteractionSurface:badPredictor');
if xParameter == yParameter
    error('doeInteractionSurface:duplicatePredictor', ...
        'xParameter and yParameter must be different.')
end
if ~isstruct(analysis) || ~isfield(analysis,'settings') || ...
        ~isfield(analysis.settings,'predictors') || ...
        ~isfield(analysis,'cleanTable') || ~istable(analysis.cleanTable)
    error('doeInteractionSurface:badAnalysis', ...
        'analysis must be the structure returned by doeAnalyze.')
end

predictors = string(analysis.settings.predictors);
if ~isfield(analysis,'preferredModel') || ...
        ~isfield(analysis.preferredModel,char(response))
    error('doeInteractionSurface:unknownResponse', ...
        'No fitted surrogate is available for response %s.',response)
end
if ~ismember(xParameter,predictors) || ~ismember(yParameter,predictors)
    error('doeInteractionSurface:unknownPredictor', ...
        'Both parameters must be DOE predictors. Available predictors: %s.', ...
        strjoin(cellstr(predictors),', '))
end

gridSize = getOption(opts,'gridSize',60);
if ~isnumeric(gridSize) || ~isscalar(gridSize) || ~isfinite(gridSize) || ...
        gridSize < 2 || gridSize ~= floor(gridSize)
    error('doeInteractionSurface:badGridSize', ...
        'opts.gridSize must be an integer of at least 2.')
end
visible = getOption(opts,'visible','on');
parent = getOption(opts,'parent',[]);

T = analysis.cleanTable;
x = linspace(min(T.(char(xParameter)),[],'omitnan'), ...
    max(T.(char(xParameter)),[],'omitnan'),gridSize);
y = linspace(min(T.(char(yParameter)),[],'omitnan'), ...
    max(T.(char(yParameter)),[],'omitnan'),gridSize);
[X,Y] = meshgrid(x,y);

Q = predictorMeans(T,predictors,numel(X));
Q.(char(xParameter)) = X(:);
Q.(char(yParameter)) = Y(:);
Z = reshape(predictDOEModel(analysis,response,Q),size(X));

ownsFigure = isempty(parent);
if ownsFigure
    f = figure('Name',sprintf('DOE interaction - %s',response), ...
        'Color','w','Visible',visible);
    ax = axes('Parent',f);
else
    if ~isgraphics(parent,'axes')
        error('doeInteractionSurface:badParent', ...
            'opts.parent must be an axes handle.')
    end
    ax = parent;
    f = ancestor(ax,'figure');
end
surf(ax,X,Y,Z,'EdgeColor','none','FaceAlpha',0.90);
xlabel(ax,xParameter,'Interpreter','none');
ylabel(ax,yParameter,'Interpreter','none');
zlabel(ax,response,'Interpreter','none');
title(ax,sprintf('%s: %s vs %s',response,xParameter,yParameter), ...
    'Interpreter','none');
colorbar(ax);
grid(ax,'on');
view(ax,-45,30);
if ownsFigure, rotate3d(f,'on'); end
setPlotFont(f);
end

function Q = predictorMeans(T,predictors,n)
Q = table();
for i = 1:numel(predictors)
    name = char(predictors(i));
    Q.(name) = repmat(mean(T.(name),'omitnan'),n,1);
end
end

function value = scalarText(value,name,identifier)
if ~(ischar(value) || (isstring(value) && isscalar(value))) || ...
        strlength(string(value)) == 0
    error(identifier,'%s must be a nonempty text scalar.',name)
end
value = string(value);
end

function value = getOption(opts,name,default)
if ~isstruct(opts) || ~isscalar(opts)
    error('doeInteractionSurface:badOptions','opts must be a scalar struct.')
end
if isfield(opts,name), value = opts.(name); else, value = default; end
end
