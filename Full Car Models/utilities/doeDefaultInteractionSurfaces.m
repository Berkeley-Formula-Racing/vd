function figs = doeDefaultInteractionSurfaces(analysis,opts)
%DOEDEFAULTINTERACTIONSURFACES Make one rotatable surface per DOE response.
%   The two inputs with the greatest finite total-order Sobol values are
%   selected for each response. Responses without two predictors are skipped.

if nargin < 2, opts = struct(); end
if ~isstruct(analysis) || ~isfield(analysis,'preferredModel') || ...
        ~isfield(analysis,'settings') || ~isfield(analysis.settings,'predictors')
    error('doeDefaultInteractionSurfaces:badAnalysis', ...
        'analysis must be the structure returned by doeAnalyze.')
end

responses = string(fieldnames(analysis.preferredModel));
pairs = strings(0,2);
keptResponses = strings(0,1);
for i = 1:numel(responses)
    response = responses(i);
    pair = pickPair(analysis,response);
    if numel(pair) < 2, continue, end
    keptResponses(end+1,1) = response; %#ok<AGROW>
    pairs(end+1,:) = pair(1:2)'; %#ok<AGROW>
end
if isempty(keptResponses)
    figs = gobjects(0);
    return
end

visible = getOption(opts,'visible','on');
f = figure('Name','DOE interaction surfaces','Color','w','Visible',visible);
tl = tiledlayout(f,'flow','TileSpacing','compact','Padding','compact');
for i = 1:numel(keptResponses)
    response = keptResponses(i);
    ax = nexttile(tl);
    try
        surfaceOpts = opts;
        surfaceOpts.parent = ax;
        doeInteractionSurface(analysis,response,pairs(i,1),pairs(i,2),surfaceOpts);
    catch ME
        warning('doeDefaultInteractionSurfaces:skippedResponse', ...
            'Skipped %s interaction surface: %s',response,ME.message)
    end
end
title(tl,'DOE interaction surfaces: top two total-order Sobol inputs');
rotate3d(f,'on');
setPlotFont(f);
figs = f;
end

function pair = pickPair(analysis,response)
predictors = string(analysis.settings.predictors);
pair = strings(0,1);
if isfield(analysis,'sobol') && isfield(analysis.sobol,char(response))
    S = analysis.sobol.(char(response));
    if istable(S) && all(ismember({'parameter','totalOrder'},S.Properties.VariableNames))
        usable = isfinite(S.totalOrder) & ismember(string(S.parameter),predictors);
        S = S(usable,:);
        if height(S) >= 2
            S = sortrows(S,'totalOrder','descend');
            pair = string(S.parameter(1:2));
        end
    end
end
if numel(pair) < 2 && numel(predictors) >= 2
    pair = predictors(1:2);
end
end

function value = getOption(opts,name,default)
if ~isstruct(opts) || ~isscalar(opts)
    error('doeDefaultInteractionSurfaces:badOptions','opts must be a scalar struct.')
end
if isfield(opts,name), value = opts.(name); else, value = default; end
end
