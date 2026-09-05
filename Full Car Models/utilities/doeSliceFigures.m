function figs = doeSliceFigures(analysis,opts)
%DOESLICEFIGURES Open legacy-style interactive 2-D plotSlice viewers.
%   Each stepwise quadratic model gets a MATLAB plotSlice figure. Its slider
%   controls set the held-constant values of the unplotted predictors.

if nargin < 2, opts = struct(); end
if ~isstruct(analysis) || ~isfield(analysis,'models') || ...
        ~isstruct(analysis.models)
    error('doeSliceFigures:badAnalysis', ...
        'analysis must be the structure returned by doeAnalyze.')
end
if ~isstruct(opts) || ~isscalar(opts)
    error('doeSliceFigures:badOptions','opts must be a scalar struct.')
end
if isfield(opts,'visible'), visible = opts.visible; else, visible = 'on'; end

responses = string(fieldnames(analysis.models));
figs = gobjects(0);
for i = 1:numel(responses)
    response = responses(i);
    mdl = analysis.models.(char(response));
    if isempty(mdl), continue, end
    try
        f = uifigure('Name',sprintf('DOE slices - %s',response), ...
            'Color','w','Visible',visible);
        plotSlice(f,mdl);
        setPlotFont(f);
        figs(end+1) = f; %#ok<AGROW>
    catch ME
        if exist('f','var') && isgraphics(f), close(f); end
        warning('doeSliceFigures:skippedResponse', ...
            'Skipped %s slice viewer: %s',response,ME.message)
    end
end
end
