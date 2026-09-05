%% DOE_Fitting - metric extraction, response surfaces, and design plots
setup_paths

modelRoot = fileparts(which('DOE_Fitting'));
study = DOEStudyConfig();
resultSource = "";  % Optional: a DOE output folder or .mat result file.
resultPath = doeFindResultPath(modelRoot,study.output.directory,resultSource);

%% Offline analysis settings (this script does not run DOE simulations)
% Leave empty to fit every response with sufficient finite data. This lets
% a no-ramp screening study run without trying to fit NaN ramp responses.
% To force a subset, use {'t_autox','t_accel','t_skid'}.
responsesWanted = {};

plotsWanted = ["quality","sensitivity","main_effects","interactions", ...
    "pareto_events","pareto_energy","correlation","speed_grip", ...
    "balance","validation","sobol"];

fitGaussianProcesses = true;
runRampMetrics = false;  % false keeps analysis offline when cache is absent
outputDirectory = fileparts(resultPath);
analysisSavePath = fullfile(outputDirectory,"DOE_analysis.mat");
figureSavePath = fullfile(outputDirectory,"figures","doe");

opts = struct();
if ~isempty(responsesWanted), opts.responsesWanted = responsesWanted; end
opts.plots = plotsWanted;
opts.visible = 'on';
opts.runRampMetrics = runRampMetrics;
opts.fitGaussianProcesses = fitGaussianProcesses;
opts.savePath = analysisSavePath;

analysis = doeAnalyze(resultPath,opts);
interactionFigures = doeDefaultInteractionSurfaces(analysis, ...
    struct('visible',opts.visible));
sensitivityViewer = doeSensitivityViewer(analysis,struct('visible',opts.visible)); %#ok<NASGU>
allFigures = [analysis.figures interactionFigures];
savedFigures = saveFigures(figureSavePath,allFigures);
metricTable = analysis.metricTable;
doeTable = analysis.cleanTable;
models = analysis.models;
wantedGraphs = analysis.graphCatalog;

% Rotatable two-parameter response surface. Choose any two DOE predictors:
doeInteractionSurface(analysis,'t_autox','mass','wheelbase');

fprintf('\nDOE metrics: %d/%d valid cases, mean g-g coverage %.2f%%.\n', ...
    analysis.quality.valid,analysis.quality.total, ...
    100*analysis.quality.meanGGCoverage);
fprintf('Analysis saved to %s\n',analysisSavePath);
fprintf('Figures saved to %s\n',figureSavePath);
disp(wantedGraphs(:,{'id','title','purpose','producer','default_enabled'}));
