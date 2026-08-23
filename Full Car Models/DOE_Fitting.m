%% DOE_Fitting - metric extraction, response surfaces, and design plots
setup_paths

modelRoot = fileparts(which('DOE_Fitting'));
resultPath = fullfile(modelRoot,"DOE_results.mat");
if ~isfile(resultPath)
    legacyPath = fullfile(fileparts(modelRoot),"DOE_results.mat");
    if isfile(legacyPath)
        resultPath = legacyPath;
    else
        error('DOE_Fitting:noResults','Run SteadyStateLapsim before DOE_Fitting.');
    end
end

%% Analysis settings
responsesWanted = { ...
    't_autox','t_accel','t_skid','total_work_kJ', ...
    'gLat_peak_g','gg_lat_10_g','gg_lat_20_g','gg_lat_30_g', ...
    'gg_accel_20_g','gg_brake_20_g','min_Fz_N', ...
    'understeer_proxy_10_deg','understeer_proxy_25_deg'};

plotsWanted = ["quality","sensitivity","main_effects", ...
    "pareto_events","pareto_energy","correlation","speed_grip", ...
    "balance","validation"];

runRampMetrics = false;  % true can be expensive
analysisSavePath = fullfile(modelRoot,"DOE_analysis.mat");
figureSavePath = fullfile(modelRoot,"figures","doe");

opts = struct();
opts.responsesWanted = responsesWanted;
opts.plots = plotsWanted;
opts.visible = 'on';
opts.runRampMetrics = runRampMetrics;
opts.savePath = analysisSavePath;

analysis = doeAnalyze(resultPath,opts);
savedFigures = saveFigures(figureSavePath,analysis.figures);
metricTable = analysis.metricTable;
doeTable = analysis.cleanTable;
models = analysis.models;
wantedGraphs = analysis.graphCatalog;

fprintf('\nDOE metrics: %d/%d valid cases, mean g-g coverage %.2f%%.\n', ...
    analysis.quality.valid,analysis.quality.total, ...
    100*analysis.quality.meanGGCoverage);
fprintf('Analysis saved to %s\n',analysisSavePath);
fprintf('Figures saved to %s\n',figureSavePath);
disp(wantedGraphs(:,{'id','title','purpose','producer','default_enabled'}));
