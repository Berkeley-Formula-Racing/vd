function tests = test_doe_metrics
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
setup_paths;
modelRoot = fileparts(which('DOE_Fitting'));
paths = [fullfile(modelRoot,"DOE_results.mat"), ...
         fullfile(fileparts(modelRoot),"DOE_results.mat")];
resultPath = paths(find(isfile(paths),1));
assert(~isempty(resultPath),'DOE_results.mat fixture is required for this test.');
testCase.TestData.result = load(resultPath,'carCell','designTable');
end

function testExtractsNamedMetrics(testCase)
S = testCase.TestData.result;
M = doeMetrics(S.carCell(1:2,:));

verifyEqual(testCase,height(M),2);
required = {'case_index','valid','t_autox','t_accel','t_skid', ...
    'total_work_kJ','v_mean_mps','v_max_mps','gLat_peak_g', ...
    'gg_lat_10_g','gg_lat_20_g','gg_lat_30_g','gg_accel_10_g', ...
    'gg_brake_20_g','min_Fz_N','wheel_lift_fraction','gg_coverage', ...
    'mechanical_balance','aero_balance','understeer_proxy_10_deg', ...
    'understeer_proxy_25_deg'};
verifyTrue(testCase,all(ismember(required,M.Properties.VariableNames)));
verifyEqual(testCase,M.t_autox(1),S.carCell{1,1}.comp.times.autocross,'AbsTol',1e-10);
verifyGreaterThan(testCase,M.gg_coverage(1),0.9);
verifyTrue(testCase,all(isfinite(M{1,required(6:end)}),'all'));
end

function testFailedCaseIsRetainedAndMarkedInvalid(testCase)
S = testCase.TestData.result;
cells = S.carCell(1,:);
cells(2,:) = {[],[]};
M = doeMetrics(cells);

verifyEqual(testCase,height(M),2);
verifyFalse(testCase,M.valid(2));
verifyTrue(testCase,isnan(M.t_autox(2)));
end

function testGraphCatalogAndSelectablePlot(testCase)
S = testCase.TestData.result;
M = doeMetrics(S.carCell(1:4,:));
C = doePlotCatalog();
required = {'quality','sensitivity','main_effects','interactions', ...
    'pareto_events','pareto_energy','correlation','speed_grip', ...
    'balance','validation'};
verifyTrue(testCase,all(ismember(required,C.id)));

opts = struct('plots',"quality",'visible','off');
figs = plotDOEMetrics(S.designTable(1:4,:),M,opts);
verifyEqual(testCase,numel(figs),1);
verifyTrue(testCase,isgraphics(figs(1),'figure'));
close(figs);
end

function testAnalysisFitsSelectedMetrics(testCase)
S = testCase.TestData.result;
carCell = S.carCell(1:6,:); %#ok<NASGU>
designTable = S.designTable(1:6,:); %#ok<NASGU>
fixture = [tempname '.mat'];
cleanup = onCleanup(@() deleteIfPresent(fixture)); %#ok<NASGU>
save(fixture,'carCell','designTable');

opts = struct();
opts.responsesWanted = {'t_autox','gg_lat_20_g'};
opts.predictors = {'mass','R_sf'};
opts.plots = "quality";
opts.visible = 'off';
opts.savePath = "";
A = doeAnalyze(fixture,opts);

verifyEqual(testCase,height(A.metricTable),6);
verifyTrue(testCase,all(isfield(A.models,opts.responsesWanted)));
verifyTrue(testCase,all(ismember(doePlotCatalog().id,A.graphCatalog.id)));
verifyEqual(testCase,numel(A.figures),1);
close(A.figures);
end

function deleteIfPresent(path)
if isfile(path), delete(path); end
end
