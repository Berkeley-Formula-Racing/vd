function tests = test_doe_model_selection
tests = functiontests(localfunctions);
end

function testSelectsAndPredictsFromCachedMetrics(testCase)
% Break caught: dropping cached metrics or model-comparison metadata makes
% offline surrogate analysis unavailable to downstream prediction.
setup_paths;
rng(3);
n = 80;
X = table(rand(n,1),rand(n,1),'VariableNames',{'x1','x2'});
metricTable = table((1:n)',true(n,1),repmat("",n,1), ...
    45 + 3*sin(4*X.x1) + 2*X.x2.^2, ...
    4*ones(n,1),40*ones(n,1),100*ones(n,1),0.95*ones(n,1), ...
    'VariableNames',{'case_index','valid','error_message','t_autox', ...
    't_skid','t_accel','min_Fz_N','gg_coverage'});
carCell = cell(n,2); %#ok<NASGU>
designTable = X; %#ok<NASGU>

folder = tempname;
mkdir(folder);
cleaner = onCleanup(@() rmdir(folder,'s')); %#ok<NASGU>
resultPath = fullfile(folder,'DOE_results.mat');
save(resultPath,'carCell','designTable','metricTable');

A = doeAnalyze(resultPath,struct('responsesWanted',{{'t_autox'}}, ...
    'plots',strings(0,1),'visible','off'));

verifyTrue(testCase,A.settings.usedCachedMetrics);
verifyTrue(testCase,isfield(A.gpModels,'t_autox'));
verifyTrue(testCase,isfield(A.validation.t_autox,'quadratic'));
verifyTrue(testCase,isfield(A.validation.t_autox,'gp'));
verifyTrue(testCase,any(A.preferredModel.t_autox == ["quadratic" "gp"]));
verifyEqual(testCase,A.predictorBounds.parameter,["x1";"x2"]);
verifyEqual(testCase,A.predictorBounds.lower,[min(X.x1);min(X.x2)],'AbsTol',1e-12);
verifyEqual(testCase,A.predictorBounds.upper,[max(X.x1);max(X.x2)],'AbsTol',1e-12);

[y,sd] = predictDOEModel(A,'t_autox',X(1:5,:));
verifySize(testCase,y,[5 1]);
verifySize(testCase,sd,[5 1]);
verifyTrue(testCase,all(isfinite(y)));
if A.preferredModel.t_autox == "gp"
    verifyTrue(testCase,all(isfinite(sd) & sd >= 0));
else
    verifyTrue(testCase,all(isnan(sd)));
end
end

function testRejectsMismatchedMetricCache(testCase)
% Break caught: accepting a stale cache silently assigns metrics to the
% wrong designs instead of using the legacy extraction path.
setup_paths;
modelRoot = fileparts(which('DOE_Fitting'));
paths = [fullfile(modelRoot,"DOE_results.mat"), ...
    fullfile(fileparts(modelRoot),"DOE_results.mat")];
fixturePath = paths(find(isfile(paths),1));
assert(~isempty(fixturePath),'DOE_results.mat fixture is required for this test.');
S = load(fixturePath,'carCell','designTable');
carCell = S.carCell(1:6,:); %#ok<NASGU>
designTable = S.designTable(1:6,:); %#ok<NASGU>
extractedMetrics = doeMetrics(carCell);
metricTable = extractedMetrics(1:5,:); %#ok<NASGU>

folder = tempname;
mkdir(folder);
cleaner = onCleanup(@() rmdir(folder,'s')); %#ok<NASGU>
resultPath = fullfile(folder,'DOE_results.mat');
save(resultPath,'carCell','designTable','metricTable');

A = doeAnalyze(resultPath,struct('responsesWanted',{{'t_autox'}}, ...
    'predictors',{{'mass','R_sf'}},'plots',strings(0,1), ...
    'visible','off','fitGaussianProcesses',false));

verifyFalse(testCase,A.settings.usedCachedMetrics);
verifyEqual(testCase,A.metricTable.t_autox,extractedMetrics.t_autox,'AbsTol',1e-12);
end

function testDefaultAnalysisIncludesAvailableRampResponses(testCase)
% Break caught: default offline analysis silently omits cached ramp metrics.
setup_paths;
[carCell,designTable,metricTable] = smallCachedData(true); %#ok<ASGLU>
folder = tempname;
mkdir(folder);
cleaner = onCleanup(@() rmdir(folder,'s')); %#ok<NASGU>
resultPath = fullfile(folder,'DOE_results.mat');
save(resultPath,'carCell','designTable','metricTable');

A = doeAnalyze(resultPath,struct('plots',strings(0,1),'visible','off', ...
    'fitGaussianProcesses',false));

expected = {'t_autox','understeer_gradient_10_deg_per_g', ...
    'understeer_gradient_25_deg_per_g','rebalance_speed_mps'};
verifyEqual(testCase,A.settings.responses,expected);
end

function testDefaultAnalysisSkipsUnavailableRampResponses(testCase)
% Break caught: requesting optional ramps makes legacy result files fail.
setup_paths;
[carCell,designTable,metricTable] = smallCachedData(false); %#ok<ASGLU>
folder = tempname;
mkdir(folder);
cleaner = onCleanup(@() rmdir(folder,'s')); %#ok<NASGU>
resultPath = fullfile(folder,'DOE_results.mat');
save(resultPath,'carCell','designTable','metricTable');

A = doeAnalyze(resultPath,struct('plots',strings(0,1),'visible','off', ...
    'fitGaussianProcesses',false));

verifyEqual(testCase,A.settings.responses,{'t_autox'});
end

function testLoadsCachedDataFromCheckpointContainer(testCase)
% Break caught: an offline checkpoint cannot be analysed without rerunning.
setup_paths;
[carCell,designTable,metricTable] = smallCachedData(false);
checkpoint = struct('carCell',{carCell},'designTable',designTable, ...
    'metricTable',metricTable); %#ok<NASGU>
folder = tempname;
mkdir(folder);
cleaner = onCleanup(@() rmdir(folder,'s')); %#ok<NASGU>
checkpointPath = fullfile(folder,'DOE_checkpoint.mat');
save(checkpointPath,'checkpoint');

A = doeAnalyze(checkpointPath,struct('responsesWanted',{{'t_autox'}}, ...
    'plots',strings(0,1),'visible','off','fitGaussianProcesses',false));

verifyTrue(testCase,A.settings.usedCachedMetrics);
verifyEqual(testCase,A.metricTable,metricTable);
verifyEqual(testCase,A.cleanTable.x1,designTable.x1,'AbsTol',1e-12);
end

function testLoadsResultsFromOutputDirectory(testCase)
% Break caught: offline analysis misses the default DOE_output result file.
setup_paths;
[carCell,designTable,metricTable] = smallCachedData(false); %#ok<ASGLU>
folder = tempname;
outputDirectory = fullfile(folder,'DOE_output');
mkdir(outputDirectory);
cleaner = onCleanup(@() rmdir(folder,'s')); %#ok<NASGU>
save(fullfile(outputDirectory,'DOE_results.mat'), ...
    'carCell','designTable','metricTable');

A = doeAnalyze(outputDirectory,struct('responsesWanted',{{'t_autox'}}, ...
    'plots',strings(0,1),'visible','off','fitGaussianProcesses',false));

verifyTrue(testCase,A.settings.usedCachedMetrics);
verifyEqual(testCase,A.cleanTable.x1,designTable.x1,'AbsTol',1e-12);
end

function testLoadsCheckpointFromOutputDirectory(testCase)
% Break caught: offline analysis misses the default DOE_output checkpoint.
setup_paths;
[carCell,designTable,metricTable] = smallCachedData(false);
checkpoint = struct('carCell',{carCell},'designTable',designTable, ...
    'metricTable',metricTable); %#ok<NASGU>
folder = tempname;
outputDirectory = fullfile(folder,'DOE_output');
mkdir(outputDirectory);
cleaner = onCleanup(@() rmdir(folder,'s')); %#ok<NASGU>
save(fullfile(outputDirectory,'DOE_checkpoint.mat'),'checkpoint');

A = doeAnalyze(outputDirectory,struct('responsesWanted',{{'t_autox'}}, ...
    'plots',strings(0,1),'visible','off','fitGaussianProcesses',false));

verifyTrue(testCase,A.settings.usedCachedMetrics);
verifyEqual(testCase,A.cleanTable.x1,designTable.x1,'AbsTol',1e-12);
end

function testDefaultAnalysisSkipsAllNanRampResponses(testCase)
% Break caught: default analysis attempts to fit unavailable cached ramps.
setup_paths;
[carCell,designTable,metricTable] = smallCachedData(true); %#ok<ASGLU>
metricTable.understeer_gradient_10_deg_per_g(:) = NaN;
metricTable.understeer_gradient_25_deg_per_g(:) = NaN;
metricTable.rebalance_speed_mps(:) = NaN;
folder = tempname;
mkdir(folder);
cleaner = onCleanup(@() rmdir(folder,'s')); %#ok<NASGU>
resultPath = fullfile(folder,'DOE_results.mat');
save(resultPath,'carCell','designTable','metricTable');

A = doeAnalyze(resultPath,struct('plots',strings(0,1),'visible','off', ...
    'fitGaussianProcesses',false));

verifyEqual(testCase,A.settings.responses,{'t_autox'});
end

function testExplicitAllNanResponseRaisesClearError(testCase)
% Break caught: explicitly requested unavailable response fails ambiguously.
setup_paths;
[carCell,designTable,metricTable] = smallCachedData(true); %#ok<ASGLU>
metricTable.understeer_gradient_10_deg_per_g(:) = NaN;
folder = tempname;
mkdir(folder);
cleaner = onCleanup(@() rmdir(folder,'s')); %#ok<NASGU>
resultPath = fullfile(folder,'DOE_results.mat');
save(resultPath,'carCell','designTable','metricTable');

verifyError(testCase,@() doeAnalyze(resultPath, ...
    struct('responsesWanted',{{'understeer_gradient_10_deg_per_g'}}, ...
    'plots',strings(0,1),'visible','off','fitGaussianProcesses',false)), ...
    'doeAnalyze:insufficientResponseData');
end

function [carCell,designTable,metricTable] = smallCachedData(includeRamps)
n = 6;
carCell = cell(n,2);
designTable = table(linspace(0,1,n)','VariableNames',{'x1'});
metricTable = table((1:n)',true(n,1),strings(n,1), ...
    45 + designTable.x1, ...
    'VariableNames',{'case_index','valid','error_message','t_autox'});
if includeRamps
    metricTable.understeer_gradient_10_deg_per_g = 1 + designTable.x1;
    metricTable.understeer_gradient_25_deg_per_g = 2 + designTable.x1;
    metricTable.rebalance_speed_mps = 20 + designTable.x1;
end
end
