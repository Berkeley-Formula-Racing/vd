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
