function tests = test_doe_fitting_helpers
tests = functiontests(localfunctions);
end

function testPrefersDOEOutputOverStaleLegacyResult(testCase)
% Break caught: DOE_Fitting silently analyses a root-level legacy file even
% when a newer study output is available.
root = tempname;
mkdir(root);
cleaner = onCleanup(@() rmdir(root,'s')); %#ok<NASGU>
legacy = fullfile(root,'DOE_results.mat');
save(legacy,'root');
output = fullfile(root,'DOE_output_new');
mkdir(output);
expected = fullfile(output,'DOE_results.mat');
save(expected,'output');

actual = doeFindResultPath(root,fullfile(root,'DOE_output'));
verifyEqual(testCase,string(actual),string(expected));
end

function testBuildsDefaultSurfaceForEveryFittedResponse(testCase)
% Break caught: pressing Run in DOE_Fitting omits the rotatable surrogate
% interaction surface even when a response and two inputs are available.
setup_paths;
rng(5);
T = table(rand(30,1),rand(30,1),'VariableNames',{'x1','x2'});
y = 40 + 2*T.x1 - T.x2;
model = fitlm(T,y);
accelModel = fitlm(T,4 - T.x1 + 0.5*T.x2);
analysis = struct();
analysis.cleanTable = T;
analysis.settings = struct('predictors',{{'x1','x2'}});
analysis.models = struct('t_autox',model,'t_accel',accelModel);
analysis.gpModels = struct('t_autox',[],'t_accel',[]);
analysis.preferredModel = struct('t_autox',"quadratic",'t_accel',"quadratic");
analysis.sobol = struct('t_autox',table(["x2";"x1"],[0.2;0.8], ...
    'VariableNames',{'parameter','totalOrder'}), ...
    't_accel',table(["x1";"x2"],[0.7;0.3], ...
    'VariableNames',{'parameter','totalOrder'}));

figs = doeDefaultInteractionSurfaces(analysis,struct('visible','off','gridSize',9));
cleaner = onCleanup(@() close(figs(isgraphics(figs)))); %#ok<NASGU>

verifyNumElements(testCase,figs,1);
verifyNumElements(testCase,findobj(figs,'Type','surface'),2);
end

function testBuildsSliderSliceViewerForEachQuadraticResponse(testCase)
% Break caught: the adaptive one-click analysis drops the old 2-D plotSlice
% slider view that lets users inspect a fitted response at chosen settings.
setup_paths;
rng(8);
T = table(rand(30,1),rand(30,1),rand(30,1), ...
    'VariableNames',{'x1','x2','x3'});
model = fitlm(T,40 + T.x1 - 2*T.x2 + 0.5*T.x3);
analysis = struct('models',struct('t_autox',model));

figs = doeSliceFigures(analysis,struct('visible','off'));
cleaner = onCleanup(@() close(figs(isgraphics(figs)))); %#ok<NASGU>

verifyNumElements(testCase,figs,1);
verifyNotEmpty(testCase,findall(figs,'Type','axes'));
end
