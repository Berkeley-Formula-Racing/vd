function tests = test_doe_sobol
tests = functiontests(localfunctions);
end

function testReportsDominantAnalyticalFactorAndAttachesResults(testCase)
% Break caught: sampling in the wrong predictor scale or omitting the
% attachment makes a response dominated by x1 report incorrect Sobol ranks.
setup_paths;
rng(4);
n = 10;
X = table(rand(n,1),rand(n,1),'VariableNames',{'x1','x2'});
metricTable = table((1:n)',true(n,1),repmat("",n,1),X.x1 + 0.1*X.x2, ...
    'VariableNames',{'case_index','valid','error_message','response'});
carCell = cell(n,2); %#ok<NASGU>
designTable = X; %#ok<NASGU>

folder = tempname;
mkdir(folder);
cleaner = onCleanup(@() rmdir(folder,'s')); %#ok<NASGU>
resultPath = fullfile(folder,'DOE_results.mat');
save(resultPath,'carCell','designTable','metricTable');

A = doeAnalyze(resultPath,struct('responsesWanted',{{'response'}}, ...
    'plots',strings(0,1),'visible','off','fitGaussianProcesses',false, ...
    'sobolSamples',64,'sobolSeed',9,'sobolBootstrapSamples',5));
verifyTrue(testCase,isfield(A.sobol,'response'));

S = doeSobol(A,'response',20000,9);
x1 = S(S.parameter == "x1",:);
x2 = S(S.parameter == "x2",:);
verifyGreaterThan(testCase,x1.firstOrder,0.95);
verifyGreaterThan(testCase,x1.totalOrder,x2.totalOrder);
verifyLessThan(testCase,abs(sum(S.firstOrder)-1),0.08);
verifyLessThanOrEqual(testCase,x1.firstOrderLow,x1.firstOrder);
verifyGreaterThanOrEqual(testCase,x1.firstOrderHigh,x1.firstOrder);
verifyFalse(testCase,any(S.materialExcursion));
end

function testCatalogProducesOfflineSobolFigure(testCase)
% Break caught: omitting Sobol from the model-plot catalog leaves fitted
% offline analysis unable to produce its sensitivity figure.
setup_paths;
S = table(["x1";"x2"],[0.8;0.2],[0.9;0.3], ...
    [0.7;0.1],[0.9;0.3],[0.8;0.2],[1.0;0.4],[false;false], ...
    'VariableNames',{'parameter','firstOrder','totalOrder', ...
    'firstOrderLow','firstOrderHigh','totalOrderLow','totalOrderHigh', ...
    'materialExcursion'});
opts = struct('plots',"sobol",'visible','off', ...
    'sobol',struct('t_autox',S), ...
    'preferredModel',struct('t_autox',"quadratic"), ...
    'validation',struct('t_autox',struct('quadratic', ...
    struct('normalizedRMSE',0.02))));
figures = plotDOEMetrics(table(),table(),opts);
catalog = doePlotCatalog();
verifyTrue(testCase,any(catalog.id == "sobol"));
verifyEqual(testCase,catalog.producer(catalog.id == "sobol"),{'DOE_Fitting'});
verifyEqual(testCase,numel(figures),1);
verifyTrue(testCase,isgraphics(figures(1)));
verifyNotEmpty(testCase,findobj(figures(1),'Type','bar'));
close(figures(ishandle(figures)));
end

function testUsesPhysicalPredictorBounds(testCase)
% Break caught: mapping normalized samples directly into the surrogate makes
% the 100 kg mass span appear much less influential than the 1-unit CLA span.
A = physicalLinearAnalysis();
S = doeSobol(A,'response',4096,12);
mass = S(S.parameter == "mass",:);
cla = S(S.parameter == "cla",:);
verifyGreaterThan(testCase,mass.firstOrder,0.9);
verifyGreaterThan(testCase,mass.totalOrder,cla.totalOrder);
end

function testRestoresCallerRngState(testCase)
% Break caught: Sobol sampling changes the caller's random sequence, making
% downstream adaptive DOE selection non-repeatable.
A = physicalLinearAnalysis();
rng(61,'twister');
expected = rand(1,5);
rng(61,'twister');
doeSobol(A,'response',64,4);
verifyEqual(testCase,rand(1,5),expected,'AbsTol',0);
end

function testUsesSelectedGaussianProcess(testCase)
% Break caught: ignoring preferredModel.gp bypasses GP predictions and loses
% both their nonlinear response surface and uncertainty path.
[A,T] = nonlinearGpAnalysis();
verifyEqual(testCase,string(A.preferredModel.response),"gp");
[prediction,sd] = predictDOEModel(A,'response',T(1:3,:));
verifyTrue(testCase,all(isfinite(prediction)));
verifyTrue(testCase,all(isfinite(sd) & sd >= 0));
S = doeSobol(A,'response',256,8);
verifyEqual(testCase,height(S),2);
verifyTrue(testCase,all(isfinite(S.firstOrder) & isfinite(S.totalOrder)));
end

function testFlagsMaterialSobolExcursions(testCase)
% Break caught: clamping all estimates conceals a finite-sample Sobol result
% that needs interpretation rather than presenting a false in-range value.
A = physicalLinearAnalysis();
S = doeSobol(A,'response',2,1);
verifyTrue(testCase,any(S.materialExcursion));
verifyTrue(testCase,any(S.firstOrder < 0 | S.totalOrder > 1));
end

function A = physicalLinearAnalysis()
T = table([200;200;300;300],[1;2;1;2], ...
    [22;24;32;34],'VariableNames',{'mass','cla','response'});
mdl = fitlm(T,'response~mass+cla');
A = struct('models',struct('response',mdl),'gpModels',struct(), ...
    'preferredModel',struct('response',"quadratic"), ...
    'predictorBounds',table(["mass";"cla"],[200;1],[300;2], ...
    'VariableNames',{'parameter','lower','upper'}), ...
    'settings',struct('predictors',{{'mass','cla'}}, ...
    'sobolBootstrapSamples',5));
end

function [A,T] = nonlinearGpAnalysis()
[mass,cla] = ndgrid(linspace(220,260,4),linspace(1,3,4));
T = table(mass(:),cla(:),'VariableNames',{'mass','cla'});
y = sin(pi*(T.mass-220)/40) + 0.2*T.cla.^2;
gp = fitrgp(T{:,:},y,'KernelFunction','ardmatern32', ...
    'Standardize',true,'FitMethod','exact','PredictMethod','exact');
A = struct('models',struct(),'gpModels',struct('response',gp), ...
    'preferredModel',struct('response',"gp"), ...
    'predictorBounds',table(["mass";"cla"],[220;1],[260;3], ...
    'VariableNames',{'parameter','lower','upper'}), ...
    'settings',struct('predictors',{{'mass','cla'}}, ...
    'sobolBootstrapSamples',5));
end
