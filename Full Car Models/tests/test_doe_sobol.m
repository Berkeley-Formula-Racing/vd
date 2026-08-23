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
