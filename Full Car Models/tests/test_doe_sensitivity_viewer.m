function tests = test_doe_sensitivity_viewer
tests = functiontests(localfunctions);
end

function testUpdatesViewerWhenResponseDropdownChanges(testCase)
% Break caught: selecting a different output does not update the ranked
% inputs or response plots in the consolidated sensitivity viewer.
setup_paths;
analysis = fixtureAnalysis();

f = doeSensitivityViewer(analysis,struct('visible','off','gridSize',9));
cleaner = onCleanup(@() close(f)); %#ok<NASGU>
S = f.UserData;

verifyEqual(testCase,string(S.responseDrop.Value),"t_autox");
verifyNumElements(testCase,findobj(S.rankAxes,'Type','bar'),1);
verifyNumElements(testCase,findobj(S.surfaceAxes,'Type','surface'),1);

S.responseDrop.Value = 't_accel';
S.refresh();
drawnow;
verifyTrue(testCase,contains(string(S.mainAxes.Title.String),"t_accel"));
verifyTrue(testCase,contains(string(S.qualityLabel.Text),"t_accel"));
end

function analysis = fixtureAnalysis()
rng(14);
T = table(rand(36,1),rand(36,1),rand(36,1), ...
    'VariableNames',{'x1','x2','x3'});
autox = 42 + 3*T.x1 - 2*T.x2 + T.x1.*T.x2;
accel = 4 - T.x1 + 0.8*T.x3;
analysis = struct();
analysis.cleanTable = [T table(autox,accel,'VariableNames',{'t_autox','t_accel'})];
analysis.settings = struct('predictors',{{'x1','x2','x3'}});
analysis.models = struct('t_autox',fitlm(T,autox), ...
    't_accel',fitlm(T,accel));
analysis.gpModels = struct('t_autox',[],'t_accel',[]);
analysis.preferredModel = struct('t_autox',"quadratic",'t_accel',"quadratic");
analysis.sobol = struct( ...
    't_autox',sobolTable(["x1";"x2";"x3"],[0.8;0.5;0.1]), ...
    't_accel',sobolTable(["x3";"x1";"x2"],[0.7;0.4;0.1]));
analysis.validation = struct( ...
    't_autox',validation(0.12,0.20), ...
    't_accel',validation(0.18,0.09));
end

function S = sobolTable(parameter,totalOrder)
S = table(parameter,totalOrder,'VariableNames',{'parameter','totalOrder'});
end

function V = validation(quadratic,GP)
V = struct('quadratic',struct('normalizedRMSE',quadratic), ...
    'gp',struct('normalizedRMSE',GP));
end
