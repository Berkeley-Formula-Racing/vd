function tests = test_doe_interaction_surface
tests = functiontests(localfunctions);
end

function testCreatesSelectedInteractionSurface(testCase)
% Break caught: the legacy rotatable interaction view is unavailable from
% the adaptive DOE analysis and users cannot inspect arbitrary factor pairs.
setup_paths;
analysis = fixtureAnalysis();

f = doeInteractionSurface(analysis,'t_autox','x1','x2', ...
    struct('visible','off','gridSize',12));
cleaner = onCleanup(@() close(f)); %#ok<NASGU>

surface = findobj(f,'Type','surface');
verifyNumElements(testCase,surface,1);
verifySize(testCase,surface.ZData,[12 12]);
ax = ancestor(surface,'axes');
verifyEqual(testCase,string(ax.XLabel.String),"x1");
verifyEqual(testCase,string(ax.YLabel.String),"x2");
verifyEqual(testCase,string(ax.ZLabel.String),"t_autox");
end

function testRejectsUnknownInteractionPredictor(testCase)
% Break caught: a misspelled factor gives an opaque table-assignment error.
setup_paths;
analysis = fixtureAnalysis();

verifyError(testCase,@() doeInteractionSurface(analysis,'t_autox','x1','bad', ...
    struct('visible','off')),'doeInteractionSurface:unknownPredictor');
end

function analysis = fixtureAnalysis()
rng(19);
n = 56;
designTable = table(rand(n,1),rand(n,1), ...
    'VariableNames',{'x1','x2'}); %#ok<NASGU>
metricTable = table((1:n)',true(n,1),strings(n,1), ...
    42 + 4*designTable.x1 - 2*designTable.x2 + 3*designTable.x1.*designTable.x2, ...
    4*ones(n,1),40*ones(n,1),100*ones(n,1),0.95*ones(n,1), ...
    'VariableNames',{'case_index','valid','error_message','t_autox', ...
    't_skid','t_accel','min_Fz_N','gg_coverage'}); %#ok<NASGU>
carCell = cell(n,2); %#ok<NASGU>
folder = tempname;
mkdir(folder);
cleaner = onCleanup(@() rmdir(folder,'s')); %#ok<NASGU>
resultPath = fullfile(folder,'DOE_results.mat');
save(resultPath,'carCell','designTable','metricTable');
analysis = doeAnalyze(resultPath,struct('responsesWanted',{{'t_autox'}}, ...
    'plots',strings(0,1),'visible','off','fitGaussianProcesses',false));
end
