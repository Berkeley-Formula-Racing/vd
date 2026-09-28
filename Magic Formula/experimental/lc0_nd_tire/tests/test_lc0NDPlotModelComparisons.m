function tests = test_lc0NDPlotModelComparisons
tests = functiontests(localfunctions);
end

function testWritesSmallComparisonArtifacts(testCase)
% Break caught: comparison plots omit the stated friction envelope or the
% numerical error summary needed to interpret the matched force cases.
addpath(fileparts(fileparts(mfilename('fullpath'))));
root = tempname;
mkdir(root);
cleaner = onCleanup(@() rmdir(root,'s')); %#ok<NASGU>

result = makeSmallResult();
artifacts = lc0NDPlotModelComparisons(result,root,struct('visible','off'));

verifyTrue(testCase,isfile(artifacts.pngFile));
verifyTrue(testCase,isfile(artifacts.figFile));
verifyTrue(testCase,isfile(artifacts.csvFile));
verifyEqual(testCase,height(artifacts.metrics),2);
verifyEqual(testCase,artifacts.summary.combined_count,4);
verifyTrue(testCase,any(artifacts.metrics.quantity == "Fx"));
verifyTrue(testCase,any(artifacts.metrics.quantity == "Fy"));
verifyEqual(testCase,artifacts.envelopes.classical_p,2);
nominal = artifacts.envelopes.scenarios( ...
    strcmp({artifacts.envelopes.scenarios.name},'nominal'));
verifyEqual(testCase,nominal.p,1.3,'AbsTol',1e-12);
verifyEqual(testCase,max(nominal.x),nominal.mu_x,'AbsTol',1e-12);
end

function testRejectsResultWithoutMatchedCases(testCase)
addpath(fileparts(fileparts(mfilename('fullpath'))));
root = tempname;
mkdir(root);
cleaner = onCleanup(@() rmdir(root,'s')); %#ok<NASGU>

verifyError(testCase,@() lc0NDPlotModelComparisons( ...
    struct('model',struct()),root,struct('visible','off')), ...
    'lc0NDPlotModelComparisons:badResult');
end

function testAcceptsOlderComparisonArtifact(testCase)
% Break caught: a saved comparison made before extrapolation diagnostics
% became columns cannot be plotted or inspected.
addpath(fileparts(fileparts(mfilename('fullpath'))));
root = tempname;
mkdir(root);
cleaner = onCleanup(@() rmdir(root,'s')); %#ok<NASGU>

result = makeSmallResult();
result.comparison = removevars(result.comparison,'experimental_extrapolated');
artifacts = lc0NDPlotModelComparisons(result,root,struct('visible','off'));

verifyEqual(testCase,artifacts.summary.combined_extrapolated_count,0);
verifyTrue(testCase,isfile(artifacts.pngFile));
end

function result = makeSmallResult()
targetLateral = struct('curve',table([-12;-6;0;6;12],[-1.1;-1.3;0;1.3;1.1], ...
    true(5,1),'VariableNames',{'slip_angle_deg','mu_y','is_qualified'}));
donorLateral = struct('curve',table([-12;-6;0;6;12],[-1;-1.15;0;1.15;1], ...
    true(5,1),'VariableNames',{'slip_angle_deg','mu_y','is_qualified'}));
donorLongitudinal = struct('curve',table([-0.16;-0.08;0;0.08;0.16], ...
    [-1.0;-1.1;0;1.1;1.0],true(5,1), ...
    'VariableNames',{'slip_ratio','mu_x','is_qualified'}));
model = lc0NDBuildModel(targetLateral,donorLongitudinal,donorLateral, ...
    struct('rhoMu',1,'rhoStiff',1,'couplingExponent',1.3));

caseType = [repmat("pure lateral",5,1); repmat("pure longitudinal",5,1); ...
    repmat("combined",4,1)];
alpha = [-12;-6;0;6;12; 0;0;0;0;0; -6;6;-6;6];
kappa = [zeros(5,1); -0.16;-0.08;0;0.08;0.16; -0.08;-0.08;0.08;0.08];
fz = 100*ones(14,1);
experimentalFx = [zeros(5,1); -100;-110;0;110;100; -50;-50;50;50];
experimentalFy = [-110;-130;0;130;110; zeros(5,1); -80;80;-80;80];
legacyFx = [zeros(5,1); -95;-105;0;105;95; -45;-55;45;55];
legacyFy = [-105;-120;0;120;105; zeros(5,1); -75;75;-75;75];
comparison = table(caseType,alpha,kappa,fz,zeros(14,1),12*ones(14,1), ...
    experimentalFx,experimentalFy,legacyFx,legacyFy,true(14,1), ...
    [zeros(10,1);0.9;0.9;0.9;0.9],false(14,1),zeros(14,1), ...
    'VariableNames',{'case_type','alpha_deg','slip_ratio','Fz_N', ...
    'camber_deg','pressure_psi','experimental_Fx_N','experimental_Fy_N', ...
    'legacy_Fx_N','legacy_Fy_N','experimental_supported', ...
    'experimental_utilization','experimental_extrapolated', ...
    'experimental_clamp_count'});
result = struct('model',model,'comparison',comparison, ...
    'config',struct('reference',struct('load_N',100,'pressure_psi',12, ...
    'camber_deg',0)));
end
