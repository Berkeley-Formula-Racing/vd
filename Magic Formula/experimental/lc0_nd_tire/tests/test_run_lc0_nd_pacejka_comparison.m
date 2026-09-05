function tests = test_run_lc0_nd_pacejka_comparison
tests = functiontests(localfunctions);
end

function testWritesMatchedComparisonArtifact(testCase)
% Break caught: a force-comparison plot is produced without persisting the
% exact matched states and experimental model configuration behind it.
addpath(fileparts(fileparts(mfilename('fullpath'))));
root = tempname; mkdir(root);
cleaner = onCleanup(@() rmdir(root,'s')); %#ok<NASGU>
targetDir = fullfile(root,'target'); mkdir(targetDir);
writeRun(fullfile(targetDir,'A1965run15.mat'),zeros(3,1),[-1;0;1], ...
    zeros(3,1),[-200;0;200]);
donorFile = fullfile(root,'donor.mat');
writeRun(donorFile,[-1;0;1],zeros(3,1),[-100;0;100],zeros(3,1));
donorLatFile = fullfile(root,'donorLat.mat');
writeRun(donorLatFile,zeros(3,1),[-1;0;1],zeros(3,1),[-100;0;100]);
common = struct('load_N',100*4.4482216152605,'loadHalfWidth_N',1, ...
    'pressure_psi',12,'pressureHalfWidth_psi',.1, ...
    'camber_deg',0,'camberHalfWidth_deg',.1);
cfg = struct('target',struct('runDirectory',targetDir,'runIds',15, ...
    'filePrefix','A1965run','maxAbsSlip',.001), ...
    'targetForceFit',merge(common,struct('maxAbsSlipRatio',.001, ...
    'alphaGrid',[-1;0;1],'alphaHalfWidth_deg',.01,'minSamplesPerBin',1)), ...
    'donor',struct('file',donorFile), 'donorLateral',struct('file',donorLatFile), ...
    'donorFit',merge(common,struct('maxAbsSlipAngle_deg',.01, ...
    'kappaGrid',[-1;0;1],'kappaHalfWidth',.01,'minSamplesPerBin',1)), ...
    'couplingLateralFit',merge(common,struct('maxAbsSlipRatio',.001, ...
    'alphaGrid',[-1;0;1],'alphaHalfWidth_deg',.01,'minSamplesPerBin',1)), ...
    'reference',struct('load_N',100*4.4482216152605,'pressure_psi',12,'camber_deg',0), ...
    'calibration',struct('rhoMu',1,'rhoStiff',1,'couplingExponent',2), ...
    'outputDirectory',fullfile(root,'results'));

result = run_lc0_nd_pacejka_comparison(cfg,@legacy,false);

verifyGreaterThan(testCase,height(result.comparison),0);
verifyTrue(testCase,isfile(fullfile(cfg.outputDirectory,'lc0_nd_pacejka_comparison.mat')));
end

function writeRun(file,SL,SA,FX,FY)
FZ = -100*ones(numel(SA),1); %#ok<NASGU>
IA = zeros(numel(SA),1); P = 12*ones(numel(SA),1); V = 25*ones(numel(SA),1); %#ok<NASGU>
save(file,'SL','SA','FX','FY','FZ','IA','P','V');
end

function [fx,fy] = legacy(alpha,kappa,fz,~)
fx = kappa*fz;
fy = alpha*fz;
end

function out = merge(a,b)
out = a;
names = fieldnames(b);
for i = 1:numel(names), out.(names{i}) = b.(names{i}); end
end
