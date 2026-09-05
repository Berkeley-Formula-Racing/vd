function tests = test_run_lc0_nd_coupling
tests = functiontests(localfunctions);
end

function testWritesReproducibleCouplingFitArtifact(testCase)
% Break caught: a coupling result loses the pure-reference curves and exact
% TTC sources that define it.
addpath(fileparts(fileparts(mfilename('fullpath'))));
root = tempname; mkdir(root);
cleaner = onCleanup(@() rmdir(root,'s')); %#ok<NASGU>
longFile = fullfile(root,'long.mat'); latFile = fullfile(root,'lat.mat');
combinedFile = fullfile(root,'combined.mat');
writeRun(longFile,[-1;1],zeros(2,1),[-100;100],zeros(2,1));
writeRun(latFile,zeros(2,1),[-1;1],zeros(2,1),[-100;100]);
force = 100/sqrt(2);
writeRun(combinedFile,[-1;1;1],[0;0;1],[-100;100;force],[0;0;force]);
common = struct('load_N',100*4.4482216152605,'loadHalfWidth_N',1, ...
    'pressure_psi',12,'pressureHalfWidth_psi',.1, ...
    'camber_deg',0,'camberHalfWidth_deg',.1);
cfg = struct('donor',struct('file',combinedFile), ...
    'donorLateral',struct('file',latFile), ...
    'outputDirectory',fullfile(root,'results'), ...
    'couplingLongitudinalFit',merge(common,struct('maxAbsSlipAngle_deg',.1, ...
        'kappaGrid',[-1;1],'kappaHalfWidth',.01,'minSamplesPerBin',1)), ...
    'couplingLateralFit',merge(common,struct('maxAbsSlipRatio',.001, ...
        'alphaGrid',[-1;1],'alphaHalfWidth_deg',.01,'minSamplesPerBin',1)), ...
    'couplingFit',merge(common,struct('minAbsSlipRatio',.01, ...
        'minAbsSlipAngle_deg',.1,'pGrid',[1;2;3],'minPoints',1)));

result = run_lc0_nd_coupling(cfg,false);

verifyEqual(testCase,result.coupling.best_exponent,2,'AbsTol',1e-12);
verifyTrue(testCase,isfile(fullfile(cfg.outputDirectory,'lc0_nd_donor_coupling.mat')));
end

function writeRun(file,SL,SA,FX,FY)
FZ = -100*ones(numel(SL),1); %#ok<NASGU>
IA = zeros(numel(SL),1); P = 12*ones(numel(SL),1); V = 25*ones(numel(SL),1); %#ok<NASGU>
save(file,'SL','SA','FX','FY','FZ','IA','P','V');
end

function out = merge(a,b)
out = a;
names = fieldnames(b);
for i = 1:numel(names), out.(names{i}) = b.(names{i}); end
end
