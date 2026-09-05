function tests = test_run_lc0_nd_donor
tests = functiontests(localfunctions);
end

function testWritesProvisionalDonorFitArtifact(testCase)
% Break caught: a donor curve can be inspected but not reproduced later
% with the exact source/fit options that created it.
addpath(fileparts(fileparts(mfilename('fullpath'))));
root = tempname;
mkdir(root);
cleaner = onCleanup(@() rmdir(root,'s')); %#ok<NASGU>

SA = zeros(3,1); SL = [-.1;0;.1]; %#ok<NASGU>
FX = [-100;0;100]; FY = zeros(3,1); FZ = -100*ones(3,1); %#ok<NASGU>
IA = zeros(3,1); P = 12*ones(3,1); V = 25*ones(3,1); %#ok<NASGU>
file = fullfile(root,'donor.mat');
save(file,'SA','SL','FX','FY','FZ','IA','P','V');
cfg = struct('donor',struct('file',file), ...
    'donorFit',struct('load_N',100*4.4482216152605,'loadHalfWidth_N',1, ...
    'pressure_psi',12,'pressureHalfWidth_psi',.1, ...
    'camber_deg',0,'camberHalfWidth_deg',.1, ...
    'maxAbsSlipAngle_deg',.1,'kappaGrid',[-.1;0;.1], ...
    'kappaHalfWidth',.001,'minSamplesPerBin',1), ...
    'outputDirectory',fullfile(root,'results'));

result = run_lc0_nd_donor(cfg,false);

verifyEqual(testCase,result.fit.peak_drive_mu,1,'AbsTol',1e-12);
verifyTrue(testCase,isfile(fullfile(cfg.outputDirectory,'lc0_nd_donor_longitudinal.mat')));
end
