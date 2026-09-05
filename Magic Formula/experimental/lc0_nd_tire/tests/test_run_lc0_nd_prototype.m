function tests = test_run_lc0_nd_prototype
tests = functiontests(localfunctions);
end

function testWritesMeasuredSummaryArtifact(testCase)
% Break caught: the analysis runner fails before preserving its exact input
% configuration and measured summary for later donor-model calibration.
addpath(fileparts(fileparts(mfilename('fullpath'))));
root = tempname;
mkdir(root);
cleaner = onCleanup(@() rmdir(root,'s')); %#ok<NASGU>

SA = [-2;-1;0;1;2]; %#ok<NASGU>
SL = zeros(5,1); %#ok<NASGU>
FY = [-100;-80;0;80;100]; %#ok<NASGU>
FZ = -100*ones(5,1); %#ok<NASGU>
IA = zeros(5,1); %#ok<NASGU>
P = 12*ones(5,1); %#ok<NASGU>
V = 25*ones(5,1); %#ok<NASGU>
save(fullfile(root,'A1965run15.mat'),'SA','SL','FY','FZ','IA','P','V');

cfg = lc0NDConfig();
cfg.target.runDirectory = root;
cfg.target.runIds = 15;
cfg.outputDirectory = fullfile(root,'results');
cfg.summary = struct('loadCenters_N',100*4.4482216152605, ...
    'loadHalfWidth_N',1,'pressureCenters_psi',12, ...
    'pressureHalfWidth_psi',0.1,'camberCenters_deg',0, ...
    'camberHalfWidth_deg',0.1,'smallSlipWindow_deg',1);

result = run_lc0_nd_prototype(cfg,false);

verifyEqual(testCase,result.summary.n_samples,5);
verifyTrue(testCase,isfile(fullfile(cfg.outputDirectory,'lc0_nd_target_summary.mat')));
saved = load(fullfile(cfg.outputDirectory,'lc0_nd_target_summary.mat'));
verifyEqual(testCase,saved.summary.peak_mu_y,1,'AbsTol',1e-12);
end
