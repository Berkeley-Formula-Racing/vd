function tests = test_lc0NDLoadFreeRolling
tests = functiontests(localfunctions);
end

function testLoadsTargetRunInSIUnits(testCase)
% Break caught: TTC USCS signs and units leak into the nondimensional model.
addpath(fileparts(fileparts(mfilename('fullpath'))));
folder = tempname;
mkdir(folder);
cleaner = onCleanup(@() rmdir(folder,'s')); %#ok<NASGU>

SA = [-2; 0; 2]; %#ok<NASGU>
SL = [0; 2e-4; -1e-4]; %#ok<NASGU>
FY = [-100; 0; 100]; %#ok<NASGU>
FZ = [-200; -200; -200]; %#ok<NASGU>
IA = [0; 0; 0]; %#ok<NASGU>
P = [12; 12; 12]; %#ok<NASGU>
V = [25; 25; 25]; %#ok<NASGU>
save(fullfile(folder,'A1965run15.mat'), ...
    'SA','SL','FY','FZ','IA','P','V');

cfg = struct('target',struct('runDirectory',folder,'runIds',15, ...
    'filePrefix','A1965run','maxAbsSlip',1e-3));
[data,manifest] = lc0NDLoadFreeRolling(cfg);

verifyEqual(testCase,data.Fy_N,[-100;0;100]*4.4482216152605,'AbsTol',1e-10);
verifyEqual(testCase,data.Fz_N,200*4.4482216152605*ones(3,1),'AbsTol',1e-10);
verifyEqual(testCase,data.speed_mps,25*0.44704*ones(3,1),'AbsTol',1e-12);
verifyEqual(testCase,data.runId,15*ones(3,1));
verifyEqual(testCase,manifest.run_id,15);
verifyEqual(testCase,manifest.n_samples,3);
end

function testRejectsNonFreeRollingData(testCase)
% Break caught: a drive/brake file is accidentally treated as target lateral data.
addpath(fileparts(fileparts(mfilename('fullpath'))));
folder = tempname;
mkdir(folder);
cleaner = onCleanup(@() rmdir(folder,'s')); %#ok<NASGU>

SA = 0; SL = 0.05; FY = 0; FZ = -200; IA = 0; P = 12; V = 25; %#ok<NASGU>
save(fullfile(folder,'A1965run15.mat'), ...
    'SA','SL','FY','FZ','IA','P','V');
cfg = struct('target',struct('runDirectory',folder,'runIds',15, ...
    'filePrefix','A1965run','maxAbsSlip',1e-3));

verifyError(testCase,@() lc0NDLoadFreeRolling(cfg), ...
    'lc0NDLoadFreeRolling:notFreeRolling');
end
