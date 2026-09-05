function tests = test_lc0NDLoadDonor
tests = functiontests(localfunctions);
end

function testConvertsCombinedDonorChannelsToSI(testCase)
% Break caught: donor FX/FZ signs or USCS units leak into a normalized fit.
addpath(fileparts(fileparts(mfilename('fullpath'))));
root = tempname;
mkdir(root);
cleaner = onCleanup(@() rmdir(root,'s')); %#ok<NASGU>

SA = [-1;0;1]; SL = [-.1;0;.1]; %#ok<NASGU>
FX = [-120;0;120]; FY = [-80;0;80]; FZ = -200*ones(3,1); %#ok<NASGU>
IA = zeros(3,1); P = 12*ones(3,1); V = 25*ones(3,1); %#ok<NASGU>
save(fullfile(root,'donor.mat'),'SA','SL','FX','FY','FZ','IA','P','V');
cfg = struct('file',fullfile(root,'donor.mat'));

donor = lc0NDLoadDonor(cfg);

verifyEqual(testCase,donor.Fx_N,FX*4.4482216152605,'AbsTol',1e-10);
verifyEqual(testCase,donor.Fz_N,200*4.4482216152605*ones(3,1),'AbsTol',1e-10);
verifyEqual(testCase,donor.slipRatio,SL,'AbsTol',1e-12);
verifyEqual(testCase,donor.slipAngle_deg,SA,'AbsTol',1e-12);
end
