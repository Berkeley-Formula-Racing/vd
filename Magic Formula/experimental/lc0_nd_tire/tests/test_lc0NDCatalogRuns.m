function tests = test_lc0NDCatalogRuns
tests = functiontests(localfunctions);
end

function testClassifiesFreeRollingAndSlipSweepFiles(testCase)
% Break caught: a candidate donor is chosen from a directory label instead
% of its actual recorded slip-ratio coverage.
addpath(fileparts(fileparts(mfilename('fullpath'))));
root = tempname;
mkdir(root);
cleaner = onCleanup(@() rmdir(root,'s')); %#ok<NASGU>

SA = [-5;0;5]; FY = [-100;0;100]; FZ = -200*ones(3,1); %#ok<NASGU>
SL = zeros(3,1); %#ok<NASGU>
save(fullfile(root,'free.mat'),'SA','SL','FY','FZ');
SL = [-0.1;0;0.1]; %#ok<NASGU>
save(fullfile(root,'combined.mat'),'SA','SL','FY','FZ');

catalog = lc0NDCatalogRuns(root);

verifyEqual(testCase,height(catalog),2);
free = catalog(catalog.file == "free.mat",:);
combined = catalog(catalog.file == "combined.mat",:);
verifyTrue(testCase,free.is_free_rolling);
verifyFalse(testCase,free.has_combined_coverage);
verifyFalse(testCase,combined.is_free_rolling);
verifyTrue(testCase,combined.has_combined_coverage);
verifyEqual(testCase,combined.max_abs_slip_ratio,0.1,'AbsTol',1e-12);
end
