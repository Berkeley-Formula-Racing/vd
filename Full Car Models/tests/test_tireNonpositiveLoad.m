function tests = test_tireNonpositiveLoad
tests = functiontests(localfunctions);
end
function testLegacyTireReturnsZeroForZeroAndNegativeLoads(testCase)
% Break caught: abs(Fz) converts contact loss into a loaded tire and zero
% load reaches df_z divisions before the old post-hoc zero assignment.
testDir = fileparts(mfilename('fullpath'));
addpath(fileparts(testDir));
setup_paths;
[cars,~] = carConfig();
tire = cars{1,1}.tire;

verifyEqual(testCase,tire.F_x(0,0,-100,0),0,'AbsTol',1e-12);
verifyEqual(testCase,tire.F_y(0,0,-100,0),0,'AbsTol',1e-12);
verifyEqual(testCase,tire.F_x(0,0,0,0),0,'AbsTol',1e-12);
verifyEqual(testCase,tire.F_y(0,0,0,0),0,'AbsTol',1e-12);
end
