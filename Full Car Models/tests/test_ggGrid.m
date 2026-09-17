function tests = test_ggGrid
tests = functiontests(localfunctions);
end

function setupOnce(~)
testDir = fileparts(mfilename('fullpath'));
addpath(fileparts(testDir));
setup_paths
end

function testFinalGridRetainsProductionResolution(testCase)
grid = ggGrid(30,struct('fastScreening',false));

verifyEqual(testCase,grid.velocityInterval,1);
verifyEqual(testCase,grid.lateralCount,20);
verifyEqual(testCase,grid.velocity,[5:30]);
end

function testFastScreeningUsesCoarserGrid(testCase)
grid = ggGrid(30,struct('fastScreening',true));

verifyEqual(testCase,grid.velocityInterval,2);
verifyEqual(testCase,grid.lateralCount,10);
verifyEqual(testCase,grid.velocity,[5:2:29 30]);
end
