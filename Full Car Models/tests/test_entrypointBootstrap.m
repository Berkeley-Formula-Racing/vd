function tests = test_entrypointBootstrap
%TEST_ENTRYPOINTBOOTSTRAP Verify entrypoints resolve the project root.
tests = functiontests(localfunctions);
end

function testBootstrapSetsRepoRootAndPaths(testCase)
fullCarModelsDir = fileparts(fileparts(mfilename('fullpath')));
entrypointDir = fullfile(fullCarModelsDir, 'entrypoints');
bootstrapPath = fullfile(entrypointDir, 'bootstrap.m');
repoRoot = fileparts(fullCarModelsDir);

verifyTrue(testCase, isfile(bootstrapPath));

oldFolder = pwd;
cleanup = onCleanup(@() cd(oldFolder)); %#ok<NASGU>
run(bootstrapPath);

verifyEqual(testCase, pwd, repoRoot);
verifyEqual(testCase, fileparts(which('setup_paths')), fullCarModelsDir);
end

function testBaselineRunnerSelectsScalarBaseline(testCase)
fullCarModelsDir = fileparts(fileparts(mfilename('fullpath')));
runnerPath = fullfile(fullCarModelsDir, 'entrypoints', 'run_single_baseline_qss.m');
source = fileread(runnerPath);

verifyTrue(testCase, contains(source, '[allCars, eventParams, designTable] = carConfig();'));
verifyTrue(testCase, contains(source, 'static_front_ride_height_in'));
verifyTrue(testCase, contains(source, 'static_rear_ride_height_in'));
verifyTrue(testCase, contains(source, 'baselineIndex'));
verifyTrue(testCase, contains(source, 'makeGG(gg2(car, 0), car)'));
verifyFalse(testCase, contains(source, 'carConfig("Explicit", baselineTable)'));
end
