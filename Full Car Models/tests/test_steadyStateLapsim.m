function tests = test_steadyStateLapsim
%TEST_STEADYSTATELAPSIM Verify the simple LC0 event-runner entry points.
tests = functiontests(localfunctions);
end

function testRootEntrypointUsesSimpleLc0Runner(testCase)
assertSimpleLc0Runner(testCase, 'SteadyStateLapsim.m');
end

function testPortableEntrypointUsesSimpleLc0Runner(testCase)
assertSimpleLc0Runner(testCase, fullfile('entrypoints', 'SteadyStateLapsim.m'));
end

function assertSimpleLc0Runner(testCase, relativePath)
fullCarModelsDir = fileparts(fileparts(mfilename('fullpath')));
source = fileread(fullfile(fullCarModelsDir, relativePath));
normalized = regexprep(source, '\s+', '');

verifyTrue(testCase, contains(normalized, ...
    '[carCell,eventParams]=carConfig("FullFactorial",[],"lc0_nd");'));
verifyTrue(testCase, contains(source, 'eventsToRun'));
verifyTrue(testCase, contains(source, 'makeGG(gg2('));
verifyTrue(testCase, contains(source, 'Events2('));
verifyTrue(testCase, contains(source, 'eventSim.(solverName)();'));
verifyFalse(testCase, contains(source, 'runAdaptiveDOE'));
verifyFalse(testCase, contains(source, 'DOEStudyConfig'));
end
