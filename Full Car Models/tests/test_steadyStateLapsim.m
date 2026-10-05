function tests = test_steadyStateLapsim
%TEST_STEADYSTATELAPSIM Verify the retained legacy event-runner source contract.
tests = functiontests(localfunctions);
end

function testPortableEntrypointUsesLegacyRunner(testCase)
fullCarModelsDir = fileparts(fileparts(mfilename('fullpath')));
sourcePath = fullfile(fullCarModelsDir, 'entrypoints', 'SteadyStateLapsim.m');
source = fileread(sourcePath);
normalized = regexprep(source, '\s+', '');

verifyTrue(testCase, contains(normalized, ...
    '[carCell,eventParams]=carConfig("FullFactorial",[],"legacy");'));
verifyTrue(testCase, contains(source, 'eventsToRun'));
verifyTrue(testCase, contains(source, 'makeGG(gg2('));
verifyTrue(testCase, contains(source, 'Events2('));
verifyTrue(testCase, contains(source, 'eventSim.(solverName)();'));
verifyFalse(testCase, contains(normalized, '"lc0_nd"'));
verifyFalse(testCase, contains(source, 'runAdaptiveDOE'));
verifyFalse(testCase, contains(source, 'DOEStudyConfig'));
end
