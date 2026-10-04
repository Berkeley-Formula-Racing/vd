function tests = test_qssRuntimeReport
%TEST_QSSRUNTIMEREPORT Contract coverage for the runtime benchmark report.
tests = functiontests(localfunctions);
end

function setupOnce(~)
testDir = fileparts(mfilename('fullpath'));
addpath(fileparts(testDir));
setup_paths
end

function testScreeningRunReportsFinitePhaseAndSolverCounts(testCase)
% Break caught: removing solver-count reporting or returning an incomplete
% runtime report leaves later optimisation decisions without a baseline.
options = struct( ...
    'repeats',1, ...
    'events',strings(1,0), ...
    'ggOptions',struct('fastScreening',true,'continuation',false, ...
        'maxVelocity',5));

report = benchmarkQssRuntime("legacy",options);

verifyEqual(testCase,report.configuration.model,"legacy");
verifyTrue(testCase,report.configuration.ggOptions.fastScreening);
verifyEqual(testCase,report.configuration.ggOptions.maxVelocity,5);
verifySize(testCase,report.totalSeconds,[1 1]);
verifyTrue(testCase,isfinite(report.totalSeconds));
verifyEqual(testCase,report.counts.requested,30);
verifyGreaterThanOrEqual(testCase,report.counts.requested,report.counts.accepted);
verifyGreaterThan(testCase,report.solvers.functionEvaluations.total,0);
verifyEqual(testCase,report.phases.configurationSeconds + ...
    report.phases.ggSolveSeconds + report.phases.ggAssemblySeconds + ...
    report.phases.eventPreparationSeconds,report.totalSeconds,'AbsTol',1e-12);
end
