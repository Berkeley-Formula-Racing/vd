function tests = test_rampSpeedRequestContract
tests = functiontests(localfunctions);
end

function testRejectsNonpositiveAndNonfiniteSpeeds(testCase)
invalidSpeeds = {0, -1, NaN, Inf};
for index = 1:numel(invalidSpeeds)
    request = makeRequest("longitudinal", invalidSpeeds{index});
    verifyError(testCase, @() rampSpeed.validateRequest(request), ...
        'rampSpeed:invalidSpeedDomain');
end
end

function testRejectsDuplicateSpeeds(testCase)
request = makeRequest("longitudinal", [5 10 10]);
verifyError(testCase, @() rampSpeed.validateRequest(request), ...
    'rampSpeed:invalidSpeedDomain');
end

function testRejectsNonmonotonicSpeeds(testCase)
request = makeRequest("longitudinal", [5 20 10]);
verifyError(testCase, @() rampSpeed.validateRequest(request), ...
    'rampSpeed:invalidSpeedDomain');
end

function testAcceptsIncreasingPositiveSpeedsForBothRampTypes(testCase)
rampTypes = ["lateral", "longitudinal"];
for rampType = rampTypes
    request = makeRequest(rampType, [5; 10; 20]);
    verifyTrue(testCase, rampSpeed.validateRequest(request));
end
end

function testNormalizesTextFieldsAndRequestVectors(testCase)
request = makeRequest(' LONGITUDINAL ', [5 10 20]);
request.settings.speedPolicy = ' FIXED ';
request.execution = ' SERIAL ';
request.caseIds = ["case-B", "case-A"];
request.solverProfile = ' Fast ';

normalized = rampSpeed.normalizeRequest(request);

verifyEqual(testCase, normalized.rampType, "longitudinal");
verifyEqual(testCase, normalized.settings.speeds_mps, [5; 10; 20]);
verifyClass(testCase, normalized.settings.speeds_mps, 'double');
verifyEqual(testCase, normalized.settings.speedPolicy, "fixed");
verifyEqual(testCase, normalized.execution.mode, "serial");
verifyEqual(testCase, normalized.caseIds, ["case-B"; "case-A"]);
verifyEqual(testCase, normalized.solverProfile, "fast");
end

function testRejectsNonfixedSpeedPolicy(testCase)
request = makeRequest("lateral", [5; 10; 20]);
request.settings.speedPolicy = "adaptive";

verifyError(testCase, @() rampSpeed.validateRequest(request), ...
    'rampSpeed:invalidSpeedPolicy');
verifyError(testCase, @() rampSpeed.normalizeRequest(request), ...
    'rampSpeed:invalidSpeedPolicy');
end

function testRejectsEmptySolverProfileAtValidationBoundary(testCase)
request = makeRequest("longitudinal", 5);
request.solverProfile = "   ";

verifyError(testCase, @() rampSpeed.validateRequest(request), ...
    'rampSpeed:invalidSolverProfile');
end

function testMakesOneStableTaskForEachRequestedSpeed(testCase)
request = rampSpeed.normalizeRequest(makeRequest("lateral", [5; 10; 20]));
plan = rampSpeed.makeSpeedPlan(request);

verifyTrue(testCase, istable(plan.tasks));
verifyEqual(testCase, plan.tasks.Properties.VariableNames, ...
    {'speedIndex', 'speed_mps', 'origin', 'passIndex', 'status'});
verifyEqual(testCase, plan.tasks.speedIndex, [1; 2; 3]);
verifyEqual(testCase, plan.tasks.speed_mps, [5; 10; 20]);
verifyEqual(testCase, plan.tasks.origin, repmat("requested", 3, 1));
verifyEqual(testCase, plan.tasks.passIndex, ones(3, 1));
verifyEqual(testCase, plan.tasks.status, repmat("planned", 3, 1));
end

function testMakesOneSpeedResultAndUsesNaNForEmptyMetrics(testCase)
plan = rampSpeed.makeSpeedPlan(rampSpeed.normalizeRequest( ...
    makeRequest("longitudinal", [5; 10; 20])));
task = plan.tasks(2, :);
metrics = struct('lateral_g', [], 'longitudinal_g', 0);
diagnostics = struct('message', "solver converged");
attempts = struct('attemptIndex', 1);

result = rampSpeed.makeSpeedResult(task, "converged", ...
    metrics, diagnostics, attempts);

verifyEqual(testCase, result.speedIndex, 2);
verifyEqual(testCase, result.speed_mps, 10);
verifyEqual(testCase, result.origin, "requested");
verifyEqual(testCase, result.passIndex, 1);
verifyEqual(testCase, result.status, "converged");
verifyTrue(testCase, isnan(result.metrics.lateral_g));
verifyEqual(testCase, result.metrics.longitudinal_g, 0);
verifyEqual(testCase, result.diagnostics, diagnostics);
verifyEqual(testCase, result.attempts, attempts);
end

function testAcceptsEveryContractStatus(testCase)
plan = rampSpeed.makeSpeedPlan(rampSpeed.normalizeRequest( ...
    makeRequest("longitudinal", 5)));
task = plan.tasks(1, :);
statuses = ["planned", "running", "converged", "near_feasible", ...
    "infeasible", "solver_failed", "cancelled"];

for status = statuses
    result = rampSpeed.makeSpeedResult(task, status);
    verifyEqual(testCase, result.status, status);
end
end

function testRejectsUnknownResultStatus(testCase)
plan = rampSpeed.makeSpeedPlan(rampSpeed.normalizeRequest( ...
    makeRequest("longitudinal", 5)));
task = plan.tasks(1, :);

verifyError(testCase, @() rampSpeed.makeSpeedResult(task, "unknown"), ...
    'rampSpeed:invalidStatus');
end

function request = makeRequest(rampType, speeds)
request = struct();
request.rampType = rampType;
request.settings = struct('speeds_mps', speeds);
request.execution = "serial";
request.caseIds = strings(0, 1);
request.solverProfile = "baseline";
end
