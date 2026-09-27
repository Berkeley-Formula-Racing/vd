function tests = test_rampSpeedLongitudinalPointSolver
% Focused contract tests for the reduced pure-longitudinal point solver.
tests = functiontests(localfunctions);
end

function testPureLongitudinalStateInvariants(testCase)
[car,~] = carConfigBaseline();
profile = rampSpeed.resolveSolverProfile("fastPreview",struct( ...
    "maxFunctionEvaluations",400));
result = rampSpeed.solveLongitudinalPoint(car,20,[],profile,struct());

verifyTrue(testCase,ismember(result.status,["converged","near_feasible"]), ...
    "20 m/s should produce a valid pure-longitudinal point.");
verifyEqual(testCase,result.state(1),0,"AbsTol",1e-12);
verifyEqual(testCase,result.state(3),20,"AbsTol",1e-12);
verifyEqual(testCase,result.state(4:5),[0 0],"AbsTol",1e-12);
verifyEqual(testCase,result.state(6:7),[0 0],"AbsTol",1e-12);
verifyEqual(testCase,result.state(8),result.state(9),"AbsTol",1e-12);
verifyLessThanOrEqual(testCase,abs(result.metrics.gLat),1e-12);
verifyTrue(testCase,isfield(result.diagnostics,"reducedState"));
verifyEqual(testCase,result.diagnostics.reducedState, ...
    [result.state(2),result.state(8)],"AbsTol",1e-12);
end

function testExplicitGearAttemptsCoverTransitionRange(testCase)
[car,~] = carConfigBaseline();
profile = rampSpeed.resolveSolverProfile("fastPreview",struct( ...
    "maxFunctionEvaluations",350));
speeds = [12.5 15 17.5 20 22.5 25];
for speed = speeds
    result = rampSpeed.solveLongitudinalPoint(car,speed,[],profile,struct());
    attempts = result.diagnostics.gearAttempts;
    verifyFalse(testCase,isempty(attempts),"Every speed must retain gear attempts.");
    verifyTrue(testCase,any([attempts.gear] == 3), ...
        "Gear 3 must be explicitly represented at every transition speed.");
end
end

function testFixedGearOverrideDoesNotMutateCar(testCase)
[car,~] = carConfigBaseline();
state = rampSpeed.makePureLongState(car,20,1,0.02);
original = car;
[~,~,~,~,~,~,~,gear] = car.equations(state,[],struct("gearOverride",3));
verifyEqual(testCase,gear,3);
verifyEqual(testCase,car.powertrain.gears,original.powertrain.gears);
verifyEqual(testCase,car.R_sf,original.R_sf,"AbsTol",0);
end

function testCancellationBeforeFirstEvaluation(testCase)
[car,~] = carConfigBaseline();
profile = rampSpeed.resolveSolverProfile("fastPreview",struct( ...
    "maxFunctionEvaluations",100));
control = struct("shouldCancel",@() true);
result = rampSpeed.solveLongitudinalPoint(car,20,[],profile,control);

verifyEqual(testCase,result.status,"cancelled");
verifyTrue(testCase,isempty(fieldnames(result.metrics)) || ...
    allNumericMetricsAreNaN(result.metrics));
verifyEqual(testCase,result.diagnostics.evaluationCount,0);
end

function testReducedSolverAgreesWithLegacyAtRepresentativeSpeeds(testCase)
[car,~] = carConfigBaseline();
profile = rampSpeed.resolveSolverProfile("fastPreview",struct( ...
    "maxFunctionEvaluations",450));
speeds = [5 10 15 20 25];
for speed = speeds
    reduced = rampSpeed.solveLongitudinalPoint(car,speed,[],profile,struct());
    [~,legacyAccel,legacyState,legacyDiagnostics] = max_long_accel( ...
        speed,car,[],profile.solverOptions);
    if isfinite(legacyAccel) && legacyDiagnostics.exitflag > 0 && ...
            ismember(reduced.status,["converged","near_feasible"])
        verifyEqual(testCase,reduced.metrics.aLong_mps2,legacyAccel, ...
            "RelTol",0.05, ...
            sprintf("Reduced/legacy acceleration mismatch at %.1f m/s.",speed));
        verifyEqual(testCase,reduced.state(3),legacyState(3),"AbsTol",1e-12);
        verifyLessThanOrEqual(testCase, ...
            abs(reduced.state(8)-legacyState(8)),0.03, ...
            sprintf("Rear slip mismatch at %.1f m/s.",speed));
    end
end
end

function tf = allNumericMetricsAreNaN(value)
tf = true;
names = fieldnames(value);
for i = 1:numel(names)
    candidate = value.(names{i});
    if isnumeric(candidate) && ~isempty(candidate) && any(isfinite(candidate(:)))
        tf = false;
        return
    elseif isstruct(candidate) && isscalar(candidate) && ...
            ~allNumericMetricsAreNaN(candidate)
        tf = false;
        return
    end
end
end
