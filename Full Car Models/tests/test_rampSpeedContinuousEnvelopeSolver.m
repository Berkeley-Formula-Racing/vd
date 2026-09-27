function tests = test_rampSpeedContinuousEnvelopeSolver
% Contract tests for the continuous-envelope ramp point solver.
tests = functiontests(localfunctions);
end

function testAllResolvedProfilesUseContinuousEnvelope(testCase)
profiles = rampSpeed.solverProfiles();
verifyTrue(testCase,all(string({profiles.powertrainModel}) == ...
    "continuousEnvelope"));
for profile = profiles
    resolved = rampSpeed.resolveSolverProfile(profile);
    verifyEqual(testCase,resolved.powertrainModel,"continuousEnvelope");
    metadata = rampSpeed.serializeSolverProfile(resolved);
    verifyEqual(testCase,string(metadata.powertrainModel),"continuousEnvelope");
end
end

function testPrimaryPointStoresEnvelopeAndNoGearAttempts(testCase)
[car,~] = carConfigBaseline();
profile = rampSpeed.resolveSolverProfile("fastPreview",struct( ...
    "maxFunctionEvaluations",350));
result = rampSpeed.solveLongitudinalPoint(car,20,[],profile,struct());

verifyTrue(testCase,ismember(result.status,["converged","near_feasible"]));
verifyEqual(testCase,result.diagnostics.powertrainModel,"continuousEnvelope");
verifyEmpty(testCase,result.diagnostics.gearAttempts);
verifyEqual(testCase,result.metrics.powertrain_model,"continuousEnvelope");
verifyTrue(testCase,isfinite(result.metrics.drivetrain_reduction));
end
