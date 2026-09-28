function tests = test_rampSpeedSolverProfiles
tests = functiontests(localfunctions);
end

function testAccurateProfileIsTheDefault(testCase)
profile = rampSpeed.resolveSolverProfile();

verifyEqual(testCase,profile.id,"accurate");
verifyEqual(testCase,profile.solverOptions.maxFunctionEvaluations,2000);
verifyEqual(testCase,profile.solverOptions.constraintTolerance,1e-2);
verifyEqual(testCase,profile.solverOptions.stepTolerance,1e-10);
verifyEqual(testCase,profile.solverOptions.display,"off");
verifyEqual(testCase,profile.aeroMode,"coupled");
verifyFalse(testCase,profile.approximate);
end

function testFastPreviewUsesReducedLongitudinalSolverBudget(testCase)
profile = rampSpeed.resolveSolverProfile("fastPreview");

verifyEqual(testCase,profile.id,"fastPreview");
verifyEqual(testCase,profile.solverOptions.maxFunctionEvaluations,1000);
verifyEqual(testCase,profile.solverOptions.constraintTolerance,1e-2);
verifyEqual(testCase,profile.solverOptions.stepTolerance,1e-8);
verifyEqual(testCase,profile.aeroMode,"coupled");
verifyFalse(testCase,profile.approximate);
end

function testApproximatePreviewSelectsStaticAero(testCase)
profile = rampSpeed.resolveSolverProfile("approximateAeroPreview");

verifyEqual(testCase,profile.id,"approximateAeroPreview");
verifyEqual(testCase,profile.aeroMode,"static");
verifyTrue(testCase,profile.approximate);
end

function testSolverOptionsCanBeOverriddenWithoutLosingProfileDefaults(testCase)
profile = rampSpeed.resolveSolverProfile("fastPreview",struct( ...
    "maxFunctionEvaluations",800,"stepTolerance",5e-8));

verifyEqual(testCase,profile.solverOptions.maxFunctionEvaluations,800);
verifyEqual(testCase,profile.solverOptions.stepTolerance,5e-8);
verifyEqual(testCase,profile.solverOptions.constraintTolerance,1e-2);
end

function testUnknownProfileIsRejected(testCase)
verifyError(testCase,@()rampSpeed.resolveSolverProfile("turbo"), ...
    "rampSpeed:unknownSolverProfile");
end

function testInvalidSolverOverrideIsRejected(testCase)
verifyError(testCase,@()rampSpeed.resolveSolverProfile("accurate", ...
    struct("maxFunctionEvaluations",0)), ...
    "rampSpeed:invalidSolverProfile");
end

function testProfileValidationRejectsAeroModeMismatch(testCase)
profile = rampSpeed.resolveSolverProfile("accurate");
profile.aeroMode = "static";
[ok,issues] = rampSpeed.validateSolverProfile(profile);

verifyFalse(testCase,ok);
verifyNotEmpty(testCase,issues);
end

function testProfileMetadataIsJsonSerializable(testCase)
profile = rampSpeed.resolveSolverProfile("approximateAeroPreview");
metadata = rampSpeed.serializeSolverProfile(profile);
decoded = jsondecode(jsonencode(metadata));

verifyEqual(testCase,string(decoded.id),"approximateAeroPreview");
verifyEqual(testCase,string(decoded.aeroMode),"static");
verifyTrue(testCase,logical(decoded.approximate));
verifyEqual(testCase,decoded.solverOptions.maxFunctionEvaluations,1000);
end

function testSettingsResolverSelectsProfileAndRetainsOverrides(testCase)
settings = struct("solverProfile","fastPreview", ...
    "solverOptions",struct("maxFunctionEvaluations",800));
[profile,normalized] = ...
    rampSpeed.resolveSolverProfileFromSettings(settings);

verifyEqual(testCase,profile.id,"fastPreview");
verifyEqual(testCase,profile.solverOptions.maxFunctionEvaluations,800);
verifyEqual(testCase,normalized.solverProfile,"fastPreview");
verifyEqual(testCase,normalized.solverOptions.constraintTolerance,1e-2);
end

function testApproximateProfileChangesOnlyLocalCarAero(testCase)
[car,~] = carConfigBaseline();
approximateCar = rampSpeed.applySolverProfileToCar(car, ...
    rampSpeed.resolveSolverProfile("approximateAeroPreview"));

verifyTrue(testCase,car.rideHeightAero.enabled);
verifyFalse(testCase,approximateCar.rideHeightAero.enabled);
end
