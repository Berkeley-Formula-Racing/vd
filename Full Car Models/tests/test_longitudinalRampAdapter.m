function tests = test_longitudinalRampAdapter
tests = functiontests(localfunctions);
end

function testPureLongitudinalStateAndDiagnostics(testCase)
[cars,~] = carConfig();
settings = struct("speeds",5,"verbose",false);
run = rampSpeed.runLongitudinalRamp(cars{1,2},settings, ...
    struct("id","accel","label","accel","carRole","acceleration"),struct());
row = run.points(1,:);
verifyEqual(testCase,row.steer_rad,0,"AbsTol",1e-12);
verifyEqual(testCase,row.lat_velocity_mps,0,"AbsTol",1e-12);
verifyEqual(testCase,row.yaw_rate_rps,0,"AbsTol",1e-12);
verifyLessThanOrEqual(testCase,abs(row.aLat_achieved_mps2),1e-9);
verifyTrue(testCase,row.pure_ay0);
verifyFalse(testCase,row.lateral_metrics_applicable);
verifyTrue(testCase,isnan(run.perSpeed.K_linear_rad_per_mps2(1)));
end

function testMaxLongAccelDiagnosticsPreserveLegacyOutputs(testCase)
[cars,~] = carConfig();
[legacyTable,legacyAccel,legacyGuess] = max_long_accel(5,cars{1,2});
[diagnosticTable,diagnosticAccel,diagnosticGuess,diagnostics] = ...
    max_long_accel(5,cars{1,2},legacyGuess,struct( ...
    "maxFunctionEvaluations",2000,"constraintTolerance",1e-2, ...
    "stepTolerance",1e-10,"display","off"));

verifyEqual(testCase,size(diagnosticTable),size(legacyTable));
verifyEqual(testCase,class(diagnosticTable),class(legacyTable));
verifyTrue(testCase,isfinite(legacyAccel));
verifyTrue(testCase,isfinite(diagnosticAccel));
verifyEqual(testCase,diagnosticGuess(3),5,"AbsTol",1e-12);
verifyEqual(testCase,diagnostics.state,diagnosticGuess,"AbsTol",1e-12);
verifyTrue(testCase,isfield(diagnostics,"c"));
verifyTrue(testCase,isfield(diagnostics,"ceq"));
verifyTrue(testCase,isfield(diagnostics,"metrics"));
verifyEqual(testCase,diagnostics.max_inequality_violation, ...
    max([diagnostics.c(:);0]),"AbsTol",1e-12);
verifyEqual(testCase,diagnostics.max_equality_residual, ...
    max(abs(diagnostics.ceq(:))),"AbsTol",1e-12);
verifyEqual(testCase,diagnostics.pure_ay0, ...
    abs(diagnostics.metrics.gLat) <= 1e-12);
end

function testCancellationLeavesInvalidRowsForUnsolvedSpeeds(testCase)
[cars,~] = carConfig();
settings = struct("speeds",[5 10],"verbose",false);
cancelRequested = false;
callbacks = struct("onProgress",@captureProgress, ...
    "isCancelled",@isCancelled);

run = rampSpeed.runLongitudinalRamp(cars{1,2},settings, ...
    struct("id","accel","label","accel","carRole","acceleration"),callbacks);

verifyEqual(testCase,run.status,"cancelled");
verifyEqual(testCase,run.perSpeed.speed_mps,settings.speeds(:), ...
    "AbsTol",1e-12);
verifyTrue(testCase,run.perSpeed.valid(1));
verifyFalse(testCase,any(run.perSpeed.valid(2:end)));
verifyEqual(testCase,run.perSpeed.status(2),"invalid");
verifyThat(testCase,run.perSpeed.reason(2), ...
    matlab.unittest.constraints.ContainsSubstring("cancel"));

    function captureProgress(event)
        if event.completedSpeeds >= 1
            cancelRequested = true;
        end
    end

    function value = isCancelled()
        value = cancelRequested;
    end
end

function testFailedSpeedRowsRemainInvalid(testCase)
[cars,~] = carConfig();
run = rampSpeed.runLongitudinalRamp(cars{1,2}, ...
    struct("speeds",[5 NaN],"verbose",false), ...
    struct("id","accel","label","accel","carRole","acceleration"),struct());

verifyTrue(testCase,run.perSpeed.valid(1));
verifyFalse(testCase,run.perSpeed.valid(2));
verifyEqual(testCase,run.perSpeed.status(2),"invalid");
verifyThat(testCase,run.perSpeed.reason(2), ...
    matlab.unittest.constraints.ContainsSubstring("failed"));
verifyGreaterThanOrEqual(testCase,numel(run.runMeta.speedErrors),1);
end
