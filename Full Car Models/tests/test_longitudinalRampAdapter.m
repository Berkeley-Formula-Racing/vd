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

function testAppResidualToleranceDoesNotRejectLongitudinalRows(testCase)
[car,~] = carConfigBaseline();
settings = struct("speeds",[5 10 15 20],"mode","coast", ...
    "residualTolerance",1e-6,"verbose",false);
caseInfo = struct("id","baseline","label","baseline", ...
    "carRole","auto");

run = rampSpeed.runLongitudinalRamp(car,settings,caseInfo,struct());

verifyTrue(testCase,all(run.perSpeed.valid));
verifyEqual(testCase,run.perSpeed.status, ...
    repmat("complete",4,1));
verifyTrue(testCase,all(isfinite(run.perSpeed.aLong_mps2)));
end

function testDefaultLongitudinalSweepRecoversContinuationFailure(testCase)
[car,~] = carConfigBaseline();
settings = struct("speeds",(5:2.5:25).',"mode","coast", ...
    "nRamp",10,"nBisect",8,"residualTolerance",1e-6, ...
    "verbose",false);
caseInfo = struct("id","baseline","label","baseline", ...
    "carRole","auto");

run = rampSpeed.runLongitudinalRamp(car,settings,caseInfo,struct());

verifyTrue(testCase,all(run.perSpeed.valid), ...
    "Every default-grid longitudinal point should be recovered.");
verifyEqual(testCase,run.perSpeed.status, ...
    repmat("complete",numel(settings.speeds),1));
verifyTrue(testCase,all(isfinite(run.perSpeed.downforce_N)));
end

function testBorderlineSpeedUsesExplicitGearSweep(testCase)
[car,~] = carConfigBaseline();
settings = struct("speeds",23.75,"solverProfile","accurate", ...
    "speedGrid",struct("mode","fixed"),"verbose",false);
run = rampSpeed.runLongitudinalRamp(car,settings, ...
    struct("id","borderline","label","borderline","carRole","auto"),struct());
verifyTrue(testCase,run.perSpeed.valid, ...
    "The borderline speed should be recovered by the tight solver retry.");
verifyLessThanOrEqual(testCase,run.perSpeed.max_constraint_residual, ...
    1e-2);
verifyEmpty(testCase,run.runMeta.retrySpeeds_mps);
verifyGreaterThanOrEqual(testCase, ...
    numel(run.raw.diagnostics(1).attempts),1);
end


function testLongitudinalRunnerRecordsSelectedSolverProfile(testCase)
[car,~] = carConfigBaseline();
settings = struct("speeds",5,"solverProfile","fastPreview", ...
    "verbose",false);
caseInfo = struct("id","baseline","label","baseline","carRole","auto");

run = rampSpeed.runLongitudinalRamp(car,settings,caseInfo,struct());

verifyEqual(testCase,run.runMeta.solverProfile.id,"fastPreview");
verifyEqual(testCase,run.runMeta.solver.maxFunctionEvaluations,1000);
verifyEqual(testCase,run.settings.solverProfile,"fastPreview");
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

function testRunStatusDistinguishesMixedAndAllFailedSpeeds(testCase)
[cars,~] = carConfig();
caseInfo = struct("id","accel","label","accel","carRole","acceleration");

mixed = rampSpeed.runLongitudinalRamp(cars{1,2}, ...
    struct("speeds",[5 NaN],"verbose",false),caseInfo,struct());
allFailed = rampSpeed.runLongitudinalRamp(cars{1,2}, ...
    struct("speeds",NaN,"verbose",false),caseInfo,struct());

verifyEqual(testCase,mixed.status,"partial");
verifyTrue(testCase,mixed.perSpeed.valid(1));
verifyFalse(testCase,mixed.perSpeed.valid(2));
verifyEqual(testCase,allFailed.status,"failed");
verifyFalse(testCase,any(allFailed.perSpeed.valid));
end

function testNonconvergedSpeedsAreFailedDiagnostics(testCase)
[cars,~] = carConfig();
settings = struct("speeds",[5 10],"verbose",false, ...
    "solverOptions",struct("maxFunctionEvaluations",1, ...
    "constraintTolerance",1e-2,"stepTolerance",1e-10,"display","off"));
run = rampSpeed.runLongitudinalRamp(cars{1,2},settings, ...
    struct("id","limited","label","limited","carRole","acceleration"),struct());

verifyEqual(testCase,run.status,"failed");
verifyFalse(testCase,any(run.perSpeed.valid));
verifyFalse(testCase,run.raw.diagnostics(1).success);
verifyFalse(testCase,run.raw.diagnostics(2).success);
verifyGreaterThanOrEqual(testCase,numel(run.runMeta.speedErrors),2);
diagnostic = run.raw.diagnostics(1);
verifyTrue(testCase,isfinite(diagnostic.exitflag));
verifyTrue(testCase,~isempty(diagnostic.state));
verifyTrue(testCase,all(isfinite(diagnostic.state(:))));
verifyTrue(testCase,~isempty(diagnostic.c));
verifyTrue(testCase,~isempty(diagnostic.ceq));
verifyTrue(testCase,isstruct(diagnostic.metrics));
verifyTrue(testCase,isfinite(diagnostic.max_equality_residual));
verifyTrue(testCase,isfinite(diagnostic.max_inequality_violation));
verifyEqual(testCase,diagnostic.error_identifier,"rampSpeed:nonconverged");
verifyThat(testCase,diagnostic.error_message, ...
    matlab.unittest.constraints.ContainsSubstring("nonconverged solver result"));
verifyEqual(testCase,run.runMeta.speedErrors(1).identifier, ...
    "rampSpeed:nonconverged");
verifyThat(testCase,run.runMeta.speedErrors(1).message, ...
    matlab.unittest.constraints.ContainsSubstring("exitflag"));
end

function testProgressCountsOnlySolvedSpeeds(testCase)
[cars,~] = carConfig();
events = struct("phase",{}, "speedIndex",{}, ...
    "speed_mps",{}, "completedSpeeds",{}, "requestedSpeeds",{});
run = rampSpeed.runLongitudinalRamp(cars{1,2}, ...
    struct("speeds",[NaN 5],"verbose",false), ...
    struct("id","progress","label","progress","carRole","acceleration"), ...
    struct("onProgress",@captureProgress));

verifyEqual(testCase,run.status,"partial");
verifyEqual(testCase,[events.speedIndex],[1 2 2]);
verifyEqual(testCase,[events.completedSpeeds],[0 0 1]);

    function captureProgress(event)
        events(end+1) = event;
    end
end
