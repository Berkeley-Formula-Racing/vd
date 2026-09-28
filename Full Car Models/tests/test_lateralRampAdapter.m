function tests = test_lateralRampAdapter
tests = functiontests(localfunctions);
end

function testAdapterMatchesDirectRamp(testCase)
[cars,~] = carConfig();
settings = struct("speeds",[5 10],"nRamp",4,"nBisect",0, ...
    "mode","coast","verbose",false);
direct = rampSweep(cars{1,1},settings);
run = rampSpeed.runLateralRamp(cars{1,1},settings, ...
    struct("id","baseline","label","baseline","carRole","lap"),struct());

verifyEqual(testCase,run.raw.perSpeed.vCar,direct.perSpeed.vCar);
verifyEqual(testCase,run.raw.perSpeed.K_linear,direct.perSpeed.K_linear, ...
    "AbsTol",1e-12);
verifyEqual(testCase,run.perSpeed.speed_mps, ...
    direct.perSpeed.vCar,"AbsTol",1e-12);
verifyEqual(testCase,run.perSpeed.mechanical_balance_front, ...
    direct.perSpeed.mech_balance,"AbsTol",1e-12);
verifyEqual(testCase,run.perSpeed.aLat_free_mps2, ...
    direct.perSpeed.gLat_max*9.80665,"AbsTol",1e-12);
verifyEqual(testCase,run.perSpeed.aLat_sustainable_mps2, ...
    direct.perSpeed.gLat_top*9.80665,"AbsTol",1e-12);
end

function testCancellationRetainsCompletedRows(testCase)
[cars,~] = carConfig();
settings = struct("speeds",[5 10 15],"nRamp",4,"nBisect",0, ...
    "mode","coast","verbose",false);
cancelRequested = false;
callbacks = struct("onProgress",@captureProgress, ...
    "isCancelled",@isCancelled);

run = rampSpeed.runLateralRamp(cars{1,1},settings, ...
    struct("id","baseline","label","baseline","carRole","lap"),callbacks);

verifyEqual(testCase,run.status,"cancelled");
verifyEqual(testCase,run.perSpeed.speed_mps,settings.speeds(:), ...
    "AbsTol",1e-12);
verifyTrue(testCase,run.perSpeed.valid(1));
verifyEqual(testCase,run.perSpeed.status(1),"converged");
verifyFalse(testCase,any(run.perSpeed.valid(2:end)));
verifyTrue(testCase,all(run.perSpeed.status(2:end) == "cancelled"));
verifyThat(testCase,run.perSpeed.reason(2), ...
    matlab.unittest.constraints.ContainsSubstring("cancel"));

    function captureProgress(event)
        if isfield(event,"completedSpeeds") && event.completedSpeeds >= 1
            cancelRequested = true;
        end
    end

    function value = isCancelled()
        value = cancelRequested;
    end
end

function testOmittedModeMatchesRampSweepDefault(testCase)
[cars,~] = carConfig();
settings = struct("speeds",5,"nRamp",3,"nBisect",0,"verbose",false);
direct = rampSweep(cars{1,1},settings);
run = rampSpeed.runLateralRamp(cars{1,1},settings, ...
    struct("id","baseline","label","baseline","carRole","lap"),struct());

verifyEqual(testCase,string(direct.settings.mode),"balanced");
verifyEqual(testCase,run.mode,"balanced");
verifyEqual(testCase,string(run.raw.settings.mode),run.mode);
verifyEqual(testCase,run.settings.mode,run.mode);
end

function testSkippedSpeedRetainsPerSpeedErrorDiagnostics(testCase)
[cars,~] = carConfig();
settings = struct("speeds",[0 5],"nRamp",3,"nBisect",0, ...
    "mode","coast","verbose",false);
run = rampSpeed.runLateralRamp(cars{1,1},settings, ...
    struct("id","baseline","label","baseline","carRole","lap"),struct());

verifyFalse(testCase,run.perSpeed.valid(1));
verifyEqual(testCase,run.perSpeed.status(1),"solver_failed");
verifyTrue(testCase,isfield(run.raw,"speedErrors"));
errors = run.raw.speedErrors;
verifyGreaterThanOrEqual(testCase,numel(errors),1);
errorIndex = find([errors.speed_index] == 1,1);
verifyNotEmpty(testCase,errorIndex);
speedError = errors(errorIndex);
verifyTrue(testCase,isfield(speedError,"stack"));
verifyThat(testCase,run.perSpeed.reason(1), ...
    matlab.unittest.constraints.ContainsSubstring(string(speedError.identifier)));
verifyThat(testCase,run.perSpeed.reason(1), ...
    matlab.unittest.constraints.ContainsSubstring(string(speedError.message)));
verifyEqual(testCase,run.runMeta.speedErrors,errors);
end

function testOmittedSettingsSpeedsRetainDefaultGridDiagnostics(testCase)
[cars,~] = carConfig();
settings = struct("nRamp",1,"nBisect",0,"mode","balanced", ...
    "ayMinFrac",0.99,"ayMaxFrac",1,"verbose",false);
run = rampSpeed.runLateralRamp(cars{1,1},settings, ...
    struct("id","baseline","label","baseline","carRole","lap"),struct());

requested = (5:2.5:30).';
verifyEqual(testCase,string(run.raw.status),"partial");
verifyEqual(testCase,run.raw.settings.speeds,requested.', ...
    "AbsTol",1e-12);
verifyEqual(testCase,run.perSpeed.speed_mps,requested,"AbsTol",1e-12);
verifyEqual(testCase,height(run.perSpeed),numel(requested));
verifyTrue(testCase,any(run.perSpeed.valid));
verifyTrue(testCase,isfield(run.raw,"speedErrors"));
errors = run.raw.speedErrors;
verifyGreaterThanOrEqual(testCase,numel(errors),1);
errorIndex = find([errors.speed_index] == 1,1);
verifyNotEmpty(testCase,errorIndex);
speedError = errors(errorIndex);
target = speedError.speed_index;
verifyFalse(testCase,run.perSpeed.valid(target));
verifyEqual(testCase,run.perSpeed.status(target),"infeasible");
verifyTrue(testCase,isfield(speedError,"stack"));
verifyThat(testCase,run.perSpeed.reason(target), ...
    matlab.unittest.constraints.ContainsSubstring(string(speedError.identifier)));
verifyThat(testCase,run.perSpeed.reason(target), ...
    matlab.unittest.constraints.ContainsSubstring(string(speedError.message)));
verifyEqual(testCase,run.runMeta.speedErrors,errors);
end
