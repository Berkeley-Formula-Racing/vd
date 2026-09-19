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
verifyEqual(testCase,run.perSpeed.status(1),"complete");
verifyFalse(testCase,any(run.perSpeed.valid(2:end)));
verifyTrue(testCase,all(run.perSpeed.status(2:end) == "failed"));
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
