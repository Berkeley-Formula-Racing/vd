function tests = test_rampSpeedLateralExecutionContract
tests = functiontests(localfunctions);
end

function testSuccessfulAndFailedInteriorRowsRemainExplicit(testCase)
[cars,~] = carConfig();
settings = struct("speeds",[5 0 10],"nRamp",3,"nBisect",0, ...
    "mode","coast","verbose",false,"solverProfile","fastPreview");
run = rampSpeed.runLateralRamp(cars{1,1},settings, ...
    struct("id","lateral-contract","label","lateral-contract", ...
    "carRole","lap"),struct());

verifyEqual(testCase,run.perSpeed.speed_mps,settings.speeds(:), ...
    "AbsTol",1e-12);
verifyEqual(testCase,height(run.perSpeed),3);
verifyTrue(testCase,run.perSpeed.valid(1));
verifyFalse(testCase,run.perSpeed.valid(2));
verifyTrue(testCase,run.perSpeed.valid(3));
verifyEqual(testCase,run.perSpeed.status(2),"failed");
verifyThat(testCase,run.perSpeed.reason(2), ...
    matlab.unittest.constraints.ContainsSubstring("speed"));
verifyTrue(testCase,isnan(run.perSpeed.aLat_free_mps2(2)));
verifyTrue(testCase,isnan(run.perSpeed.mechanical_balance_front(2)));
end

function testSuccessfulPointRetainsRawRampPayload(testCase)
[cars,~] = carConfig();
task = struct("speedIndex",7,"speed_mps",5,"origin","requested", ...
    "passIndex",1);
settings = struct("nRamp",3,"nBisect",0,"mode","coast", ...
    "verbose",false,"solverProfile","fastPreview");
profile = rampSpeed.resolveSolverProfile("fastPreview");
result = rampSpeed.solveLateralPoint(cars{1,1},task,settings,profile,struct());

verifyEqual(testCase,result.speedIndex,7);
verifyEqual(testCase,result.speed_mps,5,"AbsTol",1e-12);
verifyTrue(testCase,ismember(result.status,["converged","near_feasible"]));
verifyTrue(testCase,isfield(result.diagnostics,"raw"));
verifyTrue(testCase,istable(result.diagnostics.raw.points));
verifyGreaterThan(testCase,height(result.diagnostics.raw.points),0);
verifyTrue(testCase,isfield(result.metrics,"gLat_max"));
verifyTrue(testCase,isfield(result.metrics,"CUOsteerFromYaw_limit_deg"));
end

function testCancellationReturnsCancelledPoint(testCase)
[cars,~] = carConfig();
task = struct("speedIndex",3,"speed_mps",5,"origin","requested", ...
    "passIndex",1);
settings = struct("nRamp",3,"nBisect",0,"mode","coast", ...
    "verbose",false,"solverProfile","fastPreview");
control = struct("shouldCancel",@() true);
profile = rampSpeed.resolveSolverProfile("fastPreview");
result = rampSpeed.solveLateralPoint(cars{1,1},task,settings,profile,control);

verifyEqual(testCase,result.status,"cancelled");
verifyThat(testCase,result.diagnostics.reason, ...
    matlab.unittest.constraints.ContainsSubstring("cancel"));
end
