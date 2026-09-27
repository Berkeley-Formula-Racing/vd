function tests = test_rampSpeedAdaptiveSpeed
tests = functiontests(localfunctions);
end

function testAccuracyModesSetBaseAndLocalSpacing(testCase)
accurate = rampSpeed.adaptiveSpeedPolicy("accurate");
preview = rampSpeed.adaptiveSpeedPolicy("preview");
highAccuracy = rampSpeed.adaptiveSpeedPolicy("highAccuracy");

verifyEqual(testCase,accurate.baseSpacing_mps,2.5);
verifyEqual(testCase,preview.baseSpacing_mps,5);
verifyEqual(testCase,accurate.balanceTolerance_fraction,0.0025);
verifyEqual(testCase,accurate.relativeTolerance,0.005);
verifyEqual(testCase,accurate.minRefinementSpacing_mps,1.25);
verifyEqual(testCase,highAccuracy.minRefinementSpacing_mps,0.625);
verifyEqual(testCase,accurate.maxPasses,3);
verifyEqual(testCase,accurate.maxPoints,33);
end

function testLongitudinalEntryPointRunsAndRecordsAdaptiveSeedGrid(testCase)
[car,~] = carConfigBaseline();
settings = struct("speeds",[5 10],"speedGrid", ...
    struct("mode","accurate","maxPasses",1,"maxPoints",6),"verbose",false);
caseInfo = struct("id","baseline","label","baseline","carRole","auto");

run = rampSpeed.runLongitudinalRamp(car,settings,caseInfo,struct());

verifyEqual(testCase,run.perSpeed.speed_mps, ...
    sort(unique(run.perSpeed.speed_mps)));
verifyGreaterThan(testCase,numel(run.perSpeed.speed_mps),3);
verifyEqual(testCase,run.settings.speeds,run.perSpeed.speed_mps.');
verifyEqual(testCase,run.runMeta.requestedSpeeds_mps,[5; 10]);
verifyEqual(testCase,run.runMeta.speedGrid.passes,1);
verifyEqual(testCase,run.runMeta.speedGrid.stopReason,"max_passes");
verifyEqual(testCase,run.runMeta.speedGrid.refinementHistory(end).stopReason, ...
    "max_passes");
verifyTrue(testCase,all(run.perSpeed.valid));
end

function testAdaptiveRequestedOffGridSpeedRunsOnce(testCase)
[car,~] = carConfigBaseline();
settings = struct("solverProfile","accurate","speeds",26.25, ...
    "speedGrid",struct("mode","accurate","maxPasses",0), ...
    "verbose",false);
run = rampSpeed.runLongitudinalRamp(car,settings, ...
    struct("id","baseline","label","baseline","carRole","auto"),struct());

verifyEqual(testCase,height(run.perSpeed),1);
verifyEqual(testCase,numel(unique(run.perSpeed.speed_mps)),1);
verifyEmpty(testCase,run.runMeta.speedGrid.exactRetrySpeeds_mps);
verifyTrue(testCase,run.perSpeed.valid);
verifyEqual(testCase,run.perSpeed.status,"complete");
verifyEqual(testCase,run.perSpeed.speedIndex,1);
end
function testPlanRetainsRequestedSpeedsOnModeSpecificBaseGrid(testCase)
accurate = rampSpeed.adaptiveSpeedPolicy("accurate");
preview = rampSpeed.adaptiveSpeedPolicy("preview");

accuratePlan = rampSpeed.planAdaptiveSpeeds([5 10 18 25],accurate);
previewPlan = rampSpeed.planAdaptiveSpeeds([5 10 18 25],preview);

verifyEqual(testCase,accuratePlan.seedSpeeds_mps,[5; 7.5; 10; 12.5; ...
    15; 17.5; 18; 20; 22.5; 25]);
verifyEqual(testCase,previewPlan.seedSpeeds_mps,[5; 10; 15; 18; 20; 25]);
requested = accuratePlan.provenance.source == "requested";
verifyTrue(testCase,all(ismember([5; 10; 18; 25], ...
    accuratePlan.provenance.speed_mps(requested))));
end

function testPlanRejectsPointLimitThatWouldBreakBaseGrid(testCase)
policy = rampSpeed.adaptiveSpeedPolicy("accurate",struct("maxPoints",4));

verifyError(testCase,@() rampSpeed.planAdaptiveSpeeds(0:10,policy), ...
    "rampSpeed:adaptivePointLimit");
end

function testBalanceCurvatureInsertsLocalMidpointsWithProvenance(testCase)
policy = rampSpeed.adaptiveSpeedPolicy("accurate");
plan = rampSpeed.planAdaptiveSpeeds(0:2.5:10,policy);
run = longitudinalRun([0; 2.5; 5; 7.5; 10], ...
    [0.50; 0.50; 0.53; 0.54; 0.55], ...
    [100; 102.5; 105; 107.5; 110], ...
    [2; 2; 2; 2; 2],true(5,1),[1; 1; 1; 1; 1]);

[refined,report] = rampSpeed.refineAdaptiveSpeeds(plan,run);

verifyTrue(testCase,all(ismember([1.25; 3.75; 6.25], ...
    report.insertedSpeeds_mps)));
verifyTrue(testCase,any(contains(refined.provenance.reason, ...
    "balance_error")));
verifyEqual(testCase,refined.provenance.source( ...
    ismember(refined.provenance.speed_mps,report.insertedSpeeds_mps)), ...
    repmat("adaptive",numel(report.insertedSpeeds_mps),1));
end

function testForceAndAccelerationInterpolationErrorsRefine(testCase)
policy = rampSpeed.adaptiveSpeedPolicy("accurate");
plan = rampSpeed.planAdaptiveSpeeds(0:2.5:10,policy);
run = longitudinalRun([0; 2.5; 5; 7.5; 10], ...
    repmat(0.50,5,1),[100; 101; 110; 110; 120], ...
    [2; 2; 2.4; 2; 2],true(5,1),ones(5,1));

[~,report] = rampSpeed.refineAdaptiveSpeeds(plan,run);

verifyTrue(testCase,all(ismember([3.75; 6.25], ...
    report.insertedSpeeds_mps)));
verifyTrue(testCase,any(contains(report.insertReasons,"force_error")));
verifyTrue(testCase,any(contains(report.insertReasons,"acceleration_error")));
end

function testGearTransitionAndRetryWarningsRefineNeighbors(testCase)
policy = rampSpeed.adaptiveSpeedPolicy("accurate");
plan = rampSpeed.planAdaptiveSpeeds(0:2.5:10,policy);
run = longitudinalRun([0; 2.5; 5; 7.5; 10], ...
    repmat(0.50,5,1),[100; 102.5; 105; 107.5; 110], ...
    repmat(2,5,1),true(5,1),[1; 1; 2; 2; 2]);
run.runMeta.retrySpeeds_mps = 7.5;

[~,report] = rampSpeed.refineAdaptiveSpeeds(plan,run);

verifyTrue(testCase,ismember(3.75,report.insertedSpeeds_mps));
verifyTrue(testCase,ismember(6.25,report.insertedSpeeds_mps));
verifyTrue(testCase,any(contains(report.insertReasons,"state_transition")));
verifyTrue(testCase,any(contains(report.insertReasons,"solver_retry")));
end

function testInvalidRowsRemainGapsAndDoNotContributeToInterpolation(testCase)
policy = rampSpeed.adaptiveSpeedPolicy("accurate");
plan = rampSpeed.planAdaptiveSpeeds(0:2.5:10,policy);
valid = [true; true; false; true; true];
balance = [0.50; 0.50; NaN; 0.90; 0.90];
run = longitudinalRun([0; 2.5; 5; 7.5; 10],balance, ...
    [100; 102.5; NaN; 107.5; 110], ...
    [2; 2; NaN; 2; 2],valid,ones(5,1));

[~,report] = rampSpeed.refineAdaptiveSpeeds(plan,run);

verifyEqual(testCase,report.invalidSpeeds_mps,5);
verifyTrue(testCase,all(ismember([3.75; 6.25], ...
    report.insertedSpeeds_mps)));
verifyFalse(testCase,any(contains(report.insertReasons,"balance_error")));
verifyTrue(testCase,all(contains(report.insertReasons,"invalid_neighbor")));
end

function testPassAndPointLimitsStopRefinementExplicitly(testCase)
noPasses = rampSpeed.adaptiveSpeedPolicy("accurate", ...
    struct("maxPasses",0));
plan = rampSpeed.planAdaptiveSpeeds(0:2.5:10,noPasses);
run = longitudinalRun([0; 2.5; 5; 7.5; 10], ...
    [0.5; 0.5; 0.55; 0.5; 0.5], ...
    [100; 102.5; 105; 107.5; 110], ...
    repmat(2,5,1),true(5,1),ones(5,1));

[stopped,report] = rampSpeed.refineAdaptiveSpeeds(plan,run);

verifyEmpty(testCase,report.insertedSpeeds_mps);
verifyEqual(testCase,report.stopReason,"max_passes");
verifyEqual(testCase,stopped.status,"limited");
end

function testCacheKeysAreStableAndIncludeAllRunIdentity(testCase)
contextA = cacheContext();
contextB = struct("aeroMapId",contextA.aeroMapId, ...
    "baselineConfigVersion",contextA.baselineConfigVersion, ...
    "solverProfileId",contextA.solverProfileId, ...
    "setupProvenance",contextA.setupProvenance, ...
    "setupId",contextA.setupId,"solverOptions",contextA.solverOptions);

keyA = rampSpeed.adaptiveResultCacheKey(contextA,20);
keyB = rampSpeed.adaptiveResultCacheKey(contextB,20);
verifyEqual(testCase,keyA,keyB);
verifyNotEqual(testCase,keyA, ...
    rampSpeed.adaptiveResultCacheKey(contextA,22.5));
contextB.solverProfileId = "fastPreview";
verifyNotEqual(testCase,keyA, ...
    rampSpeed.adaptiveResultCacheKey(contextB,20));
contextB = contextA;
contextB.setupProvenance.setupRevision = 2;
verifyNotEqual(testCase,keyA, ...
    rampSpeed.adaptiveResultCacheKey(contextB,20));
end

function testCacheReturnsValidHitsAndRetainsInvalidDiagnostics(testCase)
context = cacheContext();
cache = rampSpeed.createAdaptiveResultCache();
run = longitudinalRun([5; 7.5],[0.5; 0.5],[100; NaN], ...
    [2; NaN],[true; false],[1; 1]);
run.points = table([5; 7.5],[true; false], ...
    'VariableNames',{'speed_mps','valid'});
run.raw = struct("diagnostics",[ ...
    struct("speed_mps",5,"state",1:9), ...
    struct("speed_mps",7.5,"state",[])]);
provenance = table([5; 7.5],["requested"; "seed"],[0; 0], ...
    ["user"; "base_grid"], ...
    'VariableNames',{'speed_mps','source','pass','reason'});

[cache,inserted] = rampSpeed.cacheAdaptiveRampRun( ...
    cache,context,run,provenance);
[payloads,validHit,found] = rampSpeed.lookupAdaptiveRampResults( ...
    cache,context,[5; 7.5; 10]);

verifyEqual(testCase,inserted,[true; true]);
verifyEqual(testCase,validHit,[true; false; false]);
verifyEqual(testCase,found,[true; true; false]);
verifyEqual(testCase,payloads{1}.perSpeedRow.aLong_mps2,2);
verifyEqual(testCase,payloads{1}.diagnostic.state,1:9);
verifyEqual(testCase,payloads{1}.provenance.source,"requested");
verifyEqual(testCase,payloads{2}.perSpeedRow.status,"invalid");
end

function run = longitudinalRun(speeds,aeroBalance,downforce,acceleration, ...
    valid,currentGear)
n = numel(speeds);
drag = 0.4*downforce;
frontAxle = 1500 + 0.5*downforce(:);
rearAxle = 1500 + 0.5*downforce(:);
run = struct();
run.type = "longitudinal";
run.perSpeed = table(speeds(:),valid(:),strings(n,1),aeroBalance(:), ...
    downforce(:),drag(:),acceleration(:),frontAxle,rearAxle, ...
    0.3*downforce(:),0.7*downforce(:),0.1*downforce(:), ...
    false(n,1),false(n,1),false(n,1),false(n,1), ...
    currentGear(:),zeros(n,1), ...
    'VariableNames',{'speed_mps','valid','status','aero_balance_front', ...
    'downforce_N','drag_N','aLong_mps2','Fz_front_axle_N', ...
    'Fz_rear_axle_N','LLT_front_N','LLT_rear_N','long_load_transfer_N', ...
    'power_limited','traction_limited','wheel_lift', ...
    'aero_outside_map','current_gear','max_constraint_residual'});
run.perSpeed.status(valid(:)) = "complete";
run.perSpeed.status(~valid(:)) = "invalid";
run.runMeta = struct("warnings",strings(0,1));
run.settings = struct("solverOptions",struct( ...
    "constraintTolerance",1e-2));
run.points = table();
run.raw = struct();
end

function context = cacheContext()
context = struct("setupId","setup-01", ...
    "setupProvenance",struct("setupRevision",1,"source","user"), ...
    "solverProfileId","accurate", ...
    "baselineConfigVersion","ramp-baseline-v1", ...
    "aeroMapId","aeromap_b26", ...
    "solverOptions",struct("constraintTolerance",1e-2, ...
    "maxFunctionEvaluations",2000));
end
