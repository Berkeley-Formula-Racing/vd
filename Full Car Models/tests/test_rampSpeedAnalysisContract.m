function tests = test_rampSpeedAnalysisContract
%TEST_RAMPSPEEDANALYSISCONTRACT Task 7 analysis and migration contracts.

tests = functiontests(localfunctions);
end

function testLongitudinalCatalogUsesCapabilityAndSeparateDragDeceleration(testCase)
catalog = rampSpeed.metricCatalog(struct("type","longitudinal"));
ids = string({catalog.id});
capability = catalog(ids == "capability_longitudinal");
dragDecel = catalog(ids == "drag_deceleration");

verifyEqual(testCase,numel(capability),1);
verifyEqual(testCase,string(capability.sourceLevel),"perSpeed");
verifyEqual(testCase,string(capability.source),"aLong_max_mps2");
verifyEqual(testCase,string(capability.xSource),"speed_mps");
verifyTrue(testCase,any(string(capability.supportedRampTypes) == ...
    "longitudinal"));
verifyEqual(testCase,numel(dragDecel),1);
verifyEqual(testCase,string(dragDecel.source),"drag_N");
verifyEqual(testCase,string(dragDecel.units),"m/s^2");
verifyEqual(testCase,string(dragDecel.truncationPolicy),"allow");
end

function testCornerCamberAndSlipAngleMetricsAreDeclared(testCase)
catalog = rampSpeed.metricCatalog(struct("type","lateral"));
ids = string({catalog.id});
expected = ["camber_FL","camber_FR","camber_RL","camber_RR", ...
    "slip_angle_FL","slip_angle_FR","slip_angle_RL","slip_angle_RR"];
for id = expected
    metric = catalog(ids == id);
    verifyEqual(testCase,numel(metric),1);
    verifyEqual(testCase,string(metric.sourceLevel),"points");
    verifyEqual(testCase,string(metric.xSource),"speed_mps");
    verifyTrue(testCase,any(string(metric.supportedRampTypes) == "lateral"));
end
verifyEqual(testCase,string(catalog(ids == "camber_FL").source), ...
    "gamma_FL_rad");
verifyEqual(testCase,string(catalog(ids == "slip_angle_RR").source), ...
    "alpha_RR_rad");
end

function testLegacyPointsRegainStableSpeedAndPointIndices(testCase)
speeds = (5:15).';
pointsPerSpeed = 12;
pointSpeeds = repelem(speeds,pointsPerSpeed);
raw = struct();
raw.settings = struct("speeds",speeds.');
raw.perSpeed = table(speeds,true(11,1),ones(11,1),zeros(11,1), ...
    'VariableNames',{'vCar','valid','exitflag','max_ceq'});
raw.points = table(pointSpeeds,true(numel(pointSpeeds),1), ...
    ones(numel(pointSpeeds),1),zeros(numel(pointSpeeds),1), ...
    'VariableNames',{'vCar','valid','exitflag','max_ceq'});

run = rampSpeed.normalizeRampResult(raw,"lateral",raw.settings, ...
    struct("id","legacy"),struct("source","test"));
verifyEqual(testCase,numel(unique(run.points.speed_index)),11);
for speedIndex = 1:11
    rows = run.points.speed_index == speedIndex;
    verifyEqual(testCase,run.points.point_index(rows), ...
        (1:pointsPerSpeed).');
end
end

function testRawPlotUsesActualPointCoordinateAndPreservesInvalidGap(testCase)
fixture = makeRampFixture();
run = fixture.lateralRun;
run.points.speed_mps = [0.25;0.75];
run.points.speed_index = [1;1];
run.points.point_index = [1;2];
run.points.valid = [true;false];
run.points.aLat_mps2 = [1.2;NaN];
data = rampSpeed.buildPlotData(run,"raw_aLat");

series = data.series(1);
verifyEqual(testCase,series.x,[0.25;0.75],"AbsTol",eps);
verifyEqual(testCase,series.values(1),1.2/9.80665,"AbsTol",1e-12);
verifyTrue(testCase,isnan(series.values(2)));
verifyFalse(testCase,series.valid(2));
end

function testComparisonUsesMetricSpecificTruncationPolicy(testCase)
fixture = makeRampFixture();
baseline = fixture.lateralRun;
variant = fixture.lateralRun;
baseline.caseId = "baseline";
variant.caseId = "variant";
baseline.perSpeed.truncated(2) = true;
variant.perSpeed.truncated(2) = true;

capability = rampSpeed.buildComparison([baseline,variant], ...
    "capability_sustainable","baseline", ...
    struct("comparisonGrid_mps",[5;10]));
aero = rampSpeed.buildComparison([baseline,variant], ...
    "downforce","baseline", ...
    struct("comparisonGrid_mps",[5;10]));

verifyFalse(testCase,capability.series(1).valid(2));
verifyTrue(testCase,aero.series(1).valid(2));
end

function testStatusesAndDragDecelerationRemainTyped(testCase)
speeds = [5;10];
raw = struct();
raw.settings = struct("speeds",speeds.','mass_kg',1000);
raw.perSpeed = table(speeds,["converged";"cancelled"], ...
    [true;false],[1;1],[0;0],[20;30], ...
    'VariableNames',{'vCar','status','valid','exitflag', ...
    'max_ceq','drag'});
raw.points = table();
run = rampSpeed.normalizeRampResult(raw,"longitudinal",raw.settings, ...
    struct("id","status-test"),struct("source","test"));

verifyEqual(testCase,string(run.perSpeed.status), ...
    ["converged";"cancelled"]);
verifyFalse(testCase,run.perSpeed.valid(2));
plotData = rampSpeed.buildPlotData(run,"drag_deceleration");
verifyEqual(testCase,plotData.series(1).values, ...
    [0.02;NaN],"AbsTol",1e-12);
verifyTrue(testCase,all(isnan(run.perSpeed.aLat_free_mps2)));
end
