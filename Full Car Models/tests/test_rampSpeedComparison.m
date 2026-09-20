function tests = test_rampSpeedComparison
tests = functiontests(localfunctions);
end

function testBuildComparisonUsesVariantMinusBaselineOnCommonGrid(testCase)
baseline = makeComparisonRun("baseline",[5;10;15],[10;20;30],[true;true;true],[]);
variant = makeComparisonRun("variant",[5;12.5;15],[12;27;34],[true;true;true],[]);
runs = [baseline variant];
delta = rampSpeed.buildComparison(runs,"aero_front_load","baseline");
verifyEqual(testCase,delta.baselineId,"baseline");
verifyEqual(testCase,delta.interpolationMethod,"linear");
verifyEqual(testCase,delta.comparisonGrid_mps,[5;10;12.5;15]);
series = delta.series(strcmp(string({delta.series.id}),"variant"));
verifyEqual(testCase,series.values,[2;2;2;4],"AbsTol",1e-12);
verifyEqual(testCase,series.valid,true(4,1));
end

function testBuildComparisonNeverBridgesInvalidOrTruncatedGaps(testCase)
baseline = makeComparisonRun("baseline",[5;10;15],[10;20;30], ...
    [true;false;true],[false;false;false]);
variant = makeComparisonRun("variant",[5;10;15],[12;22;34], ...
    [true;true;true],[false;false;false]);
delta = rampSpeed.buildComparison([baseline variant], ...
    "aero_front_load","baseline",struct("grid",[5;10;15]));
series = delta.series(strcmp(string({delta.series.id}),"variant"));
verifyEqual(testCase,series.valid,[true;false;true]);
verifyTrue(testCase,isnan(series.values(2)));

baseline.perSpeed.valid = [true;true;true];
baseline.perSpeed.truncated = [false;true;false];
delta = rampSpeed.buildComparison([baseline variant], ...
    "aero_front_load","baseline",struct("grid",[5;10;15]));
series = delta.series(strcmp(string({delta.series.id}),"variant"));
verifyEqual(testCase,series.valid,[true;false;true]);
verifyTrue(testCase,isnan(series.values(2)));
end

function testBuildComparisonDoesNotExtrapolate(testCase)
baseline = makeComparisonRun("baseline",[5;10;15],[10;20;30], ...
    [true;true;true],[]);
variant = makeComparisonRun("variant",[6;10;14],[12;25;28], ...
    [true;true;true],[]);
delta = rampSpeed.buildComparison([baseline variant], ...
    "aero_front_load","baseline",struct("comparisonGrid_mps",[4;6;10;14;16]));
series = delta.series(strcmp(string({delta.series.id}),"variant"));
verifyEqual(testCase,series.valid,[false;true;true;true;false]);
verifyTrue(testCase,isnan(series.values(1)));
verifyTrue(testCase,isnan(series.values(5)));
end

function testComparisonCarriesBaselineAndGridMetadata(testCase)
baseline = makeComparisonRun("baseline",[5;10],[10;20],[true;true],[]);
variant = makeComparisonRun("variant",[5;10],[12;21],[true;true],[]);
delta = rampSpeed.buildComparison([baseline variant], ...
    "aero_front_load","baseline");
verifyEqual(testCase,delta.interpolationPolicy, ...
    "linear within contiguous valid non-truncated intervals");
verifyEqual(testCase,delta.variantMinusBaseline,true);
verifyEqual(testCase,delta.series(1).baselineId,"baseline");
end

function testComparisonRendersWithResampledGridFlags(testCase)
baseline = makeComparisonRun("baseline",[5;10;15],[10;20;30], ...
    [true;true;true],[false;false;false]);
variant = makeComparisonRun("variant",[5;12.5;15],[12;27;34], ...
    [true;true;true],[false;false;false]);
delta = rampSpeed.buildComparison([baseline variant],"aero_front_load", ...
    "baseline",struct("comparisonGrid_mps",[5;10;12.5;15]));
series = delta.series(1);
verifyEqual(testCase,numel(series.x),4);
verifySize(testCase,series.valid,[4 1]);
verifySize(testCase,series.truncated,[4 1]);
verifySize(testCase,series.power_limited,[4 1]);
verifySize(testCase,series.wheel_lift,[4 1]);
fig = figure("Visible","off");
cleanup = onCleanup(@()close(fig));
ax = axes(fig);
h = rampSpeed.renderMetric(ax,delta,struct("showLegend",false));
verifyEqual(testCase,numel(h.lines),1);
verifyTrue(testCase,all(isgraphics(h.lines)));
end

function run = makeComparisonRun(id,speeds,values,valid,truncated)
speeds = double(speeds(:));
values = double(values(:));
n = numel(speeds);
if isempty(truncated)
    truncated = false(n,1);
else
    truncated = logical(truncated(:));
end
valid = logical(valid(:));
caseInfo = struct("id",string(id),"label",string(id),"source","test", ...
    "designRow",1,"sourceIndex",1,"carRole","lap","carColumn",1);
run = rampSpeed.makeRun("lateral","coast", ...
    struct("speeds",speeds.'),caseInfo);
run.perSpeed.front_downforce_N = values;
run.perSpeed.rear_downforce_N = 0.75*values;
run.perSpeed.downforce_N = 1.75*values;
run.perSpeed.drag_N = 0.1*values;
run.perSpeed.aero_balance_front = 0.5*ones(n,1);
run.perSpeed.mechanical_balance_front = 0.5*ones(n,1);
run.perSpeed.K_linear_rad_per_mps2 = ones(n,1);
run.perSpeed.valid = valid;
run.perSpeed.status(:) = "complete";
run.perSpeed.reason(:) = "";
run.perSpeed.reason(~valid) = "invalid";
run.perSpeed.truncated = truncated;
run.perSpeed.power_limited = false(n,1);
run.perSpeed.wheel_lift = false(n,1);
run.perSpeed.aero_outside_map = false(n,1);
run.status = "complete";
run.runMeta.status = "complete";
end
