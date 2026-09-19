function tests = test_rampSpeedSchema
tests = functiontests(localfunctions);
end
function testCreatesVersionedLateralRun(testCase)
fixture = makeRampFixture(); %#ok<NASGU>
run = rampSpeed.makeRun("lateral","coast",struct("speeds",[5 10]), ...
    struct("id","baseline","label","baseline","carRole","lap"));

verifyEqual(testCase,run.schemaVersion,1);
verifyEqual(testCase,run.caseId,"baseline");
verifyEqual(testCase,run.type,"lateral");
verifyEqual(testCase,run.mode,"coast");
verifyEqual(testCase,run.status,"pending");
verifyTrue(testCase,all(ismember(["perSpeed","points","runMeta", ...
    "status","raw"],string(fieldnames(run)))));
verifyEqual(testCase,height(run.perSpeed),2);
verifyTrue(testCase,all(ismember(["speed_mps","valid","status", ...
    "aLat_free_mps2","aLat_sustainable_mps2", ...
    "mechanical_balance_front","aero_balance_front"], ...
    string(run.perSpeed.Properties.VariableNames))));
verifyTrue(testCase,all(isnan(run.perSpeed.front_ride_height_m)));
verifyTrue(testCase,run.runMeta.lateralMetricsApplicable);
end

function testLongitudinalOnlyFieldsAreUnavailable(testCase)
run = rampSpeed.makeRun("longitudinal","",struct("speeds",5), ...
    struct("id","accel","label","accel","carRole","acceleration"));

verifyFalse(testCase,run.runMeta.lateralMetricsApplicable);
verifyTrue(testCase,all(isnan(run.perSpeed.K_linear_rad_per_mps2)));
verifyTrue(testCase,all(isnan(run.perSpeed.mechanical_balance_front)));
verifyTrue(testCase,all(isnan(run.perSpeed.aero_balance_front)) == false || ...
    all(isnan(run.perSpeed.aero_balance_front)));
verifyFalse(testCase,run.perSpeed.lateral_metrics_applicable(1));
end

function testMakeStudyCreatesVersionedContainer(testCase)
study = rampSpeed.makeStudy("test-app");

verifyEqual(testCase,study.schemaVersion,1);
verifyEqual(testCase,study.appVersion,"test-app");
verifyTrue(testCase,isdatetime(study.created));
verifyTrue(testCase,isempty(study.cases));
verifyTrue(testCase,isempty(study.runs));
verifyTrue(testCase,isfield(study,"displayUnits"));
end

function testDisplayConversionDoesNotMutateInput(testCase)
values = [1 2];
displayed = rampSpeed.displayUnits(values,"m","ft");

verifyEqual(testCase,displayed,values/0.3048,"AbsTol",1e-12);
verifyEqual(testCase,values,[1 2]);
end
