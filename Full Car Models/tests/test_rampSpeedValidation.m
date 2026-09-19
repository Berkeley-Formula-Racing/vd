function tests = test_rampSpeedValidation
tests = functiontests(localfunctions);
end

function testValidStudyPassesValidation(testCase)
study = makeValidStudy();
[ok,issues] = rampSpeed.validateStudy(study);

verifyTrue(testCase,ok);
verifyEmpty(testCase,issues);
end

function testValidationReportsDuplicateCaseIds(testCase)
study = makeValidStudy();
study.cases(2) = study.cases(1);

[ok,issues] = rampSpeed.validateStudy(study);

verifyFalse(testCase,ok);
verifyTrue(testCase,contains(string(issues),"duplicate"));
end

function testValidationRejectsUnsupportedRunType(testCase)
study = makeValidStudy();
study.runs(1).type = "combined";

[ok,issues] = rampSpeed.validateStudy(study);

verifyFalse(testCase,ok);
verifyTrue(testCase,contains(string(issues),"unsupported"));
end

function testValidationRejectsNonSiCanonicalColumn(testCase)
study = makeValidStudy();
study.runs(1).perSpeed.front_ride_height_in = ...
    ones(height(study.runs(1).perSpeed),1);

[ok,issues] = rampSpeed.validateStudy(study);

verifyFalse(testCase,ok);
verifyTrue(testCase,contains(string(issues),"non-SI"));
end

function testNormalizationConvertsLegacyUnitsToSi(testCase)
fixture = makeRampFixture();
raw = fixture.legacyRampResult;
run = rampSpeed.normalizeRampResult(raw,"lateral",raw.settings, ...
    fixture.cases(1),struct("source","fixture"));

verifyEqual(testCase,run.perSpeed.front_ride_height_m(1), ...
    fixture.frontRideHeightIn*0.0254,"AbsTol",1e-12);
verifyEqual(testCase,run.perSpeed.aLat_sustainable_mps2(1), ...
    raw.perSpeed.gLat_top(1)*9.80665,"AbsTol",1e-12);
verifyEqual(testCase,run.perSpeed.K_linear_rad_per_mps2(1), ...
    raw.perSpeed.K_linear(1)*pi/180/9.80665,"AbsTol",1e-12);
verifyTrue(testCase,all(run.perSpeed.valid));
end

function testNormalizationLeavesLongitudinalLateralMetricsUnavailable(testCase)
fixture = makeRampFixture();
run = rampSpeed.normalizeRampResult(struct("perSpeed",table(5, ...
    'VariableNames',{'vCar'})),"longitudinal", ...
    struct("speeds",5),fixture.cases(2),struct());

verifyFalse(testCase,run.runMeta.lateralMetricsApplicable);
verifyTrue(testCase,isnan(run.perSpeed.K_linear_rad_per_mps2(1)));
verifyTrue(testCase,isnan(run.perSpeed.mechanical_balance_front(1)));
end

function testNormalizationRetainsMissingRequestedSpeedAsInvalidRow(testCase)
fixture = makeRampFixture();
raw = struct("perSpeed",table(5,'VariableNames',{'vCar'}));
run = rampSpeed.normalizeRampResult(raw,"lateral", ...
    struct("speeds",[5 10]),fixture.cases(1),struct());

verifyEqual(testCase,run.perSpeed.speed_mps,[5;10],"AbsTol",1e-12);
verifyTrue(testCase,run.perSpeed.valid(1));
verifyFalse(testCase,run.perSpeed.valid(2));
verifyEqual(testCase,run.perSpeed.status(2),"missing");
end

function study = makeValidStudy()
run = makeFixtureRun(struct("id","baseline","label","baseline", ...
    "carRole","lap"));
study = rampSpeed.makeStudy("test-app");
study.cases = struct("id","baseline","label","baseline", ...
    "source","fixture","designRow",1,"carRole","lap");
study.runs = run;
end
