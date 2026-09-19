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

function testValidationReportsMissingRequiredStudyField(testCase)
study = makeValidStudy();
study = rmfield(study,"displayUnits");

[ok,issues] = rampSpeed.validateStudy(study);

verifyFalse(testCase,ok);
verifyTrue(testCase,contains(string(issues),"missing study field: displayUnits"));
end

function testValidationReportsMissingRequiredRunField(testCase)
study = makeValidStudy();
study.runs = rmfield(study.runs,"raw");

[ok,issues] = rampSpeed.validateStudy(study);

verifyFalse(testCase,ok);
verifyTrue(testCase,contains(string(issues),"missing run field: raw"));
end

function testValidationRejectsInvalidRunStatus(testCase)
study = makeValidStudy();
study.runs.status = "not-a-status";

[ok,issues] = rampSpeed.validateStudy(study);

verifyFalse(testCase,ok);
verifyTrue(testCase,contains(string(issues),"invalid run status"));
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
raw = struct("perSpeed",table(5,1,0,0, ...
    'VariableNames',{'vCar','exitflag','max_ceq', ...
    'max_inequality_violation'}));
run = rampSpeed.normalizeRampResult(raw,"lateral", ...
    struct("speeds",[5 10]),fixture.cases(1),struct());

verifyEqual(testCase,run.perSpeed.speed_mps,[5;10],"AbsTol",1e-12);
verifyTrue(testCase,run.perSpeed.valid(1));
verifyFalse(testCase,run.perSpeed.valid(2));
verifyEqual(testCase,run.perSpeed.status(2),"missing");
end

function testNormalizationMatchesSparseRequestedSpeeds(testCase)
fixture = makeRampFixture();
raw = fixture.legacyRampResult;
raw.perSpeed = raw.perSpeed([1 2],:);
raw.perSpeed.vCar = [5.0000000005;15];
raw.perSpeed.gLat_top = [1.52;1.30];

run = rampSpeed.normalizeRampResult(raw,"lateral", ...
    struct("speeds",[5 10 15],"speedMatchTolerance_mps",1e-6), ...
    fixture.cases(1),struct());

verifyEqual(testCase,run.perSpeed.speed_mps,[5;10;15],"AbsTol",1e-12);
verifyEqual(testCase,run.perSpeed.valid,[true;false;true]);
verifyEqual(testCase,run.perSpeed.aLat_sustainable_mps2([1 3]), ...
    [1.52;1.30]*9.80665,"AbsTol",1e-12);
verifyTrue(testCase,isnan(run.perSpeed.aLat_sustainable_mps2(2)));
verifyEqual(testCase,run.perSpeed.status([1 3]),["complete";"complete"]);
verifyEqual(testCase,run.perSpeed.status(2),"missing");
verifyThat(testCase,run.perSpeed.reason(2), ...
    matlab.unittest.constraints.ContainsSubstring("requested speed"));
end

function testNormalizationKeepsRequestedGridWhenRawHasExtraSpeed(testCase)
fixture = makeRampFixture();
raw = fixture.legacyRampResult;
raw.perSpeed.vCar = [5;15];

run = rampSpeed.normalizeRampResult(raw,"lateral", ...
    struct("speeds",5),fixture.cases(1),struct());

verifyEqual(testCase,height(run.perSpeed),1);
verifyEqual(testCase,run.perSpeed.speed_mps,5,"AbsTol",1e-12);
verifyTrue(testCase,run.perSpeed.valid);
verifyTrue(testCase,any(contains(string(run.runMeta.warnings), ...
    "unused raw speed row(s): 15")));
end

function testNormalizationMapsLongitudinalLegacyAliases(testCase)
fixture = makeRampFixture();
raw = struct();
raw.perSpeed = table(12,2.5,0,0,0,0,0.08,0.10,1,0.001,0, ...
    'VariableNames',{'long_vel','long_accel','lat_accel','steer_angle', ...
    'lat_vel','yaw_rate','kappa_3','kappa_4','exitflag','max_ceq', ...
    'max_inequality_violation'});
raw.points = raw.perSpeed(:,{'long_vel','long_accel','lat_accel', ...
    'steer_angle','lat_vel','yaw_rate','kappa_3','kappa_4', ...
    'exitflag','max_ceq','max_inequality_violation'});

run = rampSpeed.normalizeRampResult(raw,"longitudinal", ...
    struct("speeds",12,"ceqTol",1e-2),fixture.cases(2),struct());

verifyEqual(testCase,run.perSpeed.speed_mps,12,"AbsTol",1e-12);
verifyEqual(testCase,run.perSpeed.aLong_max_mps2,2.5,"AbsTol",1e-12);
verifyEqual(testCase,run.perSpeed.aLong_mps2,2.5,"AbsTol",1e-12);
verifyEqual(testCase,run.perSpeed.aLat_achieved_mps2,0,"AbsTol",1e-12);
verifyEqual(testCase,run.perSpeed.rear_slip_ratio,0.09,"AbsTol",1e-12);
verifyTrue(testCase,run.perSpeed.pure_ay0);
verifyTrue(testCase,run.perSpeed.steer_zero);
verifyTrue(testCase,run.perSpeed.lat_velocity_zero);
verifyTrue(testCase,run.perSpeed.yaw_rate_zero);
verifyTrue(testCase,run.perSpeed.valid);
verifyEqual(testCase,run.points.speed_mps,12,"AbsTol",1e-12);
verifyEqual(testCase,run.points.aLong_mps2,2.5,"AbsTol",1e-12);
verifyTrue(testCase,run.points.pure_ay0);
verifyTrue(testCase,run.points.steer_zero);
verifyTrue(testCase,run.points.lat_velocity_zero);
verifyTrue(testCase,run.points.yaw_rate_zero);
verifyEqual(testCase,run.points.kappa_RL,0.08,"AbsTol",1e-12);
verifyEqual(testCase,run.points.kappa_RR,0.10,"AbsTol",1e-12);
end

function testPointNormalizationConvertsLegacyAnglesAndForces(testCase)
fixture = makeRampFixture();
raw = fixture.legacyRampResult;
run = rampSpeed.normalizeRampResult(raw,"lateral",raw.settings, ...
    fixture.cases(1),struct());

verifyEqual(testCase,run.points.alpha_FL_rad(1),pi/180,"AbsTol",1e-12);
verifyEqual(testCase,run.points.gamma_FL_rad(1),0.01*pi/180,"AbsTol",1e-12);
verifyEqual(testCase,run.points.aero_residual_m(1), ...
    raw.points.aero_residual_in(1)*0.0254,"AbsTol",1e-12);
verifyEqual(testCase,run.points.Fz_FL_N(1),raw.points.Fz_1(1),"AbsTol",1e-12);
verifyTrue(testCase,all(run.points.valid));
end

function testNormalizationRequiresFeasibilityEvidence(testCase)
fixture = makeRampFixture();
bad = struct("perSpeed",table(5,-2,4,0, ...
    'VariableNames',{'vCar','exitflag','max_ceq', ...
    'max_inequality_violation'}));
badRun = rampSpeed.normalizeRampResult(bad,"lateral", ...
    struct("speeds",5,"ceqTol",1e-2),fixture.cases(1),struct());
verifyFalse(testCase,badRun.perSpeed.valid);
verifyEqual(testCase,badRun.perSpeed.status,"invalid");

good = bad;
good.perSpeed.exitflag = 2;
good.perSpeed.max_ceq = 1e-3;
goodRun = rampSpeed.normalizeRampResult(good,"lateral", ...
    struct("speeds",5,"ceqTol",1e-2),fixture.cases(1),struct());
verifyTrue(testCase,goodRun.perSpeed.valid);
verifyEqual(testCase,goodRun.perSpeed.status,"complete");

unknown = struct("perSpeed",table(5,1, ...
    'VariableNames',{'vCar','exitflag'}));
unknownRun = rampSpeed.normalizeRampResult(unknown,"lateral", ...
    struct("speeds",5),fixture.cases(1),struct());
verifyFalse(testCase,unknownRun.perSpeed.valid);
verifyEqual(testCase,unknownRun.perSpeed.status,"unknown");
end

function testPointEqualityResidualGatesValidity(testCase)
fixture = makeRampFixture();
raw = struct();
raw.perSpeed = table(5,1,0,0, ...
    'VariableNames',{'vCar','exitflag','max_ceq', ...
    'max_inequality_violation'});
raw.points = table(5,1,1e-3, ...
    'VariableNames',{'vCar','exitflag','max_equality_residual'});

accepted = rampSpeed.normalizeRampResult(raw,"lateral", ...
    struct("speeds",5,"ceqTol",1e-2),fixture.cases(1),struct());
verifyTrue(testCase,accepted.points.valid);
verifyEqual(testCase,accepted.points.status,"complete");

raw.points.max_equality_residual = 1;
rejected = rampSpeed.normalizeRampResult(raw,"lateral", ...
    struct("speeds",5,"ceqTol",1e-2),fixture.cases(1),struct());
verifyFalse(testCase,rejected.points.valid);
verifyEqual(testCase,rejected.points.status,"invalid");

raw.points = table(5,'VariableNames',{'vCar'});
unknown = rampSpeed.normalizeRampResult(raw,"lateral", ...
    struct("speeds",5),fixture.cases(1),struct());
verifyFalse(testCase,unknown.points.valid);
verifyEqual(testCase,unknown.points.status,"unknown");
end

function testNegativeMinFzSetsWheelLiftAndBoundInvalidity(testCase)
fixture = makeRampFixture();
raw = struct();
raw.perSpeed = table(5,1,0,0, ...
    'VariableNames',{'vCar','exitflag','max_ceq', ...
    'max_inequality_violation'});
raw.points = table(5,-2,1,0, ...
    'VariableNames',{'vCar','min_Fz','exitflag', ...
    'max_constraint_residual'});

run = rampSpeed.normalizeRampResult(raw,"lateral", ...
    struct("speeds",5,"ceqTol",1e-2),fixture.cases(1),struct());

verifyEqual(testCase,run.points.min_Fz_N,-2,"AbsTol",1e-12);
verifyTrue(testCase,run.points.wheel_lift);
verifyEqual(testCase,run.points.max_inequality_violation,2,"AbsTol",1e-12);
verifyFalse(testCase,run.points.valid);
verifyEqual(testCase,run.points.status,"invalid");
end

function testValidationRejectsLegacyPointSpellings(testCase)
study = makeValidStudy();
study.runs(1).points.alpha_1 = ones(height(study.runs(1).points),1);
study.runs(1).points.Fz_1 = ones(height(study.runs(1).points),1);
study.runs(1).points.omega_1 = ones(height(study.runs(1).points),1);

[ok,issues] = rampSpeed.validateStudy(study);

verifyFalse(testCase,ok);
verifyGreaterThanOrEqual(testCase,nnz(contains(string(issues),"non-SI")),1);
end

function testFixtureRunIsDeterministic(testCase)
caseInfo = struct("id","baseline","label","baseline", ...
    "carRole","lap");
first = makeFixtureRun(caseInfo);
second = makeFixtureRun(caseInfo);

verifyEqual(testCase,first,second);
end

function testMissingFeasibilityFieldsRemainUnknown(testCase)
fixture = makeRampFixture();
raw = struct("perSpeed",table(5,'VariableNames',{'vCar'}));
run = rampSpeed.normalizeRampResult(raw,"lateral", ...
    struct("speeds",5),fixture.cases(1),struct());

verifyFalse(testCase,run.perSpeed.valid);
verifyEqual(testCase,run.perSpeed.status,"unknown");
verifyNotEmpty(testCase,run.perSpeed.reason);
end

function study = makeValidStudy()
run = makeFixtureRun(struct("id","baseline","label","baseline", ...
    "carRole","lap"));
study = rampSpeed.makeStudy("test-app");
study.cases = struct("id","baseline","label","baseline", ...
    "source","fixture","designRow",1,"carRole","lap");
study.runs = run;
end
