# Visual Ramp-Speed Analysis App Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (- [ ]) syntax for tracking.

**Goal:** Build a MATLAB App Designer application that runs lateral-limit and pure-longitudinal ramp studies, normalizes their results, overlays multiple vehicle setups, and exposes capability, balance, aero, suspension, tire, and solver-health metrics.

**Architecture:** Keep the existing vehicle model and solvers as the numerical source of truth. Add a +rampSpeed package containing pure schema/cache/runner/plot-data functions, thin adapters around rampSweep and max_long_accel, and a thin RampSpeedApp.mlapp shell that owns controls, background execution, and rendering.

**Tech Stack:** MATLAB App Designer, MATLAB tables/structs, Optimization Toolbox, optional Parallel Computing Toolbox (parfeval, DataQueue), existing Car/aero/tire code, MATLAB unit tests.

**Spec:** docs/superpowers/specs/2026-09-19-ramp-speed-app-design.md

## Global Constraints

- Preserve all existing user changes in C:\VD; never use blanket git add, git clean, git reset, or checkout operations.
- App-specific commits must stage only the exact files named by the task.
- Do not modify Car.m, Tire2.m, carConfig.m, or the tire-model files for this feature.
- Existing rampSweep, max_long_accel, long_accel_sweep, plotRampSpeedStudy, and current tests must keep their existing default call behavior.
- The app stores canonical results in SI units: metres, seconds, Newtons, radians, and m/s^2; display conversion happens only in plot/table builders.
- The lateral mode is coast by default for compatibility with runRampSpeedStudy.m; balanced remains selectable and visible in run metadata.
- Pure longitudinal mode means steer = 0, lateral velocity = 0, yaw rate = 0, and achieved Ay = 0; lateral-only metrics are unavailable, never numeric zero.
- A solver exit flag is not sufficient for validity; feasibility residuals and bound violations are recorded and gate plotting/interpolation.
- Use TDD: each behavioral change starts with a failing MATLAB test, then the smallest implementation, then focused and full verification.
- Keep solver functions UI-independent. UI callbacks may consume progress events but solver/adaptor code must not reference matlab.apps.AppBase.
- Use multiple Luna Max agents only on file-disjoint tasks; the integration pass owns cross-task conflicts and final verification.

---

### Task 1: Establish the versioned result schema and deterministic fixtures

**Files:**
- Create: Full Car Models/+rampSpeed/makeStudy.m
- Create: Full Car Models/+rampSpeed/makeRun.m
- Create: Full Car Models/+rampSpeed/normalizeRampResult.m
- Create: Full Car Models/+rampSpeed/validateStudy.m
- Create: Full Car Models/+rampSpeed/displayUnits.m
- Create: Full Car Models/tests/helpers/makeRampFixture.m
- Create: Full Car Models/tests/helpers/makeFixtureRun.m
- Create: Full Car Models/tests/test_rampSpeedSchema.m
- Create: Full Car Models/tests/test_rampSpeedValidation.m

**Interfaces:**
- Produces study = rampSpeed.makeStudy(appVersion).
- Produces run = rampSpeed.makeRun(type,mode,settings,caseInfo).
- Produces run = rampSpeed.normalizeRampResult(raw,type,settings,caseInfo,runMeta).
- Produces [ok,issues] = rampSpeed.validateStudy(study).
- Produces display = rampSpeed.displayUnits(values,fromUnits,toUnits) for presentation-only conversion.
- The normalized run always contains schemaVersion, caseId, type, mode, settings, perSpeed, points, runMeta, status, and raw.

**Canonical fields:**

The makeRampFixture helper returns fixture.cars, fixture.cases, fixture.frontRideHeightIn, and fixture.legacyRampResult. The makeFixtureRun helper returns a complete normalized run for a supplied caseInfo and is used by runner tests as the injected fake case function.

- Common per-speed fields: speed_mps, valid, status, reason, aLat_mps2, aLong_mps2, engine_rpm, current_gear, throttle, downforce_N, drag_N, ClA_m2, CdA_m2, LoD, aero_balance_front, aero_outside_map, aero_residual_m, Fz_front_axle_N, Fz_rear_axle_N, min_Fz_N, wheel_lift, LLTD, LLT_front_N, LLT_rear_N, and long_load_transfer_N.
- Lateral fields: aLat_free_mps2, aLat_sustainable_mps2, ramp_complete_fraction, truncated, power_limited, K_linear_rad_per_mps2, K_linear_r2, K_at_limit_rad_per_mps2, cuo_steer_linear_rad, cuo_steer_limit_rad, mechanical_balance_front, grip_balance_mid, grip_balance_limit, alpha_balance_mid_rad, alpha_balance_limit_rad, LLT_norm_balance_mid, LLT_norm_balance_limit, front_Fz_fraction_mid, front_Fz_fraction_limit, front_downforce_N, rear_downforce_N, front_ride_height_m, rear_ride_height_m, pitch_rad, front_shock_travel_m, rear_shock_travel_m, front_camber_rad, rear_camber_rad, min_Fz_limit_N, max_constraint_residual, n_exitflag1, and n_exitflag2.
- Longitudinal fields: aLong_max_mps2, aLat_achieved_mps2, aLat_force_residual_mps2, pure_ay0, steer_zero, lat_velocity_zero, yaw_rate_zero, rear_slip_ratio, throttle_upper_active, rear_slip_upper_active, power_limited, traction_limited, and lateral_metrics_applicable.
- Detailed point fields use speed_mps, speed_index, point_index, valid, status, exitflag, max_constraint_residual, max_equality_residual, max_inequality_violation, plus common forces and per-corner fields Fz_FL_N, Fz_FR_N, Fz_RL_N, Fz_RR_N, Fx_*_N, Fy_*_N, alpha_*_rad, gamma_*_rad, kappa_*, T_*_Nm, and omega_*_rps.
- raw retains the original rampSweep result or longitudinal diagnostic payload for debugging and legacy plotting; canonical app data never relies on legacy inch/degree/g fields.

- [ ] **Step 1: Write the failing schema-construction test**

Create test_rampSpeedSchema.m with a synthetic run and assert the exact top-level fields, types, SI names, and type-specific availability.

~~~matlab
function tests = test_rampSpeedSchema
tests = functiontests(localfunctions);
end

function testCreatesVersionedLateralRun(testCase)
fixture = makeRampFixture();
run = rampSpeed.makeRun("lateral","coast",struct("speeds",[5 10]), ...
    struct("id","baseline","label","baseline","carRole","lap"));
verifyEqual(testCase,run.schemaVersion,1);
verifyEqual(testCase,run.type,"lateral");
verifyTrue(testCase,all(ismember(["perSpeed","points","runMeta","status","raw"], ...
    string(fieldnames(run)))));
verifyTrue(testCase,all(ismember(["speed_mps","aLat_free_mps2", ...
    "aLat_sustainable_mps2","mechanical_balance_front", ...
    "aero_balance_front"], ...
    string(run.perSpeed.Properties.VariableNames))));
verifyEqual(testCase,run.perSpeed.front_ride_height_m(1), ...
    fixture.frontRideHeightIn*0.0254,"AbsTol",1e-12);
end

function testLongitudinalOnlyFieldsAreUnavailable(testCase)
run = rampSpeed.makeRun("longitudinal","",struct("speeds",5), ...
    struct("id","accel","label","accel","carRole","acceleration"));
verifyFalse(testCase,run.runMeta.lateralMetricsApplicable);
verifyTrue(testCase,all(isnan(run.perSpeed.K_linear_rad_per_mps2)));
end
~~~

- [ ] **Step 2: Run the focused test and verify it fails**

Run:

~~~powershell
matlab -batch "cd('C:\VD\Full Car Models'); setup_paths; r=runtests('tests/test_rampSpeedSchema.m'); assert(any([r.Failed]));"
~~~

Expected: FAIL because the rampSpeed package and fixture helpers do not exist.

- [ ] **Step 3: Implement the schema constructors and fixture**

Implement makeStudy with schemaVersion = 1, created = datetime('now'), appVersion, empty cases/runs, and a displayUnits record. Implement makeRun with empty typed tables and explicit status = "pending". Implement makeRampFixture with two setup labels, a lateral run, a longitudinal run, and representative raw legacy fields.

- [ ] **Step 4: Implement normalization and validation**

Implement normalizeRampResult as a pure mapping function. Convert inches to metres, degrees to radians, g-based accelerations to m/s^2, and K_linear from deg/g to rad/(m/s^2). Preserve missing fields as NaN plus a reason/validity flag. Implement validateStudy to return ok = false with named issues for missing schema fields, unsupported type, duplicate case IDs, invalid status, or non-SI canonical column names.

- [ ] **Step 5: Run the focused tests and verify they pass**

Run:

~~~powershell
matlab -batch "cd('C:\VD\Full Car Models'); setup_paths; r=runtests({'tests/test_rampSpeedSchema.m','tests/test_rampSpeedValidation.m'}); assert(all(~[r.Failed]) && all(~[r.Incomplete]));"
~~~

Expected: PASS with all new schema and validation tests green.

- [ ] **Step 6: Commit only the schema files**

~~~powershell
git -c safe.directory='C:/VD' -C 'C:\VD' add -- "Full Car Models/+rampSpeed/makeStudy.m" "Full Car Models/+rampSpeed/makeRun.m" "Full Car Models/+rampSpeed/normalizeRampResult.m" "Full Car Models/+rampSpeed/validateStudy.m" "Full Car Models/+rampSpeed/displayUnits.m" "Full Car Models/tests/helpers/makeRampFixture.m" "Full Car Models/tests/helpers/makeFixtureRun.m" "Full Car Models/tests/test_rampSpeedSchema.m" "Full Car Models/tests/test_rampSpeedValidation.m"
git -c safe.directory='C:/VD' -C 'C:\VD' commit -m "feat: add ramp speed result schema"
~~~

---

### Task 2: Add legacy-cache migration and atomic study persistence

**Files:**
- Create: Full Car Models/+rampSpeed/migrateLegacyStudy.m
- Create: Full Car Models/+rampSpeed/loadStudy.m
- Create: Full Car Models/+rampSpeed/saveStudy.m
- Create: Full Car Models/tests/test_rampStudyMigration.m
- Create: Full Car Models/tests/test_rampStudyRoundTrip.m

**Interfaces:**
- study = rampSpeed.migrateLegacyStudy(fileName,appVersion)
- study = rampSpeed.loadStudy(fileName,appVersion)
- rampSpeed.saveStudy(fileName,study)

- [ ] **Step 1: Write the failing migration and round-trip tests**

Construct a temporary legacy file containing:

~~~matlab
study.results = {legacyRampResult};
study.labels = "baseline";
study.rampOptions = struct("speeds",[5 10],"nRamp",4,"mode","coast");
study.numWorkers = 0;
study.created = datetime("now");
save(legacyFile,"study","-v7.3");
~~~

Assert that migration creates schemaVersion = 1, one case, one run, canonical SI fields, preserved labels/settings, and no change to the source file timestamp or bytes. Assert that a current-schema study saves and loads with equal case IDs, status, and per-speed values.

- [ ] **Step 2: Run the migration tests and verify they fail**

Run:

~~~powershell
matlab -batch "cd('C:\VD\Full Car Models'); setup_paths; r=runtests({'tests/test_rampStudyMigration.m','tests/test_rampStudyRoundTrip.m'}); assert(any([r.Failed]));"
~~~

Expected: FAIL because the cache functions do not exist.

- [ ] **Step 3: Implement migrateLegacyStudy**

Load only the study variable from the source file. Map study.results{i} to study.runs(i), study.labels{i} to study.cases(i).label, retain rampOptions, numWorkers, and created in runMeta.legacy, and call normalizeRampResult for every result. Insert invalid rows for requested speeds omitted by the legacy rampSweep; do not interpolate them. Never save to or overwrite the source filename.

- [ ] **Step 4: Implement loadStudy and saveStudy**

loadStudy must accept current schema files and route legacy study.results files through migration. Reject schemaVersion > 1 with identifier rampSpeed:unsupportedSchema. saveStudy must validate before saving, write to tempname(fileparts(fileName)) with -v7.3, then replace the requested target with movefile(...,'f') only after a successful save. Record the final absolute path in study.runMeta.fileName.

- [ ] **Step 5: Run the migration and round-trip tests**

Run:

~~~powershell
matlab -batch "cd('C:\VD\Full Car Models'); setup_paths; r=runtests({'tests/test_rampStudyMigration.m','tests/test_rampStudyRoundTrip.m'}); assert(all(~[r.Failed]) && all(~[r.Incomplete]));"
~~~

Expected: PASS, with the legacy source file unchanged and current-schema studies round-tripping.

- [ ] **Step 6: Commit only the cache files**

~~~powershell
git -c safe.directory='C:/VD' -C 'C:\VD' add -- "Full Car Models/+rampSpeed/migrateLegacyStudy.m" "Full Car Models/+rampSpeed/loadStudy.m" "Full Car Models/+rampSpeed/saveStudy.m" "Full Car Models/tests/test_rampStudyMigration.m" "Full Car Models/tests/test_rampStudyRoundTrip.m"
git -c safe.directory='C:/VD' -C 'C:\VD' commit -m "feat: migrate and persist ramp studies"
~~~

---

### Task 3: Wrap the lateral solver without changing its defaults

**Files:**
- Create: Full Car Models/+rampSpeed/runLateralRamp.m
- Create: Full Car Models/tests/test_lateralRampAdapter.m
- Modify: Full Car Models/sweeps/rampSweep.m
- Modify: Full Car Models/tests/test_rampSweepAdditionalMetrics.m

**Interfaces:**
- run = rampSpeed.runLateralRamp(car,settings,caseInfo,callbacks)
- callbacks.onProgress(event) is optional.
- callbacks.isCancelled() is optional and returns a scalar logical.
- settings contains the existing rampSweep options plus speeds, mode, nRamp, nBisect, ceqTol, ayMinFrac, ayMaxFrac, linearFrac, and verbose.

- [ ] **Step 1: Write the failing adapter equivalence test**

Use a small deterministic case and compare the adapter raw result against a direct call:

~~~matlab
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
end
~~~

Add a second test with a callback that returns true after the first speed and assert run.status = "cancelled" with completed rows retained.

- [ ] **Step 2: Run the lateral adapter tests and verify they fail**

Run:

~~~powershell
matlab -batch "cd('C:\VD\Full Car Models'); setup_paths; r=runtests('tests/test_lateralRampAdapter.m'); assert(any([r.Failed]));"
~~~

Expected: FAIL because the adapter and callback hooks do not exist.

- [ ] **Step 3: Add optional progress/cancel hooks to rampSweep**

Keep all existing defaults and the one-output signature. Add empty defaults for opts.progressFcn and opts.cancelFcn. Add speed_index = iv to each raw per-speed and point record. At the start of each speed and after each completed speed, call:

~~~matlab
if ~isempty(cancelFcn) && cancelFcn()
    cancelled = true;
    break
end
if ~isempty(progressFcn)
    progressFcn(struct("phase","speed","speedIndex",iv, ...
        "speed_mps",v,"completedSpeeds",iv-1, ...
        "requestedSpeeds",numel(speeds)));
end
~~~

When cancellation occurs after at least one solved speed, return the partial result with R.status = "cancelled"; preserve the existing rampSpeed:noSolution error for the default no-solution path when no callback is supplied.

- [ ] **Step 4: Implement runLateralRamp**

Build raw options from settings, inject the callback hooks, call rampSweep, catch only the explicit cancellation/no-solution cases, and call normalizeRampResult. Reconstruct a full requested-speed perSpeed table using speed_index; a skipped speed receives valid=false, status="failed", and its captured error reason. Preserve gLat_max, gLat_top, ramp_complete, power_limited, max_ceq, exit counts, front/rear aero loads, balance metrics, ride height, camber, pitch, shock travel, and wheel-lift fields.

- [ ] **Step 5: Run focused lateral and regression tests**

Run:

~~~powershell
matlab -batch "cd('C:\VD\Full Car Models'); setup_paths; r=runtests({'tests/test_lateralRampAdapter.m','tests/test_rampSweepAdditionalMetrics.m','tests/test_plotRampSpeedStudy.m'}); assert(all(~[r.Failed]) && all(~[r.Incomplete]));"
~~~

Expected: PASS, including existing raw-field and plot behavior.

- [ ] **Step 6: Commit only the lateral files**

~~~powershell
git -c safe.directory='C:/VD' -C 'C:\VD' add -- "Full Car Models/+rampSpeed/runLateralRamp.m" "Full Car Models/sweeps/rampSweep.m" "Full Car Models/tests/test_lateralRampAdapter.m" "Full Car Models/tests/test_rampSweepAdditionalMetrics.m"
git -c safe.directory='C:/VD' -C 'C:\VD' commit -m "feat: add lateral ramp adapter"
~~~

---

### Task 4: Add the pure-longitudinal adapter and diagnostics

**Files:**
- Create: Full Car Models/+rampSpeed/runLongitudinalRamp.m
- Create: Full Car Models/tests/test_longitudinalRampAdapter.m
- Modify: Full Car Models/events/max_long_accel.m

**Interfaces:**
- Preserve existing calls: [x_accel,long_accel,long_accel_guess] = max_long_accel(vCar,car,x0).
- Add optional calls: [x_accel,long_accel,long_accel_guess,diagnostics] = max_long_accel(vCar,car,x0,solverOptions).
- Produce run = rampSpeed.runLongitudinalRamp(car,settings,caseInfo,callbacks).

- [ ] **Step 1: Write the failing pure-Ay=0 test**

~~~matlab
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
verifyTrue(testCase,isnan(row.K_linear_rad_per_mps2));
end
~~~

- [ ] **Step 2: Run the test and verify it fails**

Run:

~~~powershell
matlab -batch "cd('C:\VD\Full Car Models'); setup_paths; r=runtests('tests/test_longitudinalRampAdapter.m'); assert(any([r.Failed]));"
~~~

Expected: FAIL because the fourth diagnostic output and adapter do not exist.

- [ ] **Step 3: Add optional diagnostics to max_long_accel**

Keep the current three outputs and solver defaults. Accept an optional solverOptions struct with maxFunctionEvaluations, constraintTolerance, stepTolerance, and display. After the existing solve, populate:

~~~matlab
diagnostics.state = x;
diagnostics.exitflag = exitflag;
[diagnostics.c,diagnostics.ceq] = car.constraint1(x);
diagnostics.max_inequality_violation = max([diagnostics.c(:);0]);
diagnostics.max_equality_residual = max(abs(diagnostics.ceq(:)));
diagnostics.metrics = car.metrics(x);
diagnostics.pure_ay0 = abs(diagnostics.metrics.gLat) <= 1e-12;
~~~

Do not change x_accel, long_accel, or the third output used by long_accel_sweep.

- [ ] **Step 4: Implement runLongitudinalRamp**

Loop over the explicit requested speed vector rather than calling long_accel_sweep unchanged. Warm-start from the prior state, call max_long_accel, use diagnostics.state with Car.metrics, and create one normalized row per requested speed. Record steer_zero, lat_velocity_zero, yaw_rate_zero, pure_ay0, residuals, throttle/rear-slip bounds, aero/load/suspension/tire fields, and explicit power_limited/traction_limited indicators. Insert invalid rows for failed speeds and keep lateral_metrics_applicable=false.

- [ ] **Step 5: Run the adapter and backward-compatibility tests**

Run:

~~~powershell
matlab -batch "cd('C:\VD\Full Car Models'); setup_paths; r=runtests({'tests/test_longitudinalRampAdapter.m','tests/test_rampSweepAdditionalMetrics.m'}); assert(all(~[r.Failed]) && all(~[r.Incomplete]));"
~~~

Expected: PASS, including existing three-output callers.

- [ ] **Step 6: Commit only the longitudinal files**

~~~powershell
git -c safe.directory='C:/VD' -C 'C:\VD' add -- "Full Car Models/+rampSpeed/runLongitudinalRamp.m" "Full Car Models/events/max_long_accel.m" "Full Car Models/tests/test_longitudinalRampAdapter.m"
git -c safe.directory='C:/VD' -C 'C:\VD' commit -m "feat: add pure longitudinal ramp adapter"
~~~

---

### Task 5: Add the study runner, setup catalog, progress, cancellation, and checkpoints

**Files:**
- Create: Full Car Models/+rampSpeed/setupCaseCatalog.m
- Create: Full Car Models/+rampSpeed/runStudy.m
- Create: Full Car Models/+rampSpeed/runCase.m
- Create: Full Car Models/tests/test_rampStudyRunner.m
- Create: Full Car Models/tests/test_rampStudyCancellation.m

**Interfaces:**
- cases = rampSpeed.setupCaseCatalog(carCell,designTable,carRole)
- run = rampSpeed.runCase(car,caseInfo,request,callbacks)
- [study,events] = rampSpeed.runStudy(cars,cases,request,callbacks)
- request fields: rampType, settings, parallelRequested, numWorkers, checkpointPath, appVersion, runCaseFcn, and progressQueue.

- [ ] **Step 1: Write the failing serial runner test**

Use makeRampFixture and inject request.runCaseFcn so the runner is testable without a solver:

~~~matlab
function testRunnerRetainsCaseMetadata(testCase)
fixture = makeRampFixture();
cases = fixture.cases;
request = struct("rampType","lateral","settings",struct(), ...
    "parallelRequested",false,"numWorkers",0, ...
    "checkpointPath","", "appVersion","test", ...
    "runCaseFcn",@(car,caseInfo,request,callbacks) ...
        makeFixtureRun(caseInfo));
[study,events] = rampSpeed.runStudy(fixture.cars,cases,request,struct());
verifyEqual(testCase,numel(study.runs),2);
verifyEqual(testCase,study.runs(1).caseId,cases(1).id);
verifyEqual(testCase,study.runs(1).status,"complete");
verifyGreaterThanOrEqual(testCase,numel(events),2);
end
~~~

- [ ] **Step 2: Write the failing cancellation/partial-cache test**

Inject a fake case runner that returns one complete case, then have callbacks.isCancelled return true. Assert that the first case remains complete, the study is cancelled, the checkpoint file exists, and UI-facing events contain the cancellation reason.

- [ ] **Step 3: Run the runner tests and verify they fail**

Run:

~~~powershell
matlab -batch "cd('C:\VD\Full Car Models'); setup_paths; r=runtests({'tests/test_rampStudyRunner.m','tests/test_rampStudyCancellation.m'}); assert(any([r.Failed]));"
~~~

Expected: FAIL because the runner functions do not exist.

- [ ] **Step 4: Implement setup selection and serial execution**

Map carConfig outputs into case records with stable IDs, labels, design row, source index, and carRole. Use carCell(:,1) for lateral lap cars and carCell(:,2) for pure longitudinal acceleration cars unless the user explicitly selects a different role. runStudy must initialize all cases, execute serially in order, catch per-case errors, update status, and save after every completed case when checkpointPath is nonempty.

- [ ] **Step 5: Implement progress events and cancellation**

Use callback events shaped as:

~~~matlab
event = struct("phase","case","caseId",caseInfo.id, ...
    "speedIndex",speedIndex,"speed_mps",speed_mps, ...
    "completedCases",completedCases, ...
    "totalCases",numel(cases),"message",message);
~~~

Send events through callbacks.onProgress in serial mode or request.progressQueue when running in a worker. Check callbacks.isCancelled/the cancellation token between speeds and cases. Save partial results before returning status = "cancelled".

- [ ] **Step 6: Implement optional parallel execution and serial fallback**

Use parfeval only when parallelRequested is true and the Parallel Computing Toolbox is available. Store effectiveWorkers and parallelFallbackReason in runMeta. Without the toolbox, set effectiveWorkers = 0, run serially, and emit a visible warning. Do not run nested parfor loops.

- [ ] **Step 7: Run focused runner tests**

Run:

~~~powershell
matlab -batch "cd('C:\VD\Full Car Models'); setup_paths; r=runtests({'tests/test_rampStudyRunner.m','tests/test_rampStudyCancellation.m'}); assert(all(~[r.Failed]) && all(~[r.Incomplete]));"
~~~

Expected: PASS for serial, per-case failure retention, checkpointing, cancellation, and no-toolbox fallback.

- [ ] **Step 8: Commit only the runner files**

~~~powershell
git -c safe.directory='C:/VD' -C 'C:\VD' add -- "Full Car Models/+rampSpeed/setupCaseCatalog.m" "Full Car Models/+rampSpeed/runStudy.m" "Full Car Models/+rampSpeed/runCase.m" "Full Car Models/tests/test_rampStudyRunner.m" "Full Car Models/tests/test_rampStudyCancellation.m"
git -c safe.directory='C:/VD' -C 'C:\VD' commit -m "feat: add ramp study runner"
~~~

---

### Task 6: Centralize metric plotting, overlays, deltas, and inspector data

**Files:**
- Create: Full Car Models/+rampSpeed/metricCatalog.m
- Create: Full Car Models/+rampSpeed/buildPlotData.m
- Create: Full Car Models/+rampSpeed/buildComparison.m
- Create: Full Car Models/+rampSpeed/buildInspectorTable.m
- Create: Full Car Models/+rampSpeed/renderMetric.m
- Create: Full Car Models/tests/test_rampSpeedPlotCatalog.m
- Create: Full Car Models/tests/test_rampSpeedComparison.m

**Interfaces:**
- catalog = rampSpeed.metricCatalog()
- data = rampSpeed.buildPlotData(runs,metricId,options)
- delta = rampSpeed.buildComparison(runs,metricId,baselineId,options)
- T = rampSpeed.buildInspectorTable(run,selection,units)
- h = rampSpeed.renderMetric(ax,data,options)

- [ ] **Step 1: Write failing catalog and overlay tests**

Create fixture runs with two setup labels, different speed grids, one invalid speed, and one truncated lateral row. Assert that the catalog contains the existing IDs aero_front_load, aero_rear_load, aero_balance, mechanical_balance, handling_balance, ride-height/camber/pitch/shock entries, plus capability, drag, downforce, wheel-load, validity, and raw-ramp entries.

- [ ] **Step 2: Run the plot tests and verify they fail**

Run:

~~~powershell
matlab -batch "cd('C:\VD\Full Car Models'); setup_paths; r=runtests({'tests/test_rampSpeedPlotCatalog.m','tests/test_rampSpeedComparison.m'}); assert(any([r.Failed]));"
~~~

Expected: FAIL because the plot-data functions do not exist.

- [ ] **Step 3: Implement the catalog**

Each catalog row must define id, tab, validTypes, sourceLevel, field or a derivation function, stored units, display units, scale, y-label, title, subtitle, zero-line flag, series style, and validity rule. Preserve current plot conventions: vCar (m/s), setPlotFont, lines() for setup colors, parula() for speed-colored raw curves, hollow orange markers for truncated ramps, and hidden zero-reference lines.

- [ ] **Step 4: Implement native-grid plot data**

Build one series per setup from native speed_mps values. Use NaN for invalid points so plots show gaps, carry valid, truncated, power_limited, wheel_lift, and reason masks for markers/tooltips, and filter catalog entries by validTypes. Do not use exact floating-point equality to group raw points; use speed_index.

- [ ] **Step 5: Implement comparison deltas**

buildComparison accepts an explicit comparison grid or constructs the common finite speed domain. Use linear interpolation only within contiguous valid intervals, never across invalid/truncated gaps, and never extrapolate. Compute variant - baseline, store baselineId, comparisonGrid_mps, and interpolationMethod, and return invalid rows when either source is invalid.

- [ ] **Step 6: Implement inspector table construction and rendering**

Return named SI and converted display columns for the selected setup, speed, and ramp point, including four-wheel forces/loads/slips/camber, residuals, aero status, and validity flags. Render line series, flags, legends, labels, zero lines, and warning markers without reading solver objects.

- [ ] **Step 7: Run the focused plot tests**

Run:

~~~powershell
matlab -batch "cd('C:\VD\Full Car Models'); setup_paths; r=runtests({'tests/test_rampSpeedPlotCatalog.m','tests/test_rampSpeedComparison.m','tests/test_plotRampSpeedStudy.m'}); assert(all(~[r.Failed]) && all(~[r.Incomplete]));"
~~~

Expected: PASS for catalog coverage, native-grid overlays, invalid gaps, truncation flags, baseline deltas, and legacy plot behavior.

- [ ] **Step 8: Commit only the plotting files**

~~~powershell
git -c safe.directory='C:/VD' -C 'C:\VD' add -- "Full Car Models/+rampSpeed/metricCatalog.m" "Full Car Models/+rampSpeed/buildPlotData.m" "Full Car Models/+rampSpeed/buildComparison.m" "Full Car Models/+rampSpeed/buildInspectorTable.m" "Full Car Models/+rampSpeed/renderMetric.m" "Full Car Models/tests/test_rampSpeedPlotCatalog.m" "Full Car Models/tests/test_rampSpeedComparison.m"
git -c safe.directory='C:/VD' -C 'C:\VD' commit -m "feat: add ramp speed plots and comparisons"
~~~

---

### Task 7: Build the App Designer shell and connect asynchronous execution

**Files:**
- Create: Full Car Models/RampSpeedApp.mlapp
- Create: Full Car Models/tests/test_rampSpeedAppSmoke.m
- Create: Full Car Models/tests/test_rampSpeedAppRequest.m

**Interfaces:**
- App callbacks call rampSpeed.setupCaseCatalog, rampSpeed.runStudy, rampSpeed.loadStudy, rampSpeed.saveStudy, rampSpeed.buildPlotData, rampSpeed.buildComparison, and rampSpeed.buildInspectorTable.
- The app owns Future, ProgressQueue, CurrentStudy, SelectedCaseIds, BaselineCaseId, and DisplayUnits.
- The app exposes request = app.buildRequestForTest() for UI tests; production callbacks call the same private request builder.

**Required UI components:**
- Setup list/table, car-role selector, ramp-type selector, lateral-mode selector, speed start/stop/step fields, explicit speed-vector field, nRamp, nBisect, tolerance fields, serial/parallel controls, unit selector, baseline selector, Run/Cancel/Load/Save/Export buttons, progress text area, and tabs for Capability, Balance, Aero & Loads, Suspension, Raw Ramp, and Inspector/Data.
- Use a persistent left run panel and tabbed right analysis area. Keep UI callbacks thin; every calculation comes from package functions.

- [ ] **Step 1: Write the failing app-request and smoke tests**

Use a fake runCaseFcn so the UI tests do not run fmincon:

~~~matlab
function testBuildRequestUsesAccelerationCarForLongitudinal(testCase)
app = RampSpeedApp("Visible","off");
testCase.addTeardown(@() delete(app));
app.RampTypeDropDown.Value = "Pure longitudinal";
app.CarRoleDropDown.Value = "Acceleration car";
app.SpeedStartField.Value = 5;
app.SpeedStopField.Value = 10;
app.SpeedStepField.Value = 5;
request = app.buildRequestForTest();
verifyEqual(testCase,request.rampType,"longitudinal");
verifyEqual(testCase,request.carRole,"acceleration");
verifyEqual(testCase,request.settings.speeds,[5 10]);
end
~~~

Add a smoke test that launches the app, selects two fixture cases, switches between ramp types, loads a fixture study, changes baseline, and verifies the plot tabs receive data.

- [ ] **Step 2: Run the UI tests and verify they fail**

Run:

~~~powershell
matlab -batch "cd('C:\VD\Full Car Models'); setup_paths; r=runtests({'tests/test_rampSpeedAppRequest.m','tests/test_rampSpeedAppSmoke.m'}); assert(any([r.Failed]));"
~~~

Expected: FAIL because RampSpeedApp.mlapp does not exist.

- [ ] **Step 3: Create the App Designer layout**

Create the named controls and tabs. Initialize the app with Visible = "on" for normal use and accept a test-only hidden/constructor option. Populate setup cases from carConfig and display case label, design-row index, and lap/acceleration-car role. Set lateral mode default to coast, pure longitudinal default role to acceleration car, and unit display default to SI.

- [ ] **Step 4: Implement request building and plot refresh**

Implement buildRequestForTest/the production equivalent to validate speed vectors, ramp type, mode, positive speed, nRamp >= 2, nonnegative bisection count, and residual tolerance. Update the baseline selector and plot tabs from CurrentStudy through buildPlotData; disable lateral-only metric controls for longitudinal studies and display “not applicable”.

- [ ] **Step 5: Implement background execution and progress**

When a background pool is available, create a parallel.pool.DataQueue, register afterEach to update the progress panel, and run:

~~~matlab
app.Future = parfeval(backgroundPool,@rampSpeed.runStudy,1, ...
    cars,cases,request,struct("progressQueue",app.ProgressQueue));
~~~

On completion, update CurrentStudy, re-enable controls, refresh plots, and terminal-log status. Without the toolbox, run serially with a visible warning and effectiveWorkers = 0. Do not pass the app object into a worker.

- [ ] **Step 6: Implement cancellation and checkpoint reload**

Cancel the future when the user presses Cancel, then reload the most recent checkpoint and display cancelled with retained completed cases. Ensure all buttons are re-enabled in both success and error callbacks. Store the future cancellation result and checkpoint path in run metadata.

- [ ] **Step 7: Run the UI smoke suite**

Run:

~~~powershell
matlab -batch "cd('C:\VD\Full Car Models'); setup_paths; r=runtests({'tests/test_rampSpeedAppRequest.m','tests/test_rampSpeedAppSmoke.m'}); assert(all(~[r.Failed]) && all(~[r.Incomplete]));"
~~~

Expected: PASS for request validation, setup/role selection, tab refresh, baseline selection, unit selection, progress, cancellation recovery, and fake-runner smoke coverage.

- [ ] **Step 8: Commit only the app files**

~~~powershell
git -c safe.directory='C:/VD' -C 'C:\VD' add -- "Full Car Models/RampSpeedApp.mlapp" "Full Car Models/tests/test_rampSpeedAppSmoke.m" "Full Car Models/tests/test_rampSpeedAppRequest.m"
git -c safe.directory='C:/VD' -C 'C:\VD' commit -m "feat: add ramp speed analysis app"
~~~

---

### Task 8: Add export, project packaging, documentation, and release verification

**Files:**
- Create: Full Car Models/RampSpeedApp.prj
- Create: Full Car Models/README_ramp_speed_app.md
- Create: Full Car Models/+rampSpeed/exportStudy.m
- Create: Full Car Models/+rampSpeed/releaseManifest.m
- Create: Full Car Models/tests/test_rampSpeedExport.m
- Create: Full Car Models/tests/test_rampSpeedReleaseManifest.m

**Interfaces:**
- paths = rampSpeed.exportStudy(study,outputDirectory,options)
- manifest = rampSpeed.releaseManifest(rootDirectory)

- [ ] **Step 1: Write the failing export and manifest tests**

Assert that export creates a normalized CSV/table, a .mat study, and requested figures without changing canonical SI values. Assert that the manifest records MATLAB release, schema/app version, file list, Git commit/dirty status, toolbox availability, and creation time.

- [ ] **Step 2: Run the tests and verify they fail**

Run:

~~~powershell
matlab -batch "cd('C:\VD\Full Car Models'); setup_paths; r=runtests({'tests/test_rampSpeedExport.m','tests/test_rampSpeedReleaseManifest.m'}); assert(any([r.Failed]));"
~~~

Expected: FAIL because export and manifest functions do not exist.

- [ ] **Step 3: Implement export and release metadata**

Use writetable for normalized per-speed and point tables, save(...,'-v7.3') for the study, and exportgraphics for visible figures. Include units in column names or a companion metadata table. Use existing simLog only from the client/app side and guarantee a terminal log row for complete, failed, and cancelled runs.

- [ ] **Step 4: Create the MATLAB Project**

Create RampSpeedApp.prj with the project root at Full Car Models, include +rampSpeed, tests, RampSpeedApp.mlapp, and the required existing model folders, and add a startup shortcut that runs setup_paths. Document supported MATLAB/toolbox requirements and the serial fallback.

- [ ] **Step 5: Write the user README**

Document:
- how to open RampSpeedApp.prj;
- how to select lap versus acceleration-car setups;
- lateral coast versus balanced meaning;
- pure-longitudinal Ay = 0 definition;
- positive/negative understeer signs;
- downforce/drag sign conventions;
- validity, truncation, power, wheel-lift, and aero-map warnings;
- where caches/exports are written;
- focused and full MATLAB test commands.

- [ ] **Step 6: Run focused export/package tests**

Run:

~~~powershell
matlab -batch "cd('C:\VD\Full Car Models'); setup_paths; r=runtests({'tests/test_rampSpeedExport.m','tests/test_rampSpeedReleaseManifest.m'}); assert(all(~[r.Failed]) && all(~[r.Incomplete]));"
~~~

Expected: PASS with a temporary export directory cleaned up by test teardown.

- [ ] **Step 7: Run the full MATLAB verification suite and static checks**

Run:

~~~powershell
matlab -batch "cd('C:\VD\Full Car Models'); setup_paths; r=runtests('tests'); assert(all(~[r.Failed]) && all(~[r.Incomplete]));"
matlab -batch "cd('C:\VD\Full Car Models'); setup_paths; F={'RampSpeedApp.mlapp','+rampSpeed/makeStudy.m','+rampSpeed/makeRun.m','+rampSpeed/normalizeRampResult.m','+rampSpeed/validateStudy.m','+rampSpeed/migrateLegacyStudy.m','+rampSpeed/loadStudy.m','+rampSpeed/saveStudy.m','+rampSpeed/runLateralRamp.m','+rampSpeed/runLongitudinalRamp.m','+rampSpeed/runStudy.m','+rampSpeed/runCase.m','+rampSpeed/metricCatalog.m','+rampSpeed/buildPlotData.m','+rampSpeed/buildComparison.m','+rampSpeed/buildInspectorTable.m','+rampSpeed/renderMetric.m','+rampSpeed/exportStudy.m','+rampSpeed/releaseManifest.m'}; disp(checkcode(F,'-id'));"
git -c safe.directory='C:/VD' -C 'C:\VD' diff --check
git -c safe.directory='C:/VD' -C 'C:\VD' status --short
~~~

Expected: all tests pass, no MATLAB parse/static errors in app package files, no whitespace errors, and only app-specific paths appear in the final diff review.

- [ ] **Step 8: Perform representative visual QA**

Run one representative setup at speeds [5 15 25] in lateral coast, lateral balanced, and pure longitudinal modes. Confirm direct-vs-app agreement for:
- free/sustainable lateral capability and ramp_complete;
- pure Ay = 0;
- mechanical/aero balance and positive/negative understeer signs;
- front/rear downforce and drag;
- ride height, pitch/roll, camber, shock travel;
- minimum wheel load and wheel-lift markers;
- invalid gaps, power limitation, aero-map warnings, legends, units, baseline deltas, and exports.

- [ ] **Step 9: Commit only packaging/documentation files**

~~~powershell
git -c safe.directory='C:/VD' -C 'C:\VD' add -- "Full Car Models/RampSpeedApp.prj" "Full Car Models/README_ramp_speed_app.md" "Full Car Models/+rampSpeed/exportStudy.m" "Full Car Models/+rampSpeed/releaseManifest.m" "Full Car Models/tests/test_rampSpeedExport.m" "Full Car Models/tests/test_rampSpeedReleaseManifest.m"
git -c safe.directory='C:/VD' -C 'C:\VD' commit -m "docs: package ramp speed app"
~~~

---

## Parallel execution assignment

Use fresh Luna Max agents only after the user selects execution mode and after the worktree/branch decision is confirmed.

- Agent A: Tasks 1–4, schema/migration and solver adapters. Owns +rampSpeed/makeStudy.m, makeRun.m, normalizeRampResult.m, validateStudy.m, cache files, runLateralRamp.m, runLongitudinalRamp.m, rampSweep.m, max_long_accel.m, and their focused tests.
- Agent B: Task 6, plot catalog/comparison/inspector. Owns only the plot package files and plot/comparison tests. It uses fixture structs and does not touch solver files.
- Agent C: Task 5, runner/checkpoint/cancellation. Owns runner package files and runner tests; it consumes the schema interfaces from Task 1 and uses injected fake runners until Agent A’s adapters land.
- Agent D: Task 7, App Designer shell. Starts after the runner and plot interfaces are reviewed; owns RampSpeedApp.mlapp and UI tests, and does not edit solver/plot internals except a narrowly reviewed queue contract.
- Agent E: Task 8, export/project/docs/QA. Starts after the app shell is integrated and owns packaging and release files.

The main integrator must review each agent summary, check that no agent edited outside its ownership list, run the focused tests after each merge, then run the full suite and visual QA. If MATLAB is unavailable, stop at source/test/manifest review and report the missing verification rather than claiming the app is complete.


