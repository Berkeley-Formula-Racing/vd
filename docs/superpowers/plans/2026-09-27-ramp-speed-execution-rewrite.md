# Ramp-Speed Execution Rewrite Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the fragile Ramp Speed execution path with a deterministic, diagnosable MATLAB pipeline while preserving the current App Designer layout, setup editor, setup overlays, and MATLAB-only sharing model.

**Architecture:** Keep `RampSpeedApp.mlapp` as a view and move state transitions into a `RampSpeedSession` controller. A single `StudyExecutor` owns serial/parallel execution and cooperative cancellation, a study coordinator owns case scheduling, and separate lateral/longitudinal point solvers return one typed result for every planned speed. Analysis, persistence, and export consume the same canonical result instead of reconstructing solver output in the UI.

**Tech Stack:** MATLAB App Designer `.mlapp`, MATLAB `classdef`/struct APIs, Optimization Toolbox (`fmincon`/bounded scalar search), optional Parallel Computing Toolbox (`parfeval`, `DataQueue`), MATLAB Unit Test, MAT-file persistence, existing `Car`, `Aero`, `AeroMap`, `Tire2`, and `rampSweep` model code.

**Spec:** `docs/superpowers/specs/2026-09-27-ramp-speed-rewrite-architecture.md`

## Global Constraints

- MATLAB remains the implementation platform.
- One Ramp Speed setup owns exactly one solver-ready `Car` and can run in either lateral or pure-longitudinal mode.
- Pure longitudinal mode means `Ay = 0`, zero steer, zero vehicle lateral velocity, zero yaw rate, zero front slip, and symmetric rear slip unless a diagnostic explicitly proves that assumption invalid.
- Every planned speed produces exactly one result row with status `planned`, `running`, `converged`, `near_feasible`, `infeasible`, `solver_failed`, or `cancelled`.
- Continuous metrics that are unavailable are `NaN`; numeric zero is reserved for a measured zero or a documented Boolean value.
- The public request validator rejects nonfinite and nonpositive speeds for both direct API and App callers.
- Serial execution must not require `parfeval`, `DataQueue`, or a parallel pool. Parallel execution has one scheduling owner.
- Cancellation and progress use one cooperative control contract for serial, parallel, and tests.
- Adaptive sampling schedules only new speed IDs and never recursively re-runs an already solved speed without recording a new attempt.
- Result normalization happens once, keyed by stable speed IDs, not approximate floating-point joins.
- Longitudinal capability plots `aLong_max_mps2`; drag force and drag deceleration are separate metrics.
- Raw ramp plots use the actual ramp coordinate, and invalid samples remain visible as gaps or diagnostic markers with reasons.
- Setup specifications contain editable values and provenance only. Derived roll stiffness, resolved paths, solver `Car` objects, and diagnostics are outputs.
- Aero, tire, and camber assets resolve relative to the project root and are fingerprinted for saved-study provenance.
- Existing studies without setup specifications remain loadable as read-only result data.

## Review Focus

- **Gear-transition speeds:** a result near a gear boundary must use an explicit gear branch and preserve all rejected gear attempts; test at 12.5, 15, 17.5, 20, 22.5, and 25 m/s.
- **No Parallel Computing Toolbox:** serial App/API execution must never construct `DataQueue`, call `parfeval`, or open a pool; test with an injected serial executor.
- **Failed and cancelled speeds:** every planned speed must remain visible with a status and reason, while continuous measurements remain `NaN`; test an interior failure, a final failure, and cancellation during a speed.
- **Legacy studies and raw ramp grouping:** an old study must load read-only, preserve speed/point grouping, and plot the true ramp coordinate rather than repeated vehicle speed; test a fixture with at least 11 speeds and 12 points per speed.
- **Asset reproducibility:** changing the current directory or replacing an aero/tire asset at the same path must not silently change a saved setup; test explicit project-root resolution and asset fingerprints.

## File and ownership map

The implementation uses one fresh implementer and one fresh reviewer per task. Tasks with shared interfaces are sequential; read-only codebase reviews and test-result reviews may run in parallel. No two implementation agents edit the same file at the same time.

### New files

- `Full Car Models/+rampSpeed/validateRequest.m` — public request-domain validation.
- `Full Car Models/+rampSpeed/normalizeRequest.m` — canonical request construction.
- `Full Car Models/+rampSpeed/makeSpeedPlan.m` — stable speed-task creation.
- `Full Car Models/+rampSpeed/makeSpeedResult.m` — one-speed result/status construction.
- `Full Car Models/+rampSpeed/normalizeSetupSpec.m` — data-only setup normalization.
- `Full Car Models/+rampSpeed/loadRampAssets.m` — project-root asset resolution and loading.
- `Full Car Models/+rampSpeed/fingerprintFile.m` — deterministic file fingerprints.
- `Full Car Models/+rampSpeed/solveLongitudinalPoint.m` — reduced pure-longitudinal point solver.
- `Full Car Models/+rampSpeed/makePureLongState.m` — full Car-state assembly from reduced variables.
- `Full Car Models/+rampSpeed/candidateGears.m` — explicit gear candidates for a speed.
- `Full Car Models/+rampSpeed/executeSpeedPlan.m` — one-pass speed-task execution.
- `Full Car Models/+rampSpeed/assembleRun.m` — result assembly with one row per task.
- `Full Car Models/+rampSpeed/StudyExecutor.m` — serial/parallel lifecycle and cancellation.
- `Full Car Models/+rampSpeed/RampSpeedSession.m` — UI-independent application state controller.
- `Full Car Models/+rampSpeed/solveLateralPoint.m` — standardized lateral point-result adapter.
- `Full Car Models/tests/test_rampSpeedRequestContract.m` — request and speed-plan tests.
- `Full Car Models/tests/test_rampSpeedAssetResolution.m` — setup/asset reproducibility tests.
- `Full Car Models/tests/test_rampSpeedLongitudinalPointSolver.m` — reduced-solver tests.
- `Full Car Models/tests/test_rampSpeedExecutionContract.m` — statuses, cancellation, and adaptive execution.
- `Full Car Models/tests/test_rampSpeedSession.m` — controller lifecycle tests.
- `Full Car Models/tests/test_rampSpeedAnalysisContract.m` — metrics, axes, gaps, and exports.
- `Full Car Models/tests/test_rampSpeedEndToEndRealCar.m` — small real-car smoke tests.

### Existing files to modify

- `Full Car Models/carComponents/Car.m` — optional fixed-gear evaluation and explicit asset data path while preserving existing callers.
- `Full Car Models/carComponents/Tire2.m` — explicit camber-ratio asset input; remove bare current-directory loading from the Ramp Speed construction path.
- `Full Car Models/carComponents/Camber_Evaluation.m` — keyed or explicit camber-model data instead of an unkeyed process-wide cache.
- `Full Car Models/utilities/buildSingleRampCar.m` — consume resolved assets once and return build provenance.
- `Full Car Models/+rampSpeed/buildCarFromSetup.m` — normalize data-only specs and avoid loading/validating the same map twice.
- `Full Car Models/+rampSpeed/buildSetupCatalog.m` — build each setup once and retain resolved provenance.
- `Full Car Models/+rampSpeed/runLongitudinalRamp.m` — consume the new point solver and speed-plan executor.
- `Full Car Models/+rampSpeed/runLateralRamp.m` — consume the standardized lateral point contract.
- `Full Car Models/+rampSpeed/runCase.m` — delegate to the unified case contract.
- `Full Car Models/+rampSpeed/runStudy.m` — become the synchronous study core and remove duplicate scheduling responsibilities.
- `Full Car Models/+rampSpeed/normalizeRampResult.m` — migrate legacy data once using stable indices and status preservation.
- `Full Car Models/+rampSpeed/makeRun.m` and `makeStudy.m` — versioned canonical result construction.
- `Full Car Models/+rampSpeed/validateStudy.m` — one authoritative validator for runtime, save, load, and export.
- `Full Car Models/+rampSpeed/metricCatalog.m`, `buildPlotData.m`, `renderMetric.m` — explicit x-axis and validity policies.
- `Full Car Models/+rampSpeed/buildComparison.m` — metric-specific truncation/validity policy.
- `Full Car Models/+rampSpeed/buildInspectorTable.m` — selectable speed/point and full detail fields.
- `Full Car Models/+rampSpeed/saveStudy.m`, `loadStudy.m`, `exportStudy.m` — setup/provenance round-trip and full-study export.
- `Full Car Models/+rampSpeed/validateAeroMapFile.m`, `aeroMapCatalog.m` — asset resolver integration and fingerprints.
- `Full Car Models/RampSpeedApp.mlapp` — thin view/controller integration while preserving the layout.
- `Full Car Models/runRampSpeedStudy.m` — compatibility wrapper into the canonical study path.
- `Full Car Models/README_ramp_speed_app.md` and active Ramp Speed design documents — reconcile the single-setup model, units, statuses, and entry points.

---

### Task 1: Freeze request, status, and speed-task contracts

**Agent lane:** contract/data-model implementer. This task owns only new request/result helpers and their tests.

**Files:**
- Create: `Full Car Models/+rampSpeed/validateRequest.m`
- Create: `Full Car Models/+rampSpeed/normalizeRequest.m`
- Create: `Full Car Models/+rampSpeed/makeSpeedPlan.m`
- Create: `Full Car Models/+rampSpeed/makeSpeedResult.m`
- Create: `Full Car Models/tests/test_rampSpeedRequestContract.m`

**Interfaces:**
- Consumes: raw App/API request structs with `rampType`, `settings`, `execution`, `caseIds`, and `solverProfile` fields.
- Produces: `request = rampSpeed.normalizeRequest(rawRequest)` and `plan = rampSpeed.makeSpeedPlan(request)`.
- `request.settings.speeds_mps` is a strictly increasing column vector of finite positive doubles.
- `plan.tasks` is a table with `speedIndex`, `speed_mps`, `origin`, `passIndex`, and `status` columns.
- `makeSpeedResult(task,status,metrics,diagnostics,attempts)` always returns the seven-field status vocabulary from the Global Constraints.

- [ ] **Step 1: Write failing request-domain tests.** Add tests that reject zero, negative, `NaN`, `Inf`, duplicate, and non-monotonic speeds through `validateRequest`, and accept `[5; 10; 20]` for both `lateral` and `longitudinal` requests.

```matlab
testCase.verifyError(@() rampSpeed.validateRequest(struct( ...
    'rampType', "longitudinal", 'settings', struct('speeds_mps', 0))), ...
    'rampSpeed:invalidSpeedDomain');
```

- [ ] **Step 2: Run the focused test and verify the new contract fails.**

Run from `Full Car Models`:

```matlab
results = runtests("tests/test_rampSpeedRequestContract.m");
assertSuccess(results);
```

Expected: FAIL because the new functions do not exist.

- [ ] **Step 3: Implement canonical request normalization.** Normalize string/character ramp types, solver profile, execution mode, case IDs, and speed-grid policy. Reject lateral settings that request an unsupported adaptive mode instead of silently ignoring the control. Do not reference `parallel.pool` from this task.

- [ ] **Step 4: Implement stable speed tasks and explicit result statuses.** Use one numeric `speedIndex` per planned speed and retain `origin` (`requested`, `seed`, `refined`, or `retry`) and `passIndex`. Initialize absent continuous metrics to `NaN`, never zero.

- [ ] **Step 5: Run the focused test and verify it passes.**

```matlab
results = runtests("tests/test_rampSpeedRequestContract.m");
assertSuccess(results);
```

- [ ] **Step 6: Commit the contract slice.**

```text
feat: define ramp speed request and speed result contracts
```

---

### Task 2: Make setup and asset construction deterministic

**Agent lane:** setup/model reproducibility implementer. Run after Task 1; it owns setup construction and asset-resolution files, not solver orchestration.

**Files:**
- Create: `Full Car Models/+rampSpeed/normalizeSetupSpec.m`
- Create: `Full Car Models/+rampSpeed/loadRampAssets.m`
- Create: `Full Car Models/+rampSpeed/fingerprintFile.m`
- Modify: `Full Car Models/carComponents/Tire2.m`
- Modify: `Full Car Models/carComponents/Camber_Evaluation.m`
- Modify: `Full Car Models/utilities/buildSingleRampCar.m`
- Modify: `Full Car Models/+rampSpeed/buildCarFromSetup.m`
- Modify: `Full Car Models/+rampSpeed/buildSetupCatalog.m`
- Create: `Full Car Models/tests/test_rampSpeedAssetResolution.m`

**Interfaces:**
- Consumes: the existing baseline config and editable setup fields.
- Produces: a data-only normalized setup, a resolved asset bundle, one solver-ready `Car`, and build provenance containing asset fingerprints.
- `normalizedSpec` contains identity, editable values, `schemaVersion`, `baselineVersion`, and `aeroMapId`; it does not contain `derived`, absolute paths, or a `Car`.
- `assets = rampSpeed.loadRampAssets(config, setup)` resolves all paths from the baseline file/project root and loads each asset once.
- `buildCarFromSetup` preserves its existing primary `Car` output and adds a provenance output without breaking current callers.

- [ ] **Step 1: Add failing tests for current-directory independence and one-time resolution.** Change the test working directory before building the baseline, verify the build succeeds, verify the selected map fingerprint is stored, and verify duplicate setup construction does not produce stale derived values.

- [ ] **Step 2: Run the focused setup tests and verify the new cases fail.**

```matlab
results = runtests(["tests/test_rampSpeedAssetResolution.m", ...
                    "tests/test_rampSpeedSetupModel.m"]);
assertSuccess(results);
```

Expected: FAIL on missing resolver/fingerprint behavior.

- [ ] **Step 3: Implement data-only setup normalization.** Strip or recompute derived fields, validate `schemaVersion`, validate discrete spring/ARB selections, preserve setup IDs, and make baseline identity derived from the protected baseline ID rather than a copied flag.

- [ ] **Step 4: Implement explicit asset loading and fingerprints.** Resolve tire, camber, and aero files relative to the baseline/config root. Use a SHA-256 fingerprint over file bytes. Pass camber-ratio/model data explicitly into the Ramp Speed construction path and key any in-memory cache by the asset fingerprint.

- [ ] **Step 5: Make car construction load the aero map once.** Validate the already-resolved map object, pass it into `buildSingleRampCar`, and retain the map ID/path/fingerprint in build provenance. Build exactly one `Car` per setup catalog entry.

- [ ] **Step 6: Run the focused setup and baseline tests.**

```matlab
results = runtests(["tests/test_rampSpeedAssetResolution.m", ...
                    "tests/test_rampSpeedSetupModel.m", ...
                    "tests/test_rampSpeedDriverWeight.m"]);
assertSuccess(results);
```

- [ ] **Step 7: Commit the setup slice.**

```text
refactor: make ramp speed setup construction deterministic
```

---

### Task 3: Implement the reduced pure-longitudinal point solver

**Agent lane:** numerical solver implementer. Run after Task 2; this task owns the new solver and the fixed-gear evaluation seam.

**Files:**
- Create: `Full Car Models/+rampSpeed/candidateGears.m`
- Create: `Full Car Models/+rampSpeed/makePureLongState.m`
- Create: `Full Car Models/+rampSpeed/solveLongitudinalPoint.m`
- Modify: `Full Car Models/carComponents/Car.m`
- Create: `Full Car Models/tests/test_rampSpeedLongitudinalPointSolver.m`

**Interfaces:**
- Consumes: one setup `Car`, one positive `speed_mps`, optional prior point seed, solver profile, and `control.shouldCancel`.
- Produces: one `SpeedResult` with `metrics`, `diagnostics`, explicit gear attempts, and `status`.
- `rampSpeed.solveLongitudinalPoint(car,speed_mps,seed,profile,control)` is the only new Ramp Speed entry point for one pure-longitudinal speed.
- `Car.equations` and `Car.constraint1` gain an optional evaluation-options struct with `gearOverride`; all existing call signatures remain valid.

- [ ] **Step 1: Write failing solver-contract tests.** Test that the result state has zero steer, zero lateral velocity, zero yaw rate, zero front slip, equal rear slip, and `abs(gLat) <= 1e-12`. Test explicit gear attempts at 12.5–25 m/s and cancellation before the first evaluation.

```matlab
result = rampSpeed.solveLongitudinalPoint(car,20,[],profile,control);
testCase.verifyTrue(any([result.diagnostics.gearAttempts.gear] == 3));
testCase.verifyLessThanOrEqual(abs(result.metrics.gLat),1e-12);
testCase.verifyEqual(result.state(1),0);
testCase.verifyEqual(result.state(4:5),[0 0]);
testCase.verifyEqual(result.state(6:7),[0 0]);
testCase.verifyEqual(result.state(8),result.state(9),'AbsTol',1e-12);
```

- [ ] **Step 2: Run the solver tests and verify they fail before implementation.**

```matlab
results = runtests("tests/test_rampSpeedLongitudinalPointSolver.m");
assertSuccess(results);
```

- [ ] **Step 3: Add fixed-gear evaluation without changing legacy callers.** Apply `gearOverride` only when supplied; otherwise preserve the existing powertrain gear-selection behavior. Ensure the override is included in diagnostics and does not mutate the source `Car`.

- [ ] **Step 4: Build the pure-longitudinal state from reduced variables.** Use the full state vector only as a model interface, with independent variables `[throttle, rearSlip]`. Set steer, lateral velocity, yaw rate, vehicle speed, and front slips directly. Set both rear slips equal and record the reduced state in diagnostics.

- [ ] **Step 5: Solve each candidate gear explicitly.** For each candidate gear, solve rear wheel equilibrium over the bounded rear-slip interval, then maximize longitudinal acceleration over bounded throttle. Reject candidates for RPM, slip-angle, wheel-load, aero, or residual violations. Store every candidate’s exit flag, residuals, and rejection reason.

- [ ] **Step 6: Retain a reference fallback.** If the reduced solver cannot produce an acceptable candidate, call the existing `max_long_accel` once as a diagnostic fallback, mark the attempt `reference_fallback`, and never convert an unsuccessful fallback into a valid result.

- [ ] **Step 7: Compare against the existing solver at representative speeds.** The reduced solution must agree with the reference solver within an explicitly tested tolerance at 5, 10, 15, 20, and 25 m/s when both are feasible. The test must compare acceleration, gear, rear slip, and residual status, not only the plotted value.

- [ ] **Step 8: Run focused solver and model tests.**

```matlab
results = runtests(["tests/test_rampSpeedLongitudinalPointSolver.m", ...
                    "tests/test_longitudinalRampAdapter.m", ...
                    "tests/test_rampSpeedSolverProfiles.m"]);
assertSuccess(results);
```

- [ ] **Step 9: Commit the solver slice.**

```text
feat: add explicit-gear reduced longitudinal ramp solver
```

---

### Task 4: Replace recursive adaptive execution with one speed-plan executor

**Agent lane:** speed-planning implementer. Run after Tasks 1 and 3; it owns longitudinal run orchestration and adaptive state, not the App.

**Files:**
- Create: `Full Car Models/+rampSpeed/executeSpeedPlan.m`
- Create: `Full Car Models/+rampSpeed/assembleRun.m`
- Modify: `Full Car Models/+rampSpeed/runLongitudinalRamp.m`
- Modify: `Full Car Models/+rampSpeed/planAdaptiveSpeeds.m`
- Modify: `Full Car Models/+rampSpeed/refineAdaptiveSpeeds.m`
- Create: `Full Car Models/tests/test_rampSpeedExecutionContract.m`

**Interfaces:**
- Consumes: `request`, `SpeedTask` table, a point-solver function handle, and the cooperative control struct.
- Produces: one `CaseRun` with one speed row per planned task, adaptive history, attempt history, and no duplicate speed evaluation.
- `runLongitudinalRamp` delegates each speed to `solveLongitudinalPoint` and no longer calls itself recursively for adaptive batches.

- [ ] **Step 1: Write failing execution tests.** Test one row per planned speed, distinct `solver_failed` versus `infeasible`, preserved cancellation, no duplicate speed IDs across adaptive passes, and `NaN` for absent continuous metrics.

- [ ] **Step 2: Run the focused execution tests and verify they fail.**

```matlab
results = runtests("tests/test_rampSpeedExecutionContract.m");
assertSuccess(results);
```

- [ ] **Step 3: Implement `executeSpeedPlan`.** Iterate over tasks once, check `control.shouldCancel()` before each task, emit structured progress after each task, call the injected point solver, and store the returned result by `speedIndex`.

- [ ] **Step 4: Implement `assembleRun`.** Preallocate all planned rows with explicit status and `NaN` continuous fields. Replace a row only with the corresponding `SpeedResult`; never infer a missing row from neighboring speeds.

- [ ] **Step 5: Rework adaptive refinement.** Keep a cache keyed by setup fingerprint, ramp type, solver profile, and `speedIndex`. The planner may request a new task only if that key is absent. Store `passIndex`, `origin`, and `refinementReason` for every task.

- [ ] **Step 6: Preserve near-feasible and invalid diagnostics.** Separate physical infeasibility, optimizer failure, cancellation, and missing results in both per-speed status and `speedErrors`.

- [ ] **Step 7: Run adapter and execution tests.**

```matlab
results = runtests(["tests/test_rampSpeedExecutionContract.m", ...
                    "tests/test_longitudinalRampAdapter.m", ...
                    "tests/test_rampSpeedAdaptiveSpeedGrid.m"]);
assertSuccess(results);
```

- [ ] **Step 8: Commit the speed execution slice.**

```text
refactor: execute adaptive ramp speeds through one task plan
```

---

### Task 5: Standardize the lateral solver adapter

**Agent lane:** lateral solver implementer. Run after Task 1; it may retain `rampSweep` internally but must expose the same speed-result contract as Task 4.

**Files:**
- Create: `Full Car Models/+rampSpeed/solveLateralPoint.m`
- Modify: `Full Car Models/+rampSpeed/runLateralRamp.m`
- Modify: `Full Car Models/sweeps/rampSweep.m`
- Create/modify: `Full Car Models/tests/test_rampSpeedLateralExecutionContract.m`

**Interfaces:**
- Consumes: one `Car`, a positive speed task, ramp policy, profile, and cooperative control.
- Produces: one `SpeedResult` with free-limit, ramp-point, balance, tire, camber, and slip-angle metrics.
- Legacy `rampSweep` remains callable, but its sparse result is adapted once at the boundary and does not define the canonical study schema.

- [ ] **Step 1: Add failing tests for a failed interior speed and a successful speed.** The returned run must retain both rows; the failed row must have a reason and no interpolated continuous metrics.

- [ ] **Step 2: Run the lateral contract test and verify the missing-result behavior fails.**

```matlab
results = runtests("tests/test_rampSpeedLateralExecutionContract.m");
assertSuccess(results);
```

- [ ] **Step 3: Wrap the existing lateral limit/ramp solve in `solveLateralPoint`.** Convert solver exceptions, omitted speeds, failed ramp points, and cancellation into explicit statuses. Preserve the detailed raw point payload for analysis.

- [ ] **Step 4: Apply the same profile settings and cancellation checks to lateral and longitudinal branches.** Remove hard-coded lateral optimizer budgets where they conflict with the resolved profile.

- [ ] **Step 5: Rebuild requested-speed rows from stable task IDs.** Do not match rows by approximate speed values and do not bridge an invalid interior speed with a line.

- [ ] **Step 6: Run lateral regression tests.**

```matlab
results = runtests(["tests/test_rampSpeedLateralExecutionContract.m", ...
                    "tests/test_rampSpeedLateralAdapter.m", ...
                    "tests/test_rampSpeedRampSweep.m"]);
assertSuccess(results);
```

- [ ] **Step 7: Commit the lateral slice.**

```text
refactor: standardize lateral ramp speed results
```

---

### Task 6: Split study coordination from execution lifecycle

**Agent lane:** integration/concurrency implementer. Run after Tasks 4 and 5; this task owns the controller-facing run lifecycle, not App callbacks.

**Files:**
- Create: `Full Car Models/+rampSpeed/StudyExecutor.m`
- Modify: `Full Car Models/+rampSpeed/runStudy.m`
- Modify: `Full Car Models/+rampSpeed/runCase.m`
- Create: `Full Car Models/+rampSpeed/RampSpeedSession.m`
- Create: `Full Car Models/tests/test_rampSpeedSession.m`
- Modify: `Full Car Models/tests/test_rampSpeedRunner.m`

**Interfaces:**
- `StudyExecutor.start(cars,cases,request,callbacks)` returns a job handle/state object.
- `StudyExecutor.cancel(job)` signals the same cooperative token read by the study core.
- `rampSpeed.runStudy` is synchronous and deterministic when called directly; it never creates an outer future.
- `RampSpeedSession` exposes `viewModel`, `editSetup`, `duplicateSetup`, `deleteSetup`, `selectCases`, `start`, `cancel`, `acceptResult`, `save`, `load`, `export`, and `clearResults`.

- [ ] **Step 1: Add failing controller tests.** Test setup selection, duplicate/delete protection, state transitions `idle → running → completed`, `idle → running → cancelled`, and rejection of setup edits while running.

- [ ] **Step 2: Add a serial executor test with parallel APIs unavailable.** Inject a point solver and verify no `parallel.pool.DataQueue`, `parfeval`, or pool creation is attempted in serial mode.

- [ ] **Step 3: Run the new tests and verify they fail.**

```matlab
results = runtests(["tests/test_rampSpeedSession.m", ...
                    "tests/test_rampSpeedRunner.m"]);
assertSuccess(results);
```

- [ ] **Step 4: Move request/case normalization and final assembly into the synchronous study core.** `runStudy` may report progress and checkpoints, but it must not own UI futures or create an additional scheduler.

- [ ] **Step 5: Implement `StudyExecutor`.** In serial mode, call the synchronous core directly. In parallel mode, create exactly one outer future and pass `parallelRequested=false` into the worker request so nested scheduling cannot occur. Use a unique checkpoint path per job.

- [ ] **Step 6: Connect cancellation and progress.** `cancel(job)` flips the shared token; the core checks it before each case, speed, and adaptive pass and writes terminal status before returning.

- [ ] **Step 7: Implement `RampSpeedSession` around setup catalog and executor.** Keep draft setup edits separate from committed setup specs; rebuild only on Apply/Duplicate/Delete, not on every control keystroke.

- [ ] **Step 8: Run runner, session, and cancellation tests.**

```matlab
results = runtests(["tests/test_rampSpeedSession.m", ...
                    "tests/test_rampSpeedRunner.m", ...
                    "tests/test_rampSpeedCancellation.m", ...
                    "tests/test_rampSpeedProgress.m"]);
assertSuccess(results);
```

- [ ] **Step 9: Commit the lifecycle slice.**

```text
refactor: separate ramp study coordination from execution lifecycle
```

---

### Task 7: Rebuild canonical normalization, metrics, comparison, and persistence

**Agent lane:** result/analysis implementer. Run after Tasks 1, 4, and 5; it owns result consumers and backward compatibility.

**Files:**
- Modify: `Full Car Models/+rampSpeed/normalizeRampResult.m`
- Modify: `Full Car Models/+rampSpeed/makeRun.m`
- Modify: `Full Car Models/+rampSpeed/makeStudy.m`
- Modify: `Full Car Models/+rampSpeed/validateStudy.m`
- Modify: `Full Car Models/+rampSpeed/metricCatalog.m`
- Modify: `Full Car Models/+rampSpeed/buildPlotData.m`
- Modify: `Full Car Models/+rampSpeed/renderMetric.m`
- Modify: `Full Car Models/+rampSpeed/buildComparison.m`
- Modify: `Full Car Models/+rampSpeed/buildInspectorTable.m`
- Modify: `Full Car Models/+rampSpeed/saveStudy.m`
- Modify: `Full Car Models/+rampSpeed/loadStudy.m`
- Modify: `Full Car Models/+rampSpeed/exportStudy.m`
- Create: `Full Car Models/tests/test_rampSpeedAnalysisContract.m`

**Interfaces:**
- Consumes: canonical `Study`, `CaseRun`, and `SpeedResult` data from Tasks 4–6 plus schema-v1 legacy studies.
- Produces: schema-v2 canonical data with migration metadata, explicit metric definitions, selectable inspector data, and full-study exports.
- `metricCatalog` entries define `source`, `xSource`, `units`, `validityRule`, `truncationPolicy`, and supported `rampTypes`.

- [ ] **Step 1: Add failing analysis tests.** Cover longitudinal capability from `aLong_max_mps2`, separate drag force/deceleration, raw ramp x-coordinate, invalid-row gaps/markers, camber FL/FR and RL/RR, slip angles, and one valid speed with a physically zero metric.

```matlab
catalog = rampSpeed.metricCatalog("longitudinal");
capability = catalog(strcmp(string({catalog.id}),"capability_longitudinal"));
testCase.verifyEqual(capability.source,"aLong_max_mps2");
```

- [ ] **Step 2: Add a legacy grouping fixture.** Use 11 speeds and 12 points per speed with no explicit legacy indices. Assert that migration reconstructs 11 groups of 12 points from source speed values and preserves within-speed order.

- [ ] **Step 3: Run the analysis tests and verify they fail.**

```matlab
results = runtests("tests/test_rampSpeedAnalysisContract.m");
assertSuccess(results);
```

- [ ] **Step 4: Make migration and validation status-preserving.** Keep `planned`, `failed`, `solver_failed`, `infeasible`, and `cancelled` distinct. Load studies without setup specifications as read-only. Validate setup specifications, asset fingerprints, result references, and required typed columns on save and load.

- [ ] **Step 5: Replace floating-point row matching with stable speed IDs.** Preserve legacy indices when present; otherwise derive indices from grouped source-speed records and record the migration rule in metadata.

- [ ] **Step 6: Expand the metric catalog and plot data.** Add longitudinal capability, force/deceleration distinctions, per-corner camber, per-corner slip angle, mechanical balance, aero balance, yaw-derived understeer/oversteer, lift, drag, and ride-height metrics. Make each metric choose its own x-axis and validity/truncation policy.

- [ ] **Step 7: Make comparison metric-specific.** A truncated sustainable-limit metric may be invalid while aero/load metrics remain comparable. Never blanket-mask an entire comparison because one unrelated metric is truncated.

- [ ] **Step 8: Make inspector and export complete.** Accept selected run, speed index, and point index. Export the full canonical study plus setup/provenance metadata through `exportStudy`, not only the first inspector table.

- [ ] **Step 9: Run persistence and analysis regression tests.**

```matlab
results = runtests(["tests/test_rampSpeedAnalysisContract.m", ...
                    "tests/test_rampSpeedSchema.m", ...
                    "tests/test_rampSpeedPersistence.m", ...
                    "tests/test_rampSpeedComparison.m", ...
                    "tests/test_rampSpeedPlotCatalog.m"]);
assertSuccess(results);
```

- [ ] **Step 10: Commit the result/analysis slice.**

```text
refactor: make ramp speed results and metrics status-aware
```

---

### Task 8: Reconnect the App Designer UI to the session/controller

**Agent lane:** App Designer implementer. Run after Task 6 and Task 7; this task owns the `.mlapp` package and UI smoke tests only.

**Files:**
- Modify: `Full Car Models/RampSpeedApp.mlapp`
- Modify: `Full Car Models/tests/test_rampSpeedAppRequest.m`
- Modify: `Full Car Models/tests/test_rampSpeedAppSmoke.m`
- Modify: `Full Car Models/tests/test_rampSpeedSetupEditorApp.m`
- Create/modify: `Full Car Models/tests/test_rampSpeedAppLifecycle.m`

**Interfaces:**
- The App owns controls and rendering; `RampSpeedSession` owns setup/run state.
- Run callbacks call `session.start`; Cancel calls `session.cancel`; completion calls `session.acceptResult`.
- Serial mode never constructs a DataQueue or Future. Parallel mode uses the executor’s one job handle.
- The existing left setup panel and right analysis tabs remain in place.

- [ ] **Step 1: Add failing App lifecycle tests with an injected session/executor.** Cover serial Run, parallel Run request construction, Cancel, progress, Clear, Save, Load, setup-edit lockout, and read-only legacy study behavior.

- [ ] **Step 2: Run App smoke tests and verify lifecycle failures.**

```matlab
results = runtests(["tests/test_rampSpeedAppRequest.m", ...
                    "tests/test_rampSpeedAppSmoke.m", ...
                    "tests/test_rampSpeedSetupEditorApp.m", ...
                    "tests/test_rampSpeedAppLifecycle.m"]);
assertSuccess(results);
```

- [ ] **Step 3: Replace App-owned run state with the session view model.** Remove direct solver scheduling from `RunButtonPushed`, `CancelButtonPushed`, and completion callbacks. Keep the App’s graphics handles and control references local to the view.

- [ ] **Step 4: Add analysis selectors.** Add metric, setup overlay, speed, and ramp-point selectors. Disable incompatible controls by ramp type instead of silently ignoring them.

- [ ] **Step 5: Add truthful progress and failure display.** Show setup ID, speed, status, reason, completed/total counts, and cancellation state. Do not render invalid numeric samples as zero.

- [ ] **Step 6: Fix clear, resize, and export integration.** Clear every graphics layer through the session analysis state, recompute layout on any window size, and route export through `session.export`/`rampSpeed.exportStudy`.

- [ ] **Step 7: Exercise the packaged `.mlapp`.** Open the saved App Designer package in MATLAB, run a serial fixture study, run a real baseline smoke study, switch lateral/longitudinal mode, duplicate a setup, and verify the persistent layout remains usable at small and large window sizes.

- [ ] **Step 8: Commit the App slice.**

```text
refactor: make Ramp Speed App a thin session-driven view
```

---

### Task 9: Retire the split entry point and update documentation

**Agent lane:** compatibility/documentation implementer. Run after Tasks 6–8.

**Files:**
- Modify: `Full Car Models/runRampSpeedStudy.m`
- Modify: `Full Car Models/README_ramp_speed_app.md`
- Modify: `Full Car Models/docs/superpowers/specs/2026-09-19-ramp-speed-app-design.md`
- Modify: `Full Car Models/docs/superpowers/plans/2026-09-19-ramp-speed-app-implementation.md`
- Modify: `Full Car Models/docs/superpowers/plans/2026-09-22-ramp-speed-baseline-setup-editor.md`
- Create: `Full Car Models/tests/test_rampSpeedEndToEndRealCar.m`

**Interfaces:**
- `runRampSpeedStudy` becomes a compatibility wrapper that builds the baseline setup catalog and calls the canonical study path.
- Documentation names `RampSpeedApp.mlapp` and the canonical `rampSpeed` API as the supported entry points.

- [ ] **Step 1: Add the real-car end-to-end smoke test.** Run one baseline setup through both ramp types at `[5 10 15 17.5 20 22.5 25]` m/s in serial mode. Assert one result row per speed, explicit status for every row, no invalid row with finite continuous zero substituted for missing data, and no exception from setup construction.

- [ ] **Step 2: Run the end-to-end test before changing the wrapper.**

```matlab
results = runtests("tests/test_rampSpeedEndToEndRealCar.m");
assertSuccess(results);
```

- [ ] **Step 3: Route the legacy script through the canonical baseline/setup catalog.** Preserve its command-line options and returned variables through a documented compatibility adapter, but remove its direct legacy plotting/data path.

- [ ] **Step 4: Reconcile all active documentation.** State the one-Car-per-setup model, pure-longitudinal definition, setup-input units versus SI result units, status vocabulary, serial/parallel behavior, map IDs, and supported save/load compatibility.

- [ ] **Step 5: Add the complete real-car and performance checks.** Measure map construction count, repeated-speed cache hits, per-speed solve time, total study time, and memory for a 7-speed baseline run. Record the measurements in the test output without asserting a machine-specific absolute runtime.

- [ ] **Step 6: Run the broad Ramp Speed regression suite.**

```matlab
results = runtests("tests", "IncludeSubfolders", true);
assertSuccess(results);
```

- [ ] **Step 7: Run repository hygiene checks.**

```text
git diff --check
git status --short
```

- [ ] **Step 8: Commit the compatibility/documentation slice.**

```text
docs: document canonical ramp speed execution path
```

---

## Sub-agent execution plan

Use `superpowers:subagent-driven-development` during implementation. Every task receives a fresh implementer and a separate reviewer. Use `gpt-6-luna` with `max` reasoning for the mechanical and isolated tasks; use the strongest available reviewer for Tasks 3, 6, 7, and the final whole-branch review.

Implementation order is sequential at the shared-interface boundaries:

```text
Task 1 ─┬─ Task 2 ─ Task 3 ─ Task 4 ─┬─ Task 6 ─ Task 8 ─ Task 9
        └─ Task 5 ───────────────────┘
                         └─ Task 7 ──┘
```

Read-only agents may inspect Tasks 2–5, test coverage, and MATLAB packaging in parallel before each implementation dispatch. Implementation agents do not spawn their own agents and do not edit files outside their task ownership. Every completed task gets a spec-compliance and code-quality review before the next dependent task starts.

## Verification commands

Focused slices:

```matlab
results = runtests([ ...
    "tests/test_rampSpeedRequestContract.m", ...
    "tests/test_rampSpeedAssetResolution.m", ...
    "tests/test_rampSpeedLongitudinalPointSolver.m", ...
    "tests/test_rampSpeedExecutionContract.m", ...
    "tests/test_rampSpeedLateralExecutionContract.m", ...
    "tests/test_rampSpeedSession.m", ...
    "tests/test_rampSpeedAnalysisContract.m", ...
    "tests/test_rampSpeedEndToEndRealCar.m"]);
assertSuccess(results);
```

App and regression checks:

```matlab
results = runtests([ ...
    "tests/test_rampSpeedAppRequest.m", ...
    "tests/test_rampSpeedAppSmoke.m", ...
    "tests/test_rampSpeedAppLifecycle.m", ...
    "tests/test_rampSpeedPersistence.m", ...
    "tests/test_rampSpeedComparison.m", ...
    "tests/test_rampSpeedPlotCatalog.m"]);
assertSuccess(results);
```

The final whole-branch review must inspect the implementation diff against `docs/superpowers/specs/2026-09-27-ramp-speed-rewrite-architecture.md`, review all deferred findings recorded in the SDD ledger, rerun the real-car smoke, and confirm that no second active Ramp Speed execution path remains.

## Self-review

- Spec coverage: request validation is Task 1; deterministic setup/assets Task 2; explicit pure-longitudinal physics Task 3; adaptive execution Task 4; lateral compatibility Task 5; lifecycle/cancellation Task 6; analysis/persistence Task 7; UI behavior Task 8; legacy/docs/performance Task 9.
- Placeholder scan: no task depends on a `TODO`, `TBD`, or unspecified edge-case step; every task names files, interfaces, tests, commands, and commit scope.
- Interface consistency: Tasks 1–5 produce `SpeedTask`/`SpeedResult` contracts; Tasks 4–6 consume them; Task 7 consumes the assembled canonical run; Tasks 8–9 consume the session and persistence APIs.
- Review focus coverage: gear transitions Task 3; no-PCT serial mode Task 6; failed/cancelled speeds Task 4 and Task 6; legacy raw grouping Task 7; asset reproducibility Task 2.
