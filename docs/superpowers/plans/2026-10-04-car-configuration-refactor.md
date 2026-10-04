# Car Configuration Refactor Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking. Recommended execution: native implementation in one coordinated checkout, with a final independent review if authorized. Do not begin implementation merely because this planning document exists.

**Goal:** Make individual vehicle configurations independently editable and consistently buildable, add named car cases, and preserve every existing car-generation entrypoint and saved-result contract.

**Architecture:** Keep the existing `Car` value class and low-level constructor. Introduce `+vehicleConfig` for setup defaults, a shared parameter registry, resolved sampling, asset loading, single-car construction and named case construction. Existing `carConfig`, `parameters_loop`, ramp builders and launchers stay compatible through adapters; named output is opt-in.

**Tech Stack:** MATLAB R2026a on Windows, existing vehicle/tyre/aero classes, MATLAB function-based tests and assertion scripts, existing App Designer archive and MAT/HDF5 persistence. No new dependencies or toolbox requirements.

**Spec:** [Car configuration refactor design](../specs/2026-10-04-car-configuration-refactor-design.md). Both documents are drafts for the user's planning review. Implementation begins only after that review and an execution instruction.

## Global Constraints

- Keep `Car`, `Aero`, `Powertrain`, `Tire2` and `NondimensionalTire` as value classes; keep the positional `Car` constructor compatible.
- Preserve every existing `carConfig` call with zero to four inputs and one to four outputs; default output remains an N-by-2 cell array with lap then acceleration columns.
- Preserve the ordered 57-column QSS design-table schema, case ordering, baseline selection, DOE random-seed behavior and existing public warning/error identifiers.
- Preserve separate QSS and frozen ramp baselines, their aero-map role policies, ride-height datums, tyre model selection and roll-stiffness policies.
- Preserve `carCell` MAT/checkpoint payloads, ramp study schemas and telemetry HDF5/CSV schemas; do not convert existing files in place.
- Configuration edits followed by rebuilding must produce consistent derived fields and a fresh unsolved car; the supported edit path must not reuse stale G-G data or event-result handles.
- Load immutable map/camber assets once per batch or distinct asset identity; do not add file reads or configuration validation to solver evaluation loops.
- Do not run simulations during planning; implementation validation must use bounded runs and record which launchers were actually executed.
- Preserve unrelated dirty files and existing user launchers; no cleanup, baseline retuning or solver-physics changes are included.

## Review Focus

1. Empty, singleton and multirow Explicit designs must preserve table column order, row-to-car correspondence and 0-by-2/1-by-2/N-by-2 legacy shapes (Tasks 1 and 4).
2. Editing wheelbase, wheel radius, final drive or springs must update every associated derived value when rebuilding, without mutating the original setup or reusing results (Task 3).
3. QSS acceleration and frozen ramp cars use different map/datum/roll policies; consolidation must preserve both, including LC0 selection and uncertainty (Tasks 3 and 5).
4. Historic `carCell` caches/checkpoints containing `Car` and `Events2` objects must remain readable; packing named cases must preserve attached results rather than rebuild them (Tasks 4 and 6).
5. Path shadowing, foreign working directories, package release contents and two same-named ramp entrypoints must not select a wrong class or script (Tasks 6 and 7).

## Starting state and planning evidence

Source inspection on 2026-10-04 found branch `master`, HEAD `edee80a`, one commit ahead of `origin/master`. There are existing modified ramp UI/session/test files and untracked launchers/tests/`OLD`. This document does not freeze or authorize altering those changes. Reinspect at execution time; do not assume a clean base or that an isolated worktree contains uncommitted inputs.

The repository ignores `docs/superpowers`; these two planning files are saved locally and are not committed. Before moving execution into another checkout, copy both the plan and linked spec into that workspace or deliberately track them as part of the approved workflow. Do not assume a worktree created from HEAD contains the plan or current dirty inputs.

No MATLAB tests or simulations were run while writing this plan. MATLAB is discoverable at `C:\Program Files\MATLAB\R2026a\bin\matlab.exe`; toolbox availability and runtime parity remain to be established during implementation.

Three baseline discrepancies are already visible in source: the SteadyState tyre/test disagreement and deleted root launcher expectation; two different `runRampSpeedStudy` entrypoints; and an older vibration script assigning undeclared `car.c`. Record their actual baseline results before changing expectations. `RH.m` is untracked and assumes paths are already configured.

## File responsibilities and public interfaces

All paths below are relative to `C:\VD`; if execution uses a managed worktree, substitute its workspace root in commands.

| File | Responsibility |
| --- | --- |
| `Full Car Models/+vehicleConfig/parameterNames.m` | Exact ordered QSS parameter-name list; replaces DOE source parsing |
| `Full Car Models/+vehicleConfig/defaultSetup.m` | Current scalar QSS numeric/fixed/map defaults and event defaults; constructs no car |
| `Full Car Models/+vehicleConfig/fromLegacyGroups.m` | Translate six legacy input groups, including vector sampling templates |
| `Full Car Models/+vehicleConfig/sampleSetups.m` | Sample/resolve numeric rows and produce independent setup structs plus the unchanged design table |
| `Full Car Models/+vehicleConfig/loadAssets.m` | Absolute-path map and camber loading for an explicit build context |
| `Full Car Models/+vehicleConfig/buildCar.m` | Validate a resolved setup; construct components and one fresh car using the low-level constructor |
| `Full Car Models/+vehicleConfig/buildCases.m` | Produce N-by-1 named case structs from resolved QSS setups |
| `Full Car Models/+vehicleConfig/toLegacyCells.m` | Pack named cars into N-by-2 cells without recalculation |
| `Full Car Models/+vehicleConfig/fromRampSetup.m` | Adapt frozen ramp source values and policy without reading QSS defaults |
| `Full Car Models/carConfig.m` | Compatible dispatch/defaults/legacy DOE behavior; optional fifth `"cases"` output format |
| `Full Car Models/utilities/parameters_loop.m` | Compatible six-group sampling/build adapter; preserve signature and public identifiers |
| `Full Car Models/utilities/buildSingleRampCar.m` | Compatible ramp builder and existing derived outputs |
| `Full Car Models/utilities/doeResolveStudy.m` | Consume the shared registry instead of reading MATLAB source text |
| `Full Car Models/entrypoints/car_setup_example.m` | Small runnable construction-only example of editing/building and named cases |

Keep `carConfigBaseline.m`'s independently owned configuration and `utilities/buildCarFromParameterSet.m`'s constructor behavior. They may need comments/import integration but not new baseline constants or constructor argument order.

The exact new signatures are:

```matlab
names = vehicleConfig.parameterNames();                       % 1-by-57 cellstr
[setup,eventParams] = vehicleConfig.defaultSetup(tireModel,uncertaintyMode);
base = vehicleConfig.fromLegacyGroups(cP,aP,eP,DTp,Bp,tP);      % scalar struct
[setups,X] = vehicleConfig.sampleSetups(base,samplingType,samples);
assets = vehicleConfig.loadAssets(setup);                     % scalar struct
[car,derived] = vehicleConfig.buildCar(setup,role,assets);     % scalar Car/struct
cases = vehicleConfig.buildCases(setups,assets);               % N-by-1 struct
carCell = vehicleConfig.toLegacyCells(cases);                  % N-by-2 cell
spec = vehicleConfig.fromRampSetup(fixed,setup,mapInput,assets);
[cars,eventParams,X,B] = carConfig(type,samples,model,uncertainty,outputFormat);
```

Defaults and struct fields are defined in the spec. `buildCar` defaults to `"lap"` and an empty asset struct; `buildCases` loads once if context omitted. New errors use the `vehicleConfig:` prefix; legacy wrappers retain `carConfig:badRequest`, `carConfig:legacyDOE`, `parameters_loop:badExplicitInput`, `parameters_loop:nonScalarBaseline`, `parameters_loop:unknownExplicitParameter`, `parameters_loop:aeroMapNotFound` and `parameters_loop:badSamplingType`. Validation of format uses `carConfig:badOutputFormat`.

## Task 1: Establish baseline evidence and replace the brittle parameter registry

**Files:**

- Create: `Full Car Models/+vehicleConfig/parameterNames.m`.
- Modify: `Full Car Models/utilities/doeResolveStudy.m` (local `parameterRegistry`).
- Create: `Full Car Models/tests/test_vehicleConfigRegistry.m`.
- Create: `Full Car Models/tests/helpers/captureCarConfigurationReference.m`.
- Create: `Full Car Models/tests/helpers/carConfigurationSnapshot.m`.
- Create: `Full Car Models/tests/car_configuration_entrypoints.csv` (validation manifest, not a vehicle dataset).
- Read: `parameters_loop.m`, `DOEStudyConfig.m`, `test_doe_study_config.m`, `test_doe_explicit_cars.m`, `test_entrypointBootstrap.m`, `test_steadyStateLapsim.m` and the matrix below.

**Interfaces:** Consumes current zero/four-output `carConfig`, `carConfigBaseline` and legacy sampling. Produces a 1-by-57 `cellstr` registry and before-change construction evidence; production construction remains unchanged in this task.

- [ ] Recheck status and attached worktrees before choosing a checkout. Preserve all dirty inputs; record root, commit, MATLAB version, resolved class paths and toolbox availability in the reference artifact. Do not reset or clean anything.
- [ ] Before altering production code, capture baseline car properties for both roles, table names/order, event defaults, a two-row Explicit design, frozen ramp setup/config/derived fields, and legacy/LC0 low/nominal/high selections where their assets are available. Store artifacts under an execution-specific ignored directory such as `artifacts/car_configuration_refactor/2026-10-04-reference`; do not overwrite existing files.
- [ ] If bounded numerical validation is authorized for implementation, capture the Task 8 G-G/event/ramp/telemetry references now, before construction routing changes. Record the original functions, dirty-state snapshot, inputs and exact commands. If those runs are deferred, retain a coordinated reproducible before-change snapshot and mark numerical comparison pending; do not manufacture a reference from the refactored code.
- [ ] The snapshot helper records every public property, recursively including component fields; represent function handles by `func2str` plus deterministic torque samples, and exclude only runtime timing/graphics handles. Preserve numeric classes, dimensions and NaN/Inf masks. Save public property lists so added/removed fields cannot be silently ignored. For solved references, record event times/arrays separately and preserve handle identity through native MAT save/load checks.
- [ ] Populate the CSV with every matrix row: `path,kind,constructionRoute,baselineStatus,postStatus,evidencePath,blocker`. Initial statuses are `unverified`, never `pass` inferred from source. Record precise missing assets/path conflicts and pre-existing failures.
- [ ] Write registry tests before adding the registry function:

```matlab
function testRegistryMatchesCurrentDesign(testCase)
    [~,~,X] = carConfig();
    names = vehicleConfig.parameterNames();
    verifyEqual(testCase,names,X.Properties.VariableNames);
    verifyEqual(testCase,numel(names),57);
    verifyEqual(testCase,numel(unique(names)),57);
end

function testDoeUsesEveryRegistryName(testCase)
    [~,~,~,B] = carConfig();
    study = DOEStudyConfig();
    resolved = doeResolveStudy(study,B);
    verifyTrue(testCase,all(ismember(cellstr(resolved.parameters.name), ...
        vehicleConfig.parameterNames())));
end
```

- [ ] Run `runtests('tests/test_vehicleConfigRegistry.m')` after model bootstrap; expect failure for missing `vehicleConfig.parameterNames`. Define names by moving the exact current list, preserving order and spelling, then change only DOE discovery:

```matlab
function knownNames = parameterRegistry()
    knownNames = string(vehicleConfig.parameterNames()).';
end
```

- [ ] Run the new registry tests and the current DOE assertion scripts explicitly with `run(fullfile(modelRoot,'tests','test_doe_study_config.m'))` and `run(fullfile(modelRoot,'tests','test_doe_explicit_cars.m'))`. Confirm current error IDs and RNG preservation. Update the unknown-parameter message to name the registry without changing its identifier.
- [ ] Commit only the registry/discovery/tests and baseline tooling if commits are part of the approved execution workflow: `refactor: expose vehicle parameter registry`. Baseline artifact data remains outside tracked source.

## Task 2: Extract editable setup defaults and preserve sampling semantics

**Files:**

- Create: `+vehicleConfig/defaultSetup.m`, `fromLegacyGroups.m`, `sampleSetups.m` under `Full Car Models`.
- Create: `Full Car Models/tests/test_vehicleConfigSetup.m`.
- Read: `carConfig.m`, `parameters_loop.m`, `doeInitialDesign.m`.

**Interfaces:** Consumes `parameterNames`. Produces the canonical source struct described in the spec and `[setups,X]`; QSS `X` has exactly 57 ordered columns. These APIs do not yet replace existing callers.

- [ ] Write setup-copy, default-parity and sampling tests. Use the before-change reference from Task 1, rather than deriving expectations only from the new functions.

```matlab
function testSetupEditsAreIndependent(testCase)
    [base,eventParams] = vehicleConfig.defaultSetup();
    variant = base;
    variant.parameters.wheelbase = 1.60;
    variant.parameters.final_drive = 3.2;
    verifyEqual(testCase,base.parameters.wheelbase,62*0.0254);
    verifyEqual(testCase,base.parameters.final_drive,33/11);
    verifyEqual(testCase,eventParams.accel_length,75);
end

function testExplicitRowsKeepTheRequestedValues(testCase)
    base = vehicleConfig.defaultSetup();
    Q = table([155;170],[0.30;0.50], ...
        'VariableNames',{'mass','R_sf'});
    [setups,X] = vehicleConfig.sampleSetups(base,"Explicit",Q);
    verifySize(testCase,setups,[2 1]);
    verifyEqual(testCase,X.mass,[155;170]);
    verifyEqual(testCase,arrayfun(@(s)s.parameters.mass,setups),[155;170]);
    verifyEqual(testCase,X.Properties.VariableNames, ...
        vehicleConfig.parameterNames());
end
```

- [ ] Add cases for a zero-row Explicit table, a one-row table, unknown variable name, non-scalar Explicit baseline, string/char sampling names and a two-factor grid. The two-factor expected order is the exact existing `ndgrid` order, not a sorted order. `X` must be a 0-by-57 table for empty Explicit input.
- [ ] Run the new tests and verify missing APIs fail before implementation. Extract every numeric/event default verbatim. Resolve tyre MAT paths from the model directory, not the caller's working directory; retain explicit injected tyre arrays and torque functions in `fixed` for legacy adapters.
- [ ] Implement `fromLegacyGroups` by mapping the exact registry to the original six groups, retaining all non-sampled engine/tyre/map fields. Implement `sampleSetups` by moving the existing numeric sampler, then expanding each sampled row into an independent setup struct. Use one canonical names list, not several duplicated field registries.
- [ ] Preserve the distinction between direct legacy `parameters_loop` LHS/Random and `carConfig`'s deprecated DOE dispatch: the latter still calls `doeResolveStudy`/`doeInitialDesign` with `DOEStudyConfig.randomSeed`. Do not merge these two sampling algorithms. Verify reproducibility, ambient RNG behavior and the existing warning ID using the old route reference.
- [ ] Run `test_vehicleConfigRegistry`, `test_vehicleConfigSetup` and the relevant DOE assertion scripts. Compare fixed/default/event values to Task 1 snapshots. Commit: `refactor: separate vehicle setup defaults and sampling`.

## Task 3: Build one consistent car from a setup and shared assets

**Files:**

- Create: `Full Car Models/+vehicleConfig/loadAssets.m`, `buildCar.m`.
- Create: `Full Car Models/tests/test_vehicleConfigBuildCar.m`.
- Read/reuse: `utilities/buildCarFromParameterSet.m`, `buildTireModel.m`, `aeroMapEnabled.m`, `Car.m`, `Powertrain.m`, `Aero.m`, `+rampSpeed/deriveRollStiffness.m`.
- Preserve: component class inheritance and positional constructor signatures.

**Interfaces:** Consumes a resolved canonical setup and optional asset context. Produces one scalar `Car` and a struct of derived construction values. Handles QSS lap/acceleration now; ramp policy is added in Task 5 without changing this signature.

- [ ] Write consistency tests that edit source setup values together and rebuild:

```matlab
function testGeometryAndGearingAreRebuilt(testCase)
    base = vehicleConfig.defaultSetup();
    edit = base;
    edit.parameters.wheelbase = 1.60;
    edit.parameters.weight_dist = 0.55;
    edit.parameters.wheel_radius = 0.205;
    edit.parameters.final_drive = 3.2;
    car = vehicleConfig.buildCar(edit,"lap");
    verifyEqual(testCase,car.W_b,1.60,'AbsTol',1e-12);
    verifyEqual(testCase,car.l_f,1.60*0.55,'AbsTol',1e-12);
    verifyEqual(testCase,car.l_r,1.60*0.45,'AbsTol',1e-12);
    verifyEqual(testCase,car.R,car.powertrain.wheel_radius);
    pt = car.powertrain;
    expected = pt.shift_point./pt.gears/pt.final_drive/ ...
        pt.primary_reduction*pi/30*car.R;
    expected(end) = pt.redline/pt.gears(end)/pt.final_drive/ ...
        pt.primary_reduction*pi/30*car.R;
    verifyEqual(testCase,pt.switch_gear_velocities,expected,'AbsTol',1e-12);
    verifyEqual(testCase,base.parameters.weight_dist,0.512);
end

function testSpringEditUpdatesWheelRateWithoutRetuningQssSplit(testCase)
    setup = vehicleConfig.defaultSetup();
    setup.parameters.spring_rate_front_lb_in = 350;
    car = vehicleConfig.buildCar(setup,"lap");
    verifyEqual(testCase,car.rideHeightAero.wheel_rate_front_Npm, ...
        350*175.126835*setup.parameters.motion_ratio_front^2,'AbsTol',1e-9);
    verifyEqual(testCase,car.R_sf,setup.parameters.R_sf);
end
```

- [ ] Add invalid cases: zero/negative wheelbase or radius, non-finite/complex/scalar-array source input, rear weight fraction outside (0,1), non-positive spring rates/MRs, unknown role/profile/schema key and vector setup supplied to single-car builder. Signed camber/toe/ride-height values and zero inertia/Crr/compliance remain valid. Use specific errors such as `vehicleConfig:invalidParameter`, `vehicleConfig:badRole`, `vehicleConfig:badSchema` and `vehicleConfig:unknownField`.
- [ ] Add construction parity for QSS lap and acceleration roles: 226/221 kg, exact aero fields, lap map present/acceleration map absent, explicit `R_sf`, both tyre implementations and LC0 uncertainty. Load and compare the entire Task 1 snapshots, not just mass and wheelbase.
- [ ] Add a rebuild test where a previously built car has sentinel `ggPoints`, lookups and `comp`; a car rebuilt from its original setup must have empty results and an empty `ggMask` struct. Do not copy any `Car` result fields into the source configuration. Verify building a second car does not mutate the first.
- [ ] Run tests and establish missing builder failures. Implement absolute-path asset resolution and explicit shared context; injected full `AeroMap` corrections and camber assets must be used without reloading or losing provenance. Do not apply a context to a setup with a different asset identity: reload for a new identity or fail explicitly with `vehicleConfig:assetMismatch`.
- [ ] Implement the builder through existing component constructors and `buildCarFromParameterSet`. Derive fields once during construction. Do not add new property setters or calls on `Car.equations`' hot path. Pass every existing grip/inertia/Crr/camber/steering input.
- [ ] Run focused builder, tyre and ride-height construction tests. For deterministic direct equation evaluations, compare against Task 1 results at preselected physically valid states; use exact finite masks and existing test tolerances. Do not run `gg2` as part of this task's construction-only tests. Commit: `feat: build independent cars from editable setups`.

## Task 4: Add named cases and route legacy generation through the new core

**Files:**

- Create: `Full Car Models/+vehicleConfig/buildCases.m`, `toLegacyCells.m`.
- Modify: `Full Car Models/carConfig.m`, `utilities/parameters_loop.m`.
- Create: `Full Car Models/tests/test_vehicleConfigCases.m`, `test_carConfigCompatibility.m`.
- Preserve: `runAdaptiveDOE.m`, `doeRunBatch.m`, `doeRunCase.m` and their public formats.

**Interfaces:** Consumes the setup/sampling/builder APIs. Produces named cases or compatible cells. The fifth `carConfig` input is `outputFormat`, default `"cells"`; the other inputs/outputs remain unchanged.

- [ ] Write conversion and opt-in tests before routing production:

```matlab
function testNamedCasesMatchTheLegacyOutput(testCase)
    Q = table([155;170],[0.30;0.50], ...
        'VariableNames',{'mass','R_sf'});
    [legacy,E,X,B] = carConfig("Explicit",Q);
    [cases,E2,X2,B2] = carConfig("Explicit",Q,"legacy","nominal","cases");
    verifySize(testCase,cases,[2 1]);
    verifyEqual(testCase,string({cases.id}).',["case-000001";"case-000002"]);
    verifyEqual(testCase,[cases.designRow].',[1;2]);
    verifyEqual(testCase,E2,E);
    verifyEqual(testCase,X2,X);
    verifyEqual(testCase,B2,B);
    packed = vehicleConfig.toLegacyCells(cases);
    verifySize(testCase,packed,[2 2]);
    for i = 1:2
        verifyEqual(testCase,carConfigurationSnapshot(packed{i,1}), ...
            carConfigurationSnapshot(legacy{i,1}));
        verifyEqual(testCase,carConfigurationSnapshot(packed{i,2}), ...
            carConfigurationSnapshot(legacy{i,2}));
    end
end
```

- [ ] Add default calls with one/two/three/four outputs; `FullFactorial`/`Grid` aliases; char/string selectors; Explicit empty and one-row tables; LHS/Random warning and seed behavior; invalid output format; two-column role ordering; and constructor failures with existing wrapper IDs. Source setups must correspond to all 57 columns in the same row.
- [ ] Add a conversion test that assigns sentinel solved G-G arrays and a live `Events2` object to a named `lapCar`; packing must preserve arrays and the same event-result handle. Constructing fresh cases must not share `comp` handles or results across rows. Save/load empty, singleton and two-row named cases and legacy cells using MAT files in a temporary test directory.
- [ ] Implement case construction with a correctly shaped empty template and shared asset context:

```matlab
template = struct('id',"",'label',"",'designRow',NaN, ...
    'setup',struct(),'lapCar',[],'accelCar',[]);
cases = repmat(template,numel(setups),1);
for i = 1:numel(setups)
    cases(i).id = "case-" + compose('%06d',i);
    cases(i).label = setups(i).label;
    cases(i).designRow = i;
    cases(i).setup = setups(i);
    cases(i).lapCar = vehicleConfig.buildCar(setups(i),"lap",assets);
    cases(i).accelCar = vehicleConfig.buildCar(setups(i),"accel",assets);
end
```

Only share a context across setups with the same asset identity. A list containing different identities requires separate contexts per identity; do not silently reuse the first map/model.

- [ ] Convert `parameters_loop` into `fromLegacyGroups -> sampleSetups -> buildCases -> toLegacyCells`, preserving its direct sampling and legacy errors. Change `carConfig` to get defaults, resolve baseline/design rows, retain its DOE warning/dispatch, build named cases, and pack only if requested. Compute the baseline table without constructing/discarding an extra pair of cars on Explicit/DOE requests; verify this removes no relied-on error or warning.
- [ ] Run all new config tests, `test_entrypointBootstrap`, `test_camberKinematicsConfig`, current aero-map star construction tests and `test_doe_explicit_cars.m`. Compare all Task 1 construction references. Commit: `feat: add named car cases with legacy-compatible configuration`.

## Task 5: Integrate the independent frozen ramp setup path

**Files:**

- Create: `Full Car Models/+vehicleConfig/fromRampSetup.m`.
- Modify: `Full Car Models/+vehicleConfig/buildCar.m`, `utilities/buildSingleRampCar.m`.
- Create: `Full Car Models/tests/test_vehicleConfigRampCompatibility.m`.
- Read/preserve public behavior: `carConfigBaseline.m`, `+rampSpeed/buildCarFromSetup.m`, `buildSetupCatalog.m`, `normalizeSetupSpec.m`, `loadRampAssets.m`, `applySolverProfileToCar.m`.
- Existing tests: `test_carConfigBaseline`, `test_rampSpeedSetupModel`, `test_rampSpeedDriverWeight`, `test_rampSpeedSetupUnits`, `test_rampSpeedSetupPersistence`, `test_rampSpeedSetupCatalog`, `test_rampSpeedAssetResolution`.

**Interfaces:** `fromRampSetup(fixed,setup,mapInput,assets)` produces the canonical `"ramp"` configuration. `buildSingleRampCar` retains `[car,derived]`; `rampSpeed.buildCarFromSetup` retains `[car,setup,derived]`. Root `carConfigBaseline` remains independent of mutable QSS defaults.

- [ ] Write parity tests with Task 1 ramp snapshots, including a changed driver weight/rear distribution, alternative allowed springs/rear ARB, selected map identity and axle datums.

```matlab
function testRampBuildRetainsTheFrozenBaseline(testCase)
    [expected,config] = carConfigBaseline();
    setup = config.defaultSetup;
    assets = rampSpeed.loadRampAssets(config,setup);
    spec = vehicleConfig.fromRampSetup(config.fixedParameters,setup,assets.map,assets);
    [actual,derived] = vehicleConfig.buildCar(spec,"accel",assets);
    verifyEqual(testCase,carConfigurationSnapshot(actual), ...
        carConfigurationSnapshot(expected));
    verifyEqual(testCase,actual.M,221);
    verifyTrue(testCase,actual.aero.hasMap());
    verifyEqual(testCase,actual.R_sf,derived.R_sf,'AbsTol',1e-12);
end
```

- [ ] Test that ramp conversion fills every required field from its own frozen inputs; change a QSS setup and verify rebuilding the ramp baseline is unchanged. Verify ramp `"lap"` is rejected, legacy QSS `"accel"` has no map, and ramp `"accel"` uses its selected map. Preserve tyre model as legacy for existing ramp baseline.
- [ ] Test derived roll split and wheel rates after spring/ARB edits, all existing `derived` output field names, transient fields (`k_rf`, `k_rr`, `k_f_b`, `k_r_b`, `k_f_r`, `k_r_r`, `rs_total`), exact kg/percent/inch/unit conversions, selected map fingerprints and setup metadata. Out-of-envelope ramp heights must preserve the current warning/nearest-extrapolation policy, not adopt the ride-height sweep's stricter rejection policy.
- [ ] Run failures for missing ramp adapter/unsupported builder policy. Implement the adapter by translating current `fixed.vehicle` plus normalized setup values; do not call QSS `defaultSetup` to fill missing ramp fields. Reuse `deriveRollStiffness` and `loadRampAssets` data. Preserve `rs_total`'s Nm/degree conversion.
- [ ] Thin `buildSingleRampCar` into the adapter plus core builder, preserving existing validation and optional `mapInput`/`assets` behavior. Leave UI-owned setup normalization/discrete choices/map warnings at `rampSpeed.buildCarFromSetup`. Do not change `RampSpeedSession` or the App Designer layout as part of construction integration.
- [ ] Run the focused new and existing ramp construction/setup/persistence tests. Compare the Task 1 reference, including protected baseline independence. Commit: `refactor: share car construction without changing ramp setup contracts`.

## Task 6: Verify launcher, sweep, checkpoint and viewer compatibility

**Files:**

- Create: `Full Car Models/tests/test_carConfigEntrypointCompatibility.m`.
- Create: `Full Car Models/tests/helpers/runCarConfigurationChecks.m` (construction-only driver).
- Update: `Full Car Models/tests/car_configuration_entrypoints.csv` with evidence.
- Create: `Full Car Models/entrypoints/car_setup_example.m`.
- Modify only where baseline triage proves necessary: `setup_paths.m`, `entrypoints/bootstrap.m`, `tests/test_entrypointBootstrap.m`, `tests/test_steadyStateLapsim.m`, and launchers with refactor-caused failures.
- Read/preserve: `sweeps/aeroStarCases.m`, `aeroMapStarCases.m`, `gripSweep.m`, `loadAeroCache.m`, `utilities/loadDOECheckpoint.m`, `saveDOECheckpoint.m`, `doeAnalyze.m`, `simViewer.m`, `exportQSSResults.m`, `runAdaptiveDOE.m`, `entrypoints/run_sim_viewer.m`, `entrypoints/LapsimSweepComparisonPlotter.m`.

**Interfaces:** Existing callers continue receiving cell arrays and saved `carCell`; no named-case migration is required in their public interfaces. New runnable example uses configuration editing and named cases but performs no solver run.

- [ ] Add tests for bootstrapping from an unrelated temporary current directory, correct `Car`/`Tire2`/`vehicleConfig.buildCar` resolution, case-role ordering, and both ramp entrypoint invocation styles. Restore path/current directory with `onCleanup`; use fresh MATLAB processes for scripts that call `clear classes` so they cannot destroy a test runner's objects.
- [ ] Add sweep construction checks showing baseline/variants remain independent, mapped sensitivity retains `AeroMap.withCorrections`, static sensitivity retains its current static `Aero` replacement policy, and untouched acceleration variants remain untouched. Do not unify these two sensitivities' aero policies.
- [ ] Build a legacy MAT test fixture with `carCell`, `plan`, `designTable` and attached solved arrays/results. Verify load/selection through cache/checkpoint/analysis helpers, packing into a named case does not discard results, and current saved variable names remain unchanged. Historic real MAT files are read-only validation inputs.
- [ ] Keep adaptive-grip injection hooks `options.carConfigFcn`, `options.baseCell` and `options.sweepFcn`; rerun `test_adaptive_grip_sweep`. Its baseline mask still uses exactly one zero/zero ride-height row. Avoid replacing it with frozen ramp datum selection or a scalar Explicit call that fails on a vector baseline.
- [ ] Run the construction-only driver against each manifest row. This verifies configuration compatibility, not full launcher execution; write status `construction-pass` where that is all the evidence supports. Analyze dependency lists and `checkcode` separately; source checks alone are not a pass for solver/UI behavior.
- [ ] Add the minimal example:

```matlab
entrypointDir = fileparts(mfilename('fullpath'));
run(fullfile(entrypointDir,'bootstrap.m'));
[setup,eventParams] = vehicleConfig.defaultSetup();
variant = setup;
variant.parameters.wheelbase = 1.60;
variant.parameters.final_drive = 3.2;
car = vehicleConfig.buildCar(variant,"lap");
Q = table([155;170],'VariableNames',{'mass'});
[cases,eventParams,designTable] = ...
    carConfig("Explicit",Q,"legacy","nominal","cases");
lapCar = cases(2).lapCar;
```

- [ ] Resolve baseline path/test discrepancies narrowly. For SteadyState, keep the runtime tyre selection unless the user explicitly directs a model change; update stale comments/test paths to the observed canonical entrypoint rather than reviving deleted root files. For same-named ramp paths, preserve both behaviors and prove function resolution and explicit `run(scriptPath)` work. For `RH.m`, preserve user settings and add bootstrap only if its expected fresh-session launch needs it, after coordinating the untracked file.
- [ ] Record external/historical failures, e.g. missing `fmincon inputs.mat` or undeclared `car.c`, as baseline issues. Repair a clearly mechanical launch error only when needed for a required supported row; do not invent a damping parameter or alter physics to make a test pass. Such an unresolved row prevents an unconditional all-entrypoints success claim.
- [ ] Run construction/preview/cache compatibility tests and existing offline DOE/telemetry/plotting tests on temporary fixtures. Commit: `test: cover vehicle configuration entrypoint and cache compatibility`.

## Task 7: Update release packaging and user documentation

**Files:**

- Modify: `Full Car Models/RampSpeedApp.prj`, `+rampSpeed/releaseManifest.m`, `tests/test_rampSpeedReleaseManifest.m`, `entrypoints/README.md`, `README_ramp_speed_app.md`.
- Create: `Full Car Models/README_car_configuration.md`.
- Inspect: `RampSpeedApp.mlapp` for builder/baseline calls and package dependency discovery; avoid editing the archive unless a dependency change requires it.

**Interfaces:** Released project/package includes `+vehicleConfig` and required root configuration functions; users retain their launch commands and ramp study schemas.

- [ ] Write a release-manifest test asserting `+vehicleConfig` plus `carConfig.m`, `carConfigBaseline.m` and root `runRampSpeedStudy.m` are inventoried as required runtime inputs. The project descriptor must include the new package without changing `RootFolder` from `.`. Existing `releaseManifest` omits root constructor/configuration inputs; establish whether dependency discovery already includes them, then make the manifest explicit.
- [ ] Implement project/manifest additions and verify a copied release directory can resolve the package and construct QSS/ramp baselines without falling back to `C:\VD`. Copy only declared release inputs; test class/function paths point inside the copy.
- [ ] Document single-car editing/rebuilding, named cases, legacy default cells, value-class copy behavior, stale-result avoidance, units, separate ramp baseline and map policy, and MAT persistence limitations. Replace instructions to edit defaults inside `carConfig.m` with editable setup/Explicit design examples where relevant; default configurations now live in `defaultSetup`.
- [ ] Inspect the `.mlapp` archive read-only for car-construction references; include the package through the project/dependency mechanism. If archive source must change, preserve unrelated edits and perform a post-package reopen/instantiate. A zip edit alone is not validation.
- [ ] Run `test_rampSpeedReleaseManifest`, `test_rampSpeedAppSmoke`, `test_rampSpeedAppRequest`, `test_rampSpeedAppLifecycle`, `test_rampSpeedSetupEditorApp` and session tests under a desktop-capable MATLAB environment. Reopen `RampSpeedApp`, confirm baseline/setup edit/duplicate/create-run behavior and close cleanly. Report a headless/desktop limitation explicitly.
- [ ] Commit documentation/package changes: `docs: document editable car setups and package the shared factory`.

## Task 8: Complete bounded numerical validation and handoff

**Files:**

- Create: `Full Car Models/tests/helpers/runCarConfigurationIntegrationChecks.m` (explicit opt-in solver driver).
- Update: manifest and execution report under `artifacts/car_configuration_refactor/<run-id>/`.
- Review: whole diff against the spec; no new physics changes or broadened tolerances.

**Interfaces:** Consumes before-change evidence and final public APIs; produces full command logs, numerical comparisons, runtime costs and per-entrypoint statuses. This is the final compatibility gate, not a performance optimization project.

- [ ] Use the before-change numerical references captured in Task 1, or the preserved coordinated before-change snapshot if baseline simulations were deferred. Separate reference/after-change runs into fresh MATLAB processes with identical inputs, solver profiles/grid settings/seeds and workers. Verify the reference predates construction routing changes. Preserve originals if capture fails; numerical parity remains unverified until a valid reference exists.
- [ ] Run one scalar QSS baseline through existing G-G generation and representative selected events; run legacy plus an LC0 construction/evaluation case when its assets exist. Retain production cap/guards and outer DOE/inner serial parallelism. Compare statuses, finite/NaN masks, G-G fields, times, points and telemetry arrays to the reference; use `1e-8` for event-time/score parity where existing DOE tests already do, `1e-12` for deterministic scalar construction and existing solver tests' tolerances for converged fields. Investigate any tolerance failure instead of changing the limit.
- [ ] Run one two-row Explicit case study and a two-worker comparison where PCT is available. Record serial fallback as skipped parallel evidence when unavailable. Adaptive DOE/checkpoint/backfill tests must retain case indices, scores, errors and resume behavior; inspect `ramp_backfill_failures`, not only `status == "complete"`.
- [ ] Run a small ramp baseline and an edited setup with one speed in coast and balanced modes, compare capability/balance/camber/ride-height fields and solver statuses/masks, including linear-region steering outputs. Use current `test_rampSpeedEndToEndRealCar`, scheduler/solver-profile tests and meaningful existing adapters; do not launch the full default multi-setup study merely for confidence.
- [ ] Execute canonical numerical launchers under their explicitly configured defaults where affordable; otherwise use actual entrypoint function options for bounded designs. A construction check or a call to a lower-level helper does not count as executing the full script. Record full-launch status separately and report any default run deliberately deferred for cost/input availability.
- [ ] Save/reload one sweep MAT file and one DOE checkpoint; open the existing MATLAB viewer on the solved car; export baseline QSS HDF5 and verify the existing loader/schema tests. Offline comparison/DOE fitting launchers need valid current fixtures and must not trigger hidden simulations.
- [ ] Complete every validation matrix row with `full-pass`, `construction-pass`, `baseline-failure`, `blocked-input` or `not-run`, exact commands and evidence paths. Canonical supported launchers require `full-pass` or an explicitly reported external blocker before claiming broad entrypoint compatibility. Preserve remaining historical failures in the report.
- [ ] Run final required tests once after the last changes; do not rerun costly simulations repeatedly without a new change/failure. Review the full branch diff for source-parsing remnants, stale duplicated defaults, unintentional tyre switches, omitted constructor inputs, cache-schema drift, dirty-file overwrites and added hot-path I/O.
- [ ] Final report: changed API, example launch commands, checkout/branch/commit/dirty state, test counts, numerical deltas, wall time and outputs, and exact incomplete/blocking rows. Integrate only through the execution method chosen by the user.

## Entrypoint compatibility matrix

Each path is explicit. A row may group closely related launchers, but the CSV must contain a separate record and evidence for each file. No entrypoint is removed or renamed in this plan.

| Path under repository root | Contract to preserve | Required evidence |
| --- | --- | --- |
| `Full Car Models/entrypoints/bootstrap.m` | Portable model/repository root and class paths | Run from a temporary foreign cwd; verify resolved functions/classes and restore cwd |
| `Full Car Models/setup_paths.m` | Correct class resolution; shadow/stale-class diagnostics | `which -all` for core classes and package; project-copy resolution |
| `Full Car Models/carConfig.m` | 0-4 input compatibility, all four outputs, 57 columns, both roles | Baseline/Explicit/DOE parity plus new opt-in named output |
| `Full Car Models/carConfigBaseline.m` | One protected acceleration-oriented ramp car/config | Exact baseline/axle-datum/setup/asset parity |
| `Full Car Models/entrypoints/SteadyStateLapsim.m` | Current tyre selection, selected events, G-G cap/workers, per-case writeback/logging | Baseline triage; bounded actual event run; no incidental legacy/LC0 switch |
| `Full Car Models/entrypoints/run_single_baseline_qss.m` | Unique zero/zero baseline, G-G before Autocross, HDF5 export | Actual single-baseline run/export/load and baseline-mask negative cases |
| `Full Car Models/entrypoints/AeroSensitivityStudy.m` | Static-aero star design and acceleration baseline, cache/results | Construction parity and one bounded study/cache cycle; preserve static-map replacement policy |
| `Full Car Models/entrypoints/AeroMapSensitivityStudy.m` | Mapped ride-height/global correction star cases; map fingerprint | Independent map-corrected cases and bounded study/cache cycle |
| `Full Car Models/entrypoints/run_grip_sweep.m` | Front/rear changes applied to both roles; balance/target/cache | Grid construction and bounded existing gripSweep solve/score; full script status recorded separately |
| `Full Car Models/entrypoints/run_adaptive_grip_sweep.m` | Baseline mask, dependency injection, refinement/confirmation | Existing injected tests plus small actual function-option run |
| `Full Car Models/entrypoints/run_ride_height_stiffness_sweep.m` | Requested override subset retained, map bounds, Explicit role pair/output | Preview checks and existing one-case MAT/CSV output test; no default 225-case run |
| `Full Car Models/entrypoints/RH.m` | User's options/output path forwarded | Path-resolution/preview of same options; preserve existing file and avoid default 225-case launch during light checks |
| `Full Car Models/entrypoints/runRampSpeedStudy.m` | Legacy script consumes lap-car list; cache/plot output | Explicit `run` path and representative ramp result/plot; shadowing check |
| `Full Car Models/runRampSpeedStudy.m` | Options function delegates to frozen ramp setup/executor | Function resolution from foreign cwd and bounded actual function run/cache reload |
| `Full Car Models/RampSpeedApp.mlapp` | Frozen setup catalogue, independent duplicates, builder/executor and UI | Archive dependency inspection; actual instantiate/edit/duplicate/run/close |
| `Full Car Models/runAdaptiveDOE.m` | Explicit design/role pair, global indices, checkpoints/resume/ramps | DOE registry/config tests; small case/checkpoint/backfill evidence |
| `Full Car Models/entrypoints/run_sim_viewer.m` | Legacy `carCell` case loading, equality cache gate and event defaults | Load real solved MAT/cache; actual viewer open and stale-cache rejection |
| `Full Car Models/entrypoints/LapsimSweepComparisonPlotter.m` | Historical `carCell` MAT concatenation and public points | Fixture plus historical-file load/plot if available; report missing default file |
| `Full Car Models/entrypoints/DOE_Fitting.m` | Offline DOE result analysis and optional ramp backfill | Existing offline result/fixture run with `runRampMetrics=false`; unchanged public columns |
| `Full Car Models/entrypoints/wheelbase_sweep.m` | Baseline extraction, tyre/aero/geometric analysis | Baseline property parity and actual bounded analysis/plot; preserve current analytical aero edits |
| `Full Car Models/entrypoints/UndersteerGradientSim.m` | Car property/method contract, case loops and plots | Construction plus representative steady-state solve; full script status separately |
| `Full Car Models/entrypoints/ModalAnalysis.m` | Independent lap car augmented with transient properties and damping assets | Current damping asset resolution and actual modal analysis/figure generation |
| `Full Car Models/entrypoints/ModalAnalysisDecoupled.m` | Same constructor/transient contract with decoupled stiffness | Actual matrix/modal calculation and plots against baseline |
| `Full Car Models/entrypoints/transientSim.m` | Lap car transient fields, mass/geometric properties and dynamics | Constructor/field parity plus bounded transient integration and output masks |
| `Full Car Models/CnAy/CnAy.m` | Existing carCell role access and constant-steer/beta APIs | Correct model path, baseline parity, representative point; full solve grid status separately |
| `Full Car Models/Phase Plane/PhasePortrait.m` | Value car copy and `PhasePlaneODE` access | Method/field parity plus short actual ODE trace |
| `Full Car Models/Suspension Transient Model/FullCarVibrationModel.m` | Historical Car field contract | Record baseline `car.c` discrepancy and dependency resolution; no invented damping fix |
| `Full Car Models/testing/test.m` | Historical fmincon-input workflow and editable powertrain | Fixture availability, independent car edit and representative solve; report absent MAT input |
| `Other/AeroCorneringSim.m` | Baseline scalar Car and optimisation methods | Resolve main `Car`/`Tire2`, build parity, representative point |
| `Other/DiffAccelSim.m` | Same with powertrain/differential edits | Same plus independent torque-transfer edit and point |
| `Other/TorqueSpeedAnalysis_aeden.m` | Existing acceleration/tyre/component methods | Correct main-model paths, constructor parity and point |
| `Other/tire_load_calculator.m` | Existing tyre/load-transfer car interface | Correct main-model paths, constructor parity and point |
| `Magic Formula/experimental/lc0_nd_tire/run_lc0_nd_pacejka_comparison.m` | Optional legacy carConfig fallback and tyre comparison | Preserve constructor fallback; existing artifact plotting/tyre tests; report absent comparison assets |
| `Full Car Models/entrypoints/car_setup_example.m` (new) | Construction-only editable setup/named-case demonstration | Run from foreign cwd; consistent geometry/gearing and correct named roles |

At execution start, repeat the repository-wide search below and extend the manifest for new/untracked scripts, `.mlx`, `.mlapp` or Simulink consumers that invoke the factory indirectly. Archived `OLD` and intentionally separate BDI classes are not silently rewritten; list them as historical/separate-model consumers if dependency inspection finds a link.

## Exact verification commands and test tiers

These commands are for implementation, not instructions to run simulations during planning. Use the chosen checkout root instead of `C:\VD` if different.

Inventory (PowerShell):

```powershell
git -C C:\VD status --short --branch
rg -n '\b(carConfig|carConfigBaseline|parameters_loop|buildSingleRampCar|buildCarFromParameterSet)\s*\(' C:\VD -g '*.m' -g '*.md' -g '!OLD/**'
rg --files C:\VD -g '*.mlapp' -g '*.mlx' -g '*.slx'
```

Construction-only driver (PowerShell; the helper is created in Task 6):

```powershell
& 'C:\Program Files\MATLAB\R2026a\bin\matlab.exe' -batch "run('C:/VD/Full Car Models/tests/helpers/runCarConfigurationChecks.m')"
```

Driver contents include model bootstrap, helper path, the new configuration test files and the verified construction/preview suites. Assert no test unexpectedly entered `gg2`, event solvers or ramp execution. Do not use a wildcard `runtests('tests')` for this tier: existing tests include real simulations, and existing `test_doe_*.m` files include assertion scripts that need `run`.

Example precise function-test invocation inside the driver:

```matlab
modelRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
run(fullfile(modelRoot,'entrypoints','bootstrap.m'));
addpath(fullfile(modelRoot,'tests','helpers'));
files = {'test_vehicleConfigRegistry','test_vehicleConfigSetup', ...
    'test_vehicleConfigBuildCar','test_vehicleConfigCases', ...
    'test_carConfigCompatibility','test_vehicleConfigRampCompatibility', ...
    'test_carConfigBaseline','test_rampSpeedSetupModel', ...
    'test_rampSpeedDriverWeight','test_rampSpeedSetupUnits'};
suite = matlab.unittest.TestSuite.empty;
for i = 1:numel(files)
    suite = [suite testsuite(fullfile(modelRoot,'tests',[files{i} '.m']))];
end
results = run(suite);
disp(table(results));
assert(~any([results.Failed]));
```

The helper's path calculation must resolve `Full Car Models` from `tests/helpers`, as shown. Additional preview tests must select only preview methods: `testRideHeightStiffnessSweep` also contains a real solver run. Restore cwd after assertion scripts with an `onCleanup` in their caller. A fresh MATLAB invocation is the boundary for launchers that clear classes.

Bounded numerical integration (explicit implementation gate):

```powershell
& 'C:\Program Files\MATLAB\R2026a\bin\matlab.exe' -batch "run('C:/VD/Full Car Models/tests/helpers/runCarConfigurationIntegrationChecks.m')"
```

The integration driver executes approved bounded numerical cases, records inputs and elapsed times, runs current `test_doe_parallel_batch.m` only in its designated expensive tier, and reports toolbox skips separately. UI tests use a desktop-capable MATLAB session. Do not classify skipped/assumed tests or source-contract assertions as full execution passes.

## Sequencing, rollback and acceptance

Tasks 1-5 are sequential because they share a registry, schema and constructor interfaces. Task 6 then verifies unchanged consumers; Task 7 protects packaged operation; Task 8 is the final numerical gate. The project remains callable through the old API after every production task. Commits are scoped so a routing change can be reverted without deleting the additive configuration APIs.

The completed refactor must demonstrate: editable setup isolation; consistent rebuilt geometry/gearing/suspension; exact default/table/role parity; named and legacy shapes; unchanged baseline/tyre/map policies; compatible saved results/checkpoints/export; correct package/path resolution; and per-entrypoint evidence. No default output-format removal is scheduled.

Implementation recommendation: **native execution** for the tightly coupled schema/builder/adapter tasks, followed by a review of the whole diff and compatibility evidence. Subagent-driven execution remains an option if the user chooses it; implementation and review must share the same spec, registry/interface contracts and coordinated dirty-state snapshot.
