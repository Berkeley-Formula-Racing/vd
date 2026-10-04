# Car configuration refactor design

Date: 2026-10-04. Status: planning artifact for user review; no implementation or simulation has been performed.

## Intent and success criteria

Provide an editable, independent configuration for one vehicle, build a consistent solver-ready `Car` from that configuration, and produce named collections of lap/acceleration cars for studies. Existing MATLAB launch commands, direct constructors, default `carConfig` callers, saved `carCell` results, ramp setup specifications and telemetry exports must keep working.

Reuse the existing `Car` value class. Scripts can edit setup structs directly and build a fresh car. This avoids changing the meaning of `variant = baseline` in sweeps. No second physics class or shared mutable handle model is required for this refactor.

## Current architecture, verified by source inspection

- `Full Car Models/carConfig.m` owns six parameter groups, event defaults, baseline construction and DOE dispatch. Its first output is an N-by-2 cell array of `Car` values, not raw cells of numbers.
- `utilities/parameters_loop.m` owns an ordered registry of 57 sampled numeric fields, sampling, asset loading and construction of both car roles.
- `utilities/buildCarFromParameterSet.m` is the shared low-level constructor. The `Car` constructor derives axle distances and roll-center height; `Powertrain` derives gear-change speeds.
- `carConfigBaseline.m` independently owns a frozen ramp baseline; `buildSingleRampCar.m` derives suspension stiffness. QSS uses a specified `R_sf`, whereas ramp derives it from springs and ARBs. These are intentionally different policies.
- Default QSS total masses are 226 kg lap and 221 kg acceleration. QSS static/reference ride heights are zero. Ramp baseline static/reference axle heights are approximately 3.86558910640757 and 5.76907823448674 inches. QSS lap cars carry the aero map; current QSS acceleration cars do not; ramp cars carry the selected map.
- `doeResolveStudy.m` discovers parameter names by parsing the source text of `parameters_loop.m`. It must migrate before that registry is moved.
- Viewers, analysis, checkpoints and sweep caches use `carCell`, table column order and `Car` properties. `Events2` is a handle object stored in `car.comp`; copying a solved `Car` does not guarantee independent event-result handles.
- `entrypoints/runRampSpeedStudy.m` is a legacy script; root `runRampSpeedStudy.m` is a different compatibility function. Both paths and their existing behavior must survive.

## Global constraints

- Keep `Car`, `Aero`, `Powertrain`, `Tire2` and `NondimensionalTire` as value classes; keep the positional `Car` constructor compatible.
- Preserve every existing `carConfig` call with zero to four inputs and one to four outputs; default output remains an N-by-2 cell array with lap then acceleration columns.
- Preserve the ordered 57-column QSS design-table schema, case ordering, baseline selection, DOE random-seed behavior and existing public warning/error identifiers.
- Preserve separate QSS and frozen ramp baselines, their aero-map role policies, ride-height datums, tyre model selection and roll-stiffness policies.
- Preserve `carCell` MAT/checkpoint payloads, ramp study schemas and telemetry HDF5/CSV schemas; do not convert existing files in place.
- Configuration edits followed by rebuilding must produce consistent derived fields and a fresh unsolved car; the supported edit path must not reuse stale G-G data or event-result handles.
- Load immutable map/camber assets once per batch or distinct asset identity; do not add file reads or configuration validation to solver evaluation loops.
- Do not run simulations during planning; implementation validation must use bounded runs and record which launchers were actually executed.
- Preserve unrelated dirty files and existing user launchers; no cleanup, baseline retuning or solver-physics changes are included.

## Public interfaces

Create `Full Car Models/+vehicleConfig` to avoid adding generic global function names or competing with `+rampSpeed` setup names.

```matlab
[setup,eventParams] = vehicleConfig.defaultSetup(tireModel,uncertaintyMode);
[car,derived] = vehicleConfig.buildCar(setup,role,assets);
[setups,designTable] = vehicleConfig.sampleSetups(baseSetup,samplingType,samples);
cases = vehicleConfig.buildCases(setups,assets);
carCell = vehicleConfig.toLegacyCells(cases);
[cases,eventParams,designTable,baselineTable] = ...
    carConfig(samplingType,samples,tireModel,uncertaintyMode,"cases");
```

Optional defaults are `"legacy"`, `"nominal"`, role `"lap"`, and an empty asset struct. `carConfig` gets a fifth optional `outputFormat` argument; `"cells"` is the default and `"cases"` opts into named cases. Its existing four outputs keep their order and meaning. A non-scalar or unknown output format is rejected as `carConfig:badOutputFormat`.

`defaultSetup` returns one scalar configuration struct and event defaults without constructing a `Car`. It is extracted from the current `carConfig` defaults without changing their values. Frozen ramp configuration stays owned by `carConfigBaseline`; it does not call `defaultSetup`.

### Configuration representation

The setup is a MATLAB struct, copyable and persistable with MAT files. It is not advertised as generally JSON portable because current torque functions and component objects have MATLAB-specific representations.

| Field | Contract |
| --- | --- |
| `schemaVersion` | Numeric scalar 1 |
| `profile` | `"qss"` or `"ramp"`; builder policy, not a solver switch |
| `id`, `label` | Non-empty scalar strings; sampled IDs generated in row order |
| `parameters` | Struct of source numeric parameters; QSS uses the exact registry names; a sampling template may contain vectors, but `buildCar` requires resolved scalars |
| `fixed` | Named non-sampled constructor inputs: engine gears, primary reduction, torque function, tyre fit arrays, friction scaling and tyre model/artifact/uncertainty |
| `aeroMap` | `enabled`, absolute `path`, `claScale`, `cdaScale`, `copOffset`; defaults true, existing map path, 1, 1, 0 |
| `suspension` | `rollSplitPolicy` (`"explicit"` or `"derived"`), front/rear ARB stiffness in Nm/rad; QSS retains explicit split; ramp retains derived split |

Do not store computed `l_f`, `l_r`, wheel rates, switching velocities or simulation results as editable source parameters. New configuration keys must be checked against the defined schema. Compatibility adapters retain their current validation order and identifiers where callers rely on them.

The supported editing pattern is:

```matlab
setup = vehicleConfig.defaultSetup();
variant = setup;
variant.parameters.wheelbase = 1.60;
variant.parameters.final_drive = 3.2;
car = vehicleConfig.buildCar(variant,"lap");
```

Existing raw `Car` assignments remain available for historical scripts, but they are not promoted as the consistent edit interface. This refactor does not install setters across the large physics class.

### Named case representation

`cases` is an N-by-1 struct array, including a correctly typed 0-by-1 empty value. Each case contains `id`, `label`, `designRow`, `setup`, `lapCar`, and `accelCar`. `setup.parameters` exactly matches the corresponding design-table row. IDs are deterministic, `"case-000001"`, `"case-000002"`, and so on, within a generated batch. They are not a replacement for adaptive DOE global case indices or ramp setup IDs.

`buildCases` builds fresh cars using an explicit shared asset context. `toLegacyCells` packs existing car values without rebuilding or removing G-G/competition data. It is one-way by design: importing legacy cells cannot reconstruct their complete source setup reliably. Existing cache readers remain legacy readers.

## Construction, assets and ramp integration

`vehicleConfig.buildCar` validates the resolved configuration, constructs role-specific `Aero`, `Powertrain` and tyre models, derives ride-height/camber context, then calls `buildCarFromParameterSet`. The low-level constructor and physics equations remain unchanged. Rebuilding therefore recalculates duplicated wheel radius, axle geometry, wheel rates and shift velocities and starts with empty `ggPoints`, steady-state information, lookups and `comp`.

Role is explicit: QSS `"lap"` uses `driver_weight`, lap aero and the configured map; QSS `"accel"` uses `accel_driver_weight`, acceleration aero and no map. Ramp uses the selected setup driver's mass, acceleration aero, selected map and derived roll stiffness. Accept ramp role `"accel"`; reject ramp `"lap"` rather than inventing a second ramp configuration.

`vehicleConfig.fromRampSetup(fixed,setup,mapInput,assets)` translates current ramp source inputs without reading mutable QSS defaults. `buildSingleRampCar` remains a public compatibility adapter and calls the new builder. It keeps its current derived-field names, assets/fingerprints and transient suspension assignments. `rampSpeed.buildCarFromSetup` continues to own discrete choices, setup normalization, map-envelope warnings and provenance.

`vehicleConfig.fromLegacyGroups(cP,aP,eP,DTp,Bp,tP)` translates direct `parameters_loop` inputs. `vehicleConfig.parameterNames` replaces source-text parsing. `vehicleConfig.loadAssets(setup)` resolves absolute file paths and loads map/camber data; builders accept an injected context so tests and ramp can reuse already loaded assets. Optional map corrections preserve the full `AeroMap` and its datum.

No persistent cross-study cache is necessary. An explicit batch asset context is sufficient; assets for different paths/models must never be silently reused. Existing `parameters_loop` missing-map identifiers and ramp asset/warning behavior remain at their wrappers.

## Migration strategy

1. Characterize current construction, entrypoint dependencies, saved formats and existing failures before editing production code.
2. Introduce the shared parameter registry and replace DOE source parsing while the old generator still works.
3. Extract defaults and add independent setup editing, single-car construction and named-case building.
4. Route `parameters_loop` and `carConfig` through the new core, retaining default cells and adding opt-in named output.
5. Route ramp construction through the core while retaining its independently owned baseline, setup schema and UI API.
6. Keep all existing launchers on the compatible cell API in this refactor. Add a documented single-car/named-case example; avoid a broad rewrite of unrelated analysis scripts.
7. Update release packaging and documentation, then complete the compatibility matrix in the plan.

Changing the default output type or rewriting every cache consumer is a possible later migration, after callers adopt named cases. It is not part of this release.

## Verification and completion claim

Every inventoried entrypoint gets dependency resolution, input/car construction, baseline/datum/model checks and a recorded validation status. Existing tests include function-based tests and executable assertion scripts; both must be accounted for. Focused construction/preview checks precede bounded real solver runs and UI checks.

Exact construction/table/role parity is required. Representative equation, G-G, event, ramp and telemetry results are compared to evidence captured before production edits; compare finite masks and solver status as well as numeric values. Do not loosen tolerances to accommodate a regression.

Pre-existing failures are recorded before editing and repaired narrowly only when necessary for the promised entrypoint compatibility. Comments or obsolete assertions must not select a different tyre model. A launcher with unavailable historical input data is marked unverified with its exact dependency, not passed. Do not claim that every launcher works until each required row has evidence or a clearly reported external blocker.

## Observed hazards requiring baseline triage

- `SteadyStateLapsim.m` currently calls `"legacy"`, while its comments and `test_steadyStateLapsim.m` expect LC0. The test also names a deleted root launcher. Establish the intended/current model before correcting test expectations; preserve the current runtime choice during this refactor.
- The two `runRampSpeedStudy.m` files have different semantics and can shadow each other. Verify `which -all` and test both invocation styles without deleting either.
- `entrypoints/RH.m` is an untracked six-line convenience script with a fixed output path and no bootstrap. Preserve it and include its current path assumptions in the matrix.
- `Suspension Transient Model/FullCarVibrationModel.m` assigns `car.c`, which is not declared in the inspected `Car` class. This is a pre-existing source-level mismatch, not a reason to make `Car` dynamically mutable.
- The working tree already contains unrelated ramp UI and test edits. Implementation must start from a coordinated snapshot/isolation strategy that includes required changes without overwriting them.
