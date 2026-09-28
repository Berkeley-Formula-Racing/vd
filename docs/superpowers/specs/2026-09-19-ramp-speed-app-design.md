# Visual Ramp-Speed Analysis App Design

**Date:** 2026-09-19  
**Status:** Implemented v1 design; canonical execution path and setup model documented below

## Goal

Build a shareable MATLAB App Designer application that runs the existing vehicle-dynamics ramp simulations and makes their speed-dependent capability, balance, aero, suspension, tire, and solver-health metrics easy to inspect and compare across multiple vehicle setups.

The first release supports two ramp types:

1. A lateral-limit ramp: solve a lateral-acceleration ramp at each `vCar` and report the free and sustainable lateral limits plus balance metrics.
2. A pure longitudinal ramp: solve maximum straight-line longitudinal acceleration at each `vCar` with `Ay = 0`.

## Context and existing interfaces

The current simulation workspace is MATLAB-first. The app should reuse these existing functions and data products rather than translating the vehicle model:

- `Full Car Models/sweeps/rampSweep.m` — lateral ramp solver and per-speed summaries.
- `Full Car Models/events/max_long_accel.m` and `long_accel_sweep.m` — pure longitudinal solver and speed sweep.
- `Full Car Models/carComponents/Car.m` — vehicle equations and named metric extraction.
- `Full Car Models/runRampSpeedStudy.m` — current multi-car lateral study/cache pattern.
- `Full Car Models/sweeps/plotRampSpeedStudy.m` and `plotRampSweep.m` — current plotting conventions.
- `Full Car Models/carConfig.m` — setup/case construction; column 1 is the lap car and column 2 is the acceleration car.
- `Full Car Models/tests/test_rampSweepAdditionalMetrics.m` and `test_plotRampSpeedStudy.m` — existing regression coverage.

The application must preserve the existing solver behavior and use adapters where the existing output shape is not yet uniform.

## Technology decision

Use MATLAB App Designer and a MATLAB Project. Do not build the first version in Python.

MATLAB is the correct foundation because the model, `Car` class, `fmincon` workflows, aero/ride-height coupling, `.mat` caches, and tests are already MATLAB code. The target users all have MATLAB, so sharing a project plus `.mlapp` is sufficient; MATLAB Compiler is an optional later packaging path, not a v1 dependency.

A Python UI over MATLAB Engine would introduce a second runtime and cross-language conversion. A full Python port would duplicate the coupled vehicle model and numerical validation suite. Both options increase risk without improving the requested MATLAB-user workflow.

## Scope

### In scope for v1

- Select one or more vehicle setups from `carConfig`/`designTable` or from a saved study.
- Run lateral-limit and pure-longitudinal ramps.
- Configure speed range, speed spacing, ramp resolution, solver tolerances, lateral coast/balanced mode, and optional parallel execution.
- Keep the UI responsive during runs with progress, elapsed time, current case/speed, and cancellation.
- Overlay results from multiple setups on common plots.
- Load and migrate the existing `ramp_speed_study.mat` result shape.
- Show raw ramp curves and a per-speed/per-point inspector.
- Export normalized tables and figures.
- Record run settings, setup identity, timestamps, MATLAB release, code/version fingerprint, and validity state.

### Explicitly out of scope for v1

- A Python rewrite of the vehicle model.
- A live Simulink control surface.
- A general optimization/DOE authoring environment.
- A combined lateral-longitudinal ramp as a third ramp type. The longitudinal mode is pure `Ay = 0` for v1; the existing combined solver remains available outside the app.
- Editing every `Car` property directly in the UI. Setup selection and existing `carConfig` cases are the initial source of truth.

## Ramp definitions

### Lateral-limit ramp

Call the existing `rampSweep` implementation through a thin adapter. For each speed, preserve both:

- `gLat_max`: the free lateral limit from `max_lat_accel`, without a longitudinal-acceleration constraint.
- `gLat_top`: the highest point actually reached by the selected lateral-ramp operating condition.

Expose the existing `coast` and `balanced` modes as an advanced operating-condition selector. `coast` is the compatibility default because the existing `runRampSpeedStudy.m` cache uses it; the UI must label the mode prominently because `balanced` holds longitudinal acceleration at zero and can terminate earlier at high speed.

Show `ramp_complete`, `power_limited`, `max_ceq`, `n_exit1`, and `n_exit2` next to the capability result. A truncated point is not silently treated as the free limit.

### Pure longitudinal ramp

Call `max_long_accel` over a configurable speed vector through a new adapter. The adapter must enforce and record:

- steer angle = 0
- lateral velocity = 0
- yaw rate = 0
- achieved `Ay = 0` within the recorded residual tolerance
- throttle and rear slip-ratio bounds from the existing solver

The result is a per-speed maximum `gLong`/`aLong` curve with powertrain, tire, aero, load-transfer, and solver-health fields. Lateral understeer/oversteer metrics are not applicable and must be represented as unavailable, never as numerical zero.

## Metric and sign conventions

All stored numerical values use SI units. Display conversions are presentation-only.

- **Mechanical balance:** front LLTD share, displayed as a percentage. This is a load-transfer distribution, not an understeer sign.
- **Aero balance:** front CoP/front downforce share, displayed as a percentage, alongside front and rear aero loads.
- **Handling balance:** `K_linear` from steer-excess versus lateral g, displayed in `deg/g`; positive is understeer and negative is oversteer.
- **Yaw-derived balance:** `CUOsteerFromYaw_linear_deg` plus selected-point and limit values; positive is understeer and negative is oversteer.
- **Aero force:** downforce is positive in the displayed downwards direction; resisting drag is positive. Use explicit labels `downforce` and `drag`, not an ambiguous signed `lift` label. If an upward-positive lift value is needed for export, expose it as `lift_N = -downforce_N` with that sign documented.
- **Capability:** lateral curves use `gLat_max` and `gLat_top`; longitudinal curves use pure `gLong_max`.
- **Loads:** show front/rear axle loads, all four wheel loads, minimum wheel load, load transfer, and a wheel-lift flag when the virtual wheel load is non-positive.
- **Validity:** show solver exit flag, maximum constraint residual, aero-map outside status, aero residual, truncated status, and power/traction limitation.

The unit selector must support at least m/s versus mph, g versus m/s², and N versus lbf. Stored SI values and exported SI columns remain unchanged.

## Application layout

Use one App Designer window with a persistent setup/run panel and analysis tabs.

### Run Setup

- setup source and case selector
- selected setup list with labels, source index, and lap/acceleration-car designation
- ramp type selector: lateral-limit or pure longitudinal
- lateral mode selector when lateral is active: coast or balanced
- speed start/stop/step or explicit vector
- ramp point count, top/bottom fractions, bisection count, and residual tolerance
- serial/parallel mode and worker count when Parallel Computing Toolbox is available
- Run, Cancel, Load Cache, Save Study, and Clear buttons
- progress log with case, speed, point count, elapsed time, and warnings

### Analysis tabs

- **Capability:** lateral free/sustainable limits or longitudinal maximum acceleration versus `vCar`.
- **Balance:** mechanical balance, aero balance, `K_linear`, yaw-derived balance, grip balance, axle-slip balance, and zero reference lines.
- **Aero & Loads:** total/front/rear downforce, drag, `ClA`, `CdA`, `L/D`, CoP, axle loads, four-wheel loads, and minimum wheel load.
- **Suspension:** front/rear ride height, pitch, roll, camber, shock travel, and aero-map status.
- **Raw Ramp:** lateral steer-excess versus `gLat`, axle-slip balance versus `gLat`, axle tire utilization, and load-transfer behavior for each selected speed.
- **Inspector/Data:** selected setup, speed, ramp point, complete named state, tire forces/slips/cambers, residuals, and exportable table.

Plots must use a central catalog so titles, units, signs, valid ramp types, and field names cannot drift between tabs. Overlay colors identify setups; line styles distinguish derived series such as free limit versus sustainable top.

## Overlay and comparison behavior

Each run carries a comparable case record containing setup label, source case, ramp type, lateral mode, speed grid, solver settings, run timestamp, and validity summary.

- Multiple setups may be selected for overlay on any compatible plot.
- Mixed lateral and longitudinal selections are allowed in the study but are filtered to type-appropriate plots; the UI shows a compatibility warning instead of drawing meaningless series together.
- Different speed grids are plotted at their native points. A numeric delta view uses an explicit comparison grid over the common speed domain and records the interpolation method.
- Invalid or truncated points remain visible as gaps or flagged markers. They are not silently interpolated into the primary result.
- The comparison view should support baseline selection and show `variant - baseline` deltas for selected metrics.

## Data model and cache

Use a versioned study struct saved as a `.mat` file:

```matlab
study.schemaVersion
study.created
study.appVersion
study.cases(i).id
study.cases(i).label
study.cases(i).source
study.cases(i).designRow
study.cases(i).carRole
study.runs(i).type
study.runs(i).mode
study.runs(i).settings
study.runs(i).perSpeed
study.runs(i).points
study.runs(i).runMeta
study.runs(i).status
```

`perSpeed` is one row per speed and contains the plot/overlay metrics. `points` contains the detailed solved operating points. `runMeta` contains solver health, provenance, warnings, and completion state.

The app must accept the existing `study.results`/`study.labels` shape from `ramp_speed_study.mat` through a migration function that fills missing schema fields without changing the original file.

## Execution and failure handling

The UI layer must not call solver internals directly. A runner owns execution and returns progress events plus result structs.

- Use serial execution as the guaranteed baseline.
- Use `parfeval`/`DataQueue` or the existing parallel capability for optional background execution when the toolbox is available.
- Keep all solver work off the UI thread where supported.
- Save completed cases incrementally so a later failure does not discard earlier setups.
- Support cancellation and mark the run `cancelled` with completed partial results retained.
- Catch per-speed and per-case failures, preserve the error message and stack in metadata, and leave a visible gap/warning in the plots.
- Reject incompatible cache files with a clear schema/version message.
- Never interpret `fmincon` exit flag 2 as automatically valid; use the recorded residual and feasibility gate already used by the existing solvers.

## Verification strategy

The implementation must add MATLAB tests for:

1. Result-schema construction and migration from the existing ramp study cache.
2. Lateral adapter equivalence to direct `rampSweep` on a small deterministic speed grid.
3. Longitudinal adapter enforcement of pure `Ay = 0`, zero steer, and zero yaw rate.
4. Preservation of named aero, balance, suspension, tire, and validity fields.
5. Overlay plotting for multiple setups, differing speed grids, invalid points, and baseline deltas.
6. Correct disabling/marking of lateral-only metrics for longitudinal runs.
7. Cache round-trip and incremental/partial-run status.
8. Unit-display conversion without changing stored SI values.

For visual QA, run one representative car at a small speed grid in both ramp types and inspect every tab. Confirm that positive/negative understeer signs, front/rear aero loads, drag/downforce, wheel-lift markers, truncated ramps, and power-limited points agree with the direct MATLAB plots.

## Acceptance criteria

The v1 app is acceptable when a MATLAB user can:

1. Select at least two different setup cases.
2. Run a lateral-limit study and see capability, mechanical/aero balance, yaw-derived understeer/oversteer, aero, suspension, tire, and validity metrics versus `vCar`.
3. Run a pure longitudinal study and verify the `Ay = 0` condition in the inspector.
4. Overlay the selected setups with clear labels and warnings for incompatible comparisons.
5. Click or select a speed and inspect the detailed named state and four-wheel loads.
6. Save the study, close/reopen the app, reload it, and reproduce the same plots.
7. Export a normalized table and figure without losing validity flags or units.
8. Identify solver failures, truncated ramps, power limitation, aero-map extrapolation, and wheel lift without reading the MATLAB command window.

## As-built canonical implementation

The implementation uses carConfigBaseline.m as the independent source of truth for Ramp Speed construction. A setup specification is serializable and contains the discrete suspension choices, numeric ride heights, driver weight, rear weight distribution, aero-map ID, and baseline/configuration version. rampSpeed.buildSetupCatalog builds exactly one solver-ready Car per setup and returns an N-by-1 catalog shared by lateral and pure-longitudinal studies.

rampSpeed.RampSpeedSession owns setup selection, duplication/editing, deletion, save/load, cancellation, and result acceptance. rampSpeed.StudyExecutor owns one study job and progress lifecycle, while rampSpeed.runStudy and rampSpeed.runCase own the numerical execution. The App Designer file is a view/controller shell over that session; it does not create solver futures or call solver internals from callbacks.

The supported function entry point is:

```matlab
[study,runs,figures,events] = runRampSpeedStudy( ...
    struct('rampType','longitudinal','speeds',[5 10 15 17.5 20 22.5 25]));
```

The compatibility wrapper delegates to the same baseline/session-independent executor path and no longer calls the legacy carConfig/rampSweep study script. Legacy study files without setup specifications remain loadable as read-only result data. The study status model preserves invalid gaps and uses explicit planned, running, converged, near_feasible, infeasible, solver_failed, and cancelled states; unavailable continuous values are stored as NaN.

The v1 metric catalog includes pure longitudinal capability and drag deceleration, lateral capability, mechanical/aero/yaw balance, lift and drag, ride height, front/rear and four-corner camber, and four-corner slip angles. The stored unit contract remains SI; app display conversion is presentation-only.

## Implementation sequencing

The implementation should be executed as independent tracks, then integrated:

- numerical result schema, metadata, cache, and longitudinal adapter
- lateral adapter and regression coverage
- plot catalog and multi-setup comparison behavior
- App Designer shell, progress, cancellation, and inspector
- packaging, visual QA, and user documentation

When implementation begins, use multiple Luna agents for independent tracks where the files and tests do not overlap, then perform an integration pass against the full MATLAB test suite before claiming completion.


