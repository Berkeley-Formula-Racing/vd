# Ramp-Speed Rewrite Architecture

## Goal

Make the Ramp Speed application numerically trustworthy and maintainable while preserving the current MATLAB App Designer layout, setup editor, setup overlays, and MATLAB-only sharing model.

## Required behavior

- MATLAB remains the implementation platform.
- One Ramp Speed setup owns exactly one solver-ready `Car` and can run in either lateral or pure-longitudinal mode.
- Pure longitudinal mode means `Ay = 0`, zero steer, zero vehicle lateral velocity, zero yaw rate, zero front slip, and symmetric rear slip unless a diagnostic explicitly proves that assumption invalid.
- Every planned speed produces exactly one result row with an explicit status: `planned`, `running`, `converged`, `near_feasible`, `infeasible`, `solver_failed`, or `cancelled`.
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

## Architecture

The App Designer class becomes a thin view over a testable `RampSpeedSession` controller. A single executor owns run lifecycle and delegates to a study coordinator; the coordinator delegates each speed task to either a lateral or longitudinal point solver. Solvers return typed speed results, the coordinator assembles one canonical run, and analysis/persistence consume that same result without reaching back into solver internals.

The current `Car` equations remain the vehicle-model boundary. The new pure-longitudinal Ramp Speed solver reduces the optimizer to the independent throttle/rear-slip variables, enumerates gear branches explicitly, and retains the existing solver as a reference/fallback during migration. The legacy `runRampSpeedStudy` path becomes a compatibility adapter and is not allowed to create a second canonical result contract.

## Compatibility

The current UI layout and setup fields are retained. Existing schema-v1 studies can be loaded as read-only data. Existing solver functions remain callable by other vehicle-model consumers while the Ramp Speed app migrates to the new coordinator and point-solver interfaces.

## Verification

Acceptance requires focused MATLAB tests for request validation, gear transitions, reduced-versus-reference longitudinal results, lateral compatibility, adaptive no-duplicate scheduling, serial execution without parallel APIs, cooperative cancellation, status-preserving checkpoint/load, setup/asset fingerprints, correct raw ramp axes, longitudinal capability metrics, overlays, and a real App-shaped end-to-end run around 15–30 m/s.
