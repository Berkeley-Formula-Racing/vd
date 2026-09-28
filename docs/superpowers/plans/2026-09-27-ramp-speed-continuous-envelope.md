# Ramp-Speed Continuous-Envelope Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task.

**Goal:** Replace the Ramp Speed primary longitudinal powertrain path with one cached continuous-ratio envelope while preserving the existing full-Car reference behavior.

**Architecture:** Compile one bounded ratio and torque/wheel-force summary per positive speed, pass that ratio through an optional `Car` evaluation seam, and solve throttle/rear slip once instead of enumerating explicit gears. Keep setup, aero, tire, result, and App contracts unchanged except for powertrain provenance.

**Tech Stack:** MATLAB `classdef`, MATLAB Unit Test, existing `Car`/`Powertrain`, `fmincon`, and the `+rampSpeed` canonical study path.

**Spec:** `docs/superpowers/specs/2026-09-27-ramp-speed-continuous-envelope-design.md`

## Global Constraints

- `continuousEnvelope` is the only Ramp Speed powertrain mode exposed to the app.
- Existing explicit-gear behavior remains callable by legacy vehicle-model consumers and is used only as an offline reference during validation.
- Continuous metrics that are unavailable remain `NaN`; no failed speed may be represented by zero.
- A setup still owns one solver-ready `Car` and may run in either ramp type.
- The continuous-ratio evaluation option must not mutate the source `Car`.
- Serial execution remains independent of Parallel Computing Toolbox APIs.

## Review Focus

- Torque-map endpoints and redline: a speed whose allowable ratio interval clips at either endpoint must remain finite and explicit.
- Low-speed ratio bounds: the builder must reject an empty feasible ratio interval instead of extrapolating engine torque.
- 20 m/s transition behavior: the continuous result must remain valid or report a reasoned status without a zero substitution.
- Legacy callers: integer gear override and automatic gear selection must retain their existing behavior.
- Setup overlays: powertrain provenance must remain distinct in canonical metadata without changing setup identity.

### Task 1: Continuous envelope and Car evaluation seam

**Files:**
- Create: `Full Car Models/+rampSpeed/buildContinuousEnvelope.m`
- Modify: `Full Car Models/carComponents/Powertrain.m`
- Modify: `Full Car Models/carComponents/Car.m`
- Test: `Full Car Models/tests/test_rampSpeedContinuousEnvelope.m`

**Interfaces:**
- `envelope = rampSpeed.buildContinuousEnvelope(car,speed_mps)` returns a scalar struct with `powertrainModel`, `speed_mps`, `drivetrainReduction`, `ratioLowerBound`, `ratioUpperBound`, `engineRpm`, `engineTorque_Nm`, `wheelForce_N`, and `exitflag`.
- `car.equations(P,rideHeightContext,struct("continuousRatio",ratio))` uses that total reduction and returns `current_gear = NaN`; existing calls remain unchanged.

- [ ] Write and run the failing continuous-envelope tests.
- [ ] Implement bounded ratio search over the torque-map domain and redline.
- [ ] Implement optional reduction handling in `Powertrain.wheel_torques` and `Car.rotatingMass`/`Car.equations` without changing legacy call signatures.
- [ ] Run the focused tests and commit the slice.

### Task 2: Switch the reduced ramp solver and profile metadata

**Files:**
- Modify: `Full Car Models/+rampSpeed/solveLongitudinalPoint.m`
- Modify: `Full Car Models/+rampSpeed/solverProfiles.m`
- Modify: `Full Car Models/+rampSpeed/validateSolverProfile.m`
- Modify: `Full Car Models/+rampSpeed/serializeSolverProfile.m`
- Modify: `Full Car Models/tests/test_rampSpeedLongitudinalPointSolver.m`
- Test: `Full Car Models/tests/test_rampSpeedContinuousEnvelopeSolver.m`

**Interfaces:**
- Every resolved profile contains `powertrainModel = "continuousEnvelope"`.
- `solveLongitudinalPoint` performs one envelope-backed reduced solve and stores the envelope in diagnostics; `gearAttempts` is empty for the primary path.

- [ ] Add failing solver assertions for continuous-envelope provenance and no explicit-gear primary loop.
- [ ] Replace the gear loop with one envelope-backed solve while retaining cancellation and reference diagnostics.
- [ ] Serialize the powertrain model in saved run metadata.
- [ ] Run focused solver, adapter, profile, and reference-comparison tests.
- [ ] Commit the solver slice.

### Task 3: Add the first compact ramp-model cache

**Files:**
- Create: `Full Car Models/+rampSpeed/buildRampModel.m`
- Create: `Full Car Models/+rampSpeed/compilePowertrainEnvelope.m`
- Create: `Full Car Models/tests/test_rampSpeedCompactModel.m`
- Modify: `Full Car Models/+rampSpeed/runCanonicalLongitudinalRamp.m`

**Interfaces:**
- `model = rampSpeed.buildRampModel(car,setupSpec,profile)` is a data-only, serializable ramp model with the setup identity, baseline version, mass/geometry, driver-adjusted values, aero map provenance, and a cached continuous-envelope evaluator.
- The first integration may continue delegating tire/aero force evaluation to the existing `Car`; it must build the compact model once per setup, not once per speed.

- [ ] Add failing tests for one model per setup and reuse across a fixed speed grid.
- [ ] Implement the data-only model and cached powertrain envelope.
- [ ] Pass model provenance through the canonical run metadata.
- [ ] Run compact-model and end-to-end ramp tests.
- [ ] Commit the cache slice.

### Task 4: Regression, App smoke, and performance verification

**Files:**
- Modify: `Full Car Models/tests/test_rampSpeedEndToEndRealCar.m`
- Modify: `Full Car Models/tests/test_rampSpeedAppSmoke.m`
- Modify: `Full Car Models/README_ramp_speed_app.md`

- [ ] Cover 5, 10, 15, 17.5, 20, 22.5, and 25 m/s in serial longitudinal mode.
- [ ] Verify status/reason fields and no invalid zero substitution.
- [ ] Verify setup overlays retain distinct setup and powertrain provenance.
- [ ] Run the focused batch, App smoke suite, broader Ramp Speed regression suite, and `git diff --check`.
- [ ] Record measured runtime versus the prior solver without asserting machine-specific absolute timing.
