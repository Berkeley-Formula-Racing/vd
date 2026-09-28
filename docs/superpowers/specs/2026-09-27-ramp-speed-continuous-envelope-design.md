# Ramp-Speed Continuous-Envelope Powertrain Design

## Decision

The Ramp Speed solver will use `continuousEnvelope` as its only ramp-speed
powertrain mode. It will not expose an explicit-gear/CVT selector in the UI.
The existing explicit-gear model remains an offline reference path only.

## Model boundary

At each positive vehicle speed, the ramp path compiles one bounded continuous
drivetrain ratio from the real powertrain's minimum and maximum total
reductions. The ratio maximizes full-throttle wheel force while respecting the
engine torque-map domain and redline. That ratio is reused during the reduced
pure-longitudinal solve for that speed.

The full `Car` remains the vehicle-equation boundary. `Car.equations` accepts
an optional `continuousRatio` evaluation option; legacy callers continue to
use automatic or integer `gearOverride` behavior. Continuous evaluation does
not mutate the `Car`, and its result reports `current_gear = NaN` plus ramp
diagnostics containing the effective ratio and engine speed.

## Numerical behavior

The continuous envelope removes discrete gear-transition branches from the
primary ramp solver. Tire, aero, rolling resistance, rotating mass, and engine
torque remain active. Invalid ratio bounds or nonfinite torque data return an
explicit solver failure/infeasible result; they never become a plotted zero.

## Rollout

1. Add the envelope builder and continuous-ratio Car evaluation seam.
2. Replace the explicit-gear loop in `solveLongitudinalPoint` with one
   envelope-backed solve.
3. Record `powertrainModel`, effective ratio, engine RPM, and envelope
   provenance in the canonical result.
4. Add a compact cached ramp-model layer around the same data boundary after
   the continuous path matches the reference model.

## Verification

The focused tests must cover ratio bounds, redline, invalid speed handling,
Car immutability, pure-longitudinal solver status, and comparison with the
legacy explicit-gear reference at representative speeds including 20 m/s.
