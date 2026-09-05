# LC0 Nondimensional Tire Prototype Design

## Goal

Create an untracked MATLAB prototype for the Hoosier 16×7.5×10 LC0 on an
8-inch rim that preserves its measured TTC free-rolling lateral behavior and
adds explicitly uncertain longitudinal and combined-slip behavior from a
separately supplied donor dataset.

## Evidence and scope

The target data is Round 8 runs 15, 16, 18, and 19 in
`Magic Formula/TTC Documentation/Round 8 (16 in Hoosiers)`. The supplied
Round 8 run table identifies these as the 16×7.5×10 LC0 on an 8-inch rim:
12 psi and speed runs (15/18) plus 8 psi and speed runs (16/19).
Their measured longitudinal slip is effectively zero, so no target-tire
longitudinal or combined-slip parameter may be fitted from them.

The prototype must not modify `Tire2`, `carConfig`, existing TTC parsers, or
any tracked file. It must not claim a fitted full MF combined-slip model.

## Architecture

All prototype code lives in `Magic Formula/experimental/lc0_nd_tire/`.

- `lc0NDConfig` centralizes target/donor paths, run metadata, reference load,
  and calibration priors.
- `lc0NDLoadFreeRolling` loads and validates target TTC runs into a common
  SI table-like struct.
- `lc0NDLateralSummary` extracts load-bin lateral peak, stiffness, and
  pressure summaries from the target free-rolling data.
- `lc0NDCombinedForce` implements a bounded anisotropic normalized-force
  coupling law. It accepts externally supplied pure Fx and Fy functions and
  explicit `rhoMu`, `rhoStiff`, and coupling-shape calibration parameters.
- `run_lc0_nd_prototype` writes a reproducible untracked MAT result and
  diagnostic figures. It does not connect to the car model.

## Data flow

```text
Round 8 LC0 TTC runs → validated SI samples → lateral summary
                                         ↓
future donor pure/combined data → normalized calibration priors
                                         ↓
external pure Fx/Fy functions → bounded combined-force law
```

## Acceptance criteria

1. The target loader accepts only configured free-rolling LC0 runs and
   confirms negligible longitudinal slip.
2. The lateral summary reports pressure, camber, load, peak Fy/Fz, and
   small-slip stiffness without silently changing TTC signs/units.
3. The combined-force law exactly returns supplied pure Fx/Fy when the other
   axis has zero demand and never exceeds its anisotropic force limit.
4. All code, tests, generated results, and figures are untracked; no existing
   vehicle-model file is edited.
