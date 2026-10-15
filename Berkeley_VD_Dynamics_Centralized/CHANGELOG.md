# Changes for review

Suggested commit title: **Centralize existing vehicle and tire inputs**

Vehicle and tire defaults were repeated in the full-car and ramp builders and
in standalone studies. This change moves the editable inputs into shared
configuration files while retaining the existing numbers, unit conversions,
fit selections and calculation code. Historical setups remain named presets.

The patch is based on current master `144b3548425d05988668fcc93e3b2b5bd0830676`.
The October 2 commit only deleted documentation; its cleanup remains intact.
The October 1 steering stop, CAD ride-height conversion, camber cache, build
helpers and solver changes are retained.

## Production changes

| Files / group | Exact change |
| --- | --- |
| `Full Car Models/config/vdVehicleConfig.m` | Moves full-car, aero, engine, drivetrain, brake, event and ramp option defaults into one file. Also contains the original suspension and standalone-study presets. The engine handle is evaluated by the builders. |
| `Full Car Models/config/vdTireConfig.m` | Centralizes default pressure, camber, grip factors, coefficient/camber paths, suspension tire stiffness and historical fit/stiffness presets. Supports settings-only reads without loading full-car coefficients. |
| `Full Car Models/config/vdRampConfig.m` | Adapts the shared values to the existing ramp configuration schema; preserves baseline version 3 and the CAD-to-axle height calculation. |
| `Full Car Models/config/vdLoadTireParameters.m` | Loads Xbestcell into an isolated structure and checks the expected 26 Fx / 42 Fy finite real coefficients. |
| `Full Car Models/config/vdApplySuspensionConfig.m` | Applies the existing named suspension settings and damping-curve selection. |
| `carConfig.m`, `carConfigBaseline.m` | Replace duplicated default blocks with configuration reads. Preserve signatures, sampling/override flow, output structures and model-selection string conversion. |
| `utilities/parameters_loop.m`, `+rampSpeed/loadRampAssets.m` | Obtain camber file paths from the tire master while preserving the current per-build cache. |
| `carComponents/Tire2.m`, `carComponents/Camber_Evaluation.m` | Route fallback camber asset reads through the shared paths. Force and camber equations are unchanged. |
| `entrypoints/ModalAnalysis.m`, `ModalAnalysisDecoupled.m`, `transientSim.m`; `Suspension Transient Model/FullCarVibrationModel.m` | Read the original suspension presets. No change to dynamics or modal calculations; the upstream vibration Car.c issue remains. |
| Nine `Bicycle Model/` scripts listed below | Replace local geometry, mass, inertia, steering ratio, stiffness and fit references with named presets where applicable. Preserve the original tire-load inputs, including the distinct toe-study load mass. |
| Eight `CamberCurveGen/` and `CamberCurveGen2/` scripts listed below | Move their original vehicle/tire input blocks into named presets and use explicit fit paths. Their optimization calculations and unit conventions are preserved. |
| `roll.m` | Reads its existing mass, roll gradient, ride rates, aero values and density from its preset. |
| `understeer_bode.m` | Reads existing geometry, rear mass fraction, mass, inertia, tire stiffness and target understeer input. The state-space output matrix remains in the script. |
| `understeer_budget.m` | Moves both front and rear inputs, including supplied force/moment values, to presets without correcting their existing values or unit comments. |
| `calcWheelForcesAndDisplacements_pitch.m` | Reads the four existing anti-lift/dive/squat coefficients from its preset. |
| `config/README.md`, `tests/shared_inputs/` | Editing guidance, default fixture and regression checks with an isolated Octave interpolation adapter. |

Bicycle scripts: `BicycleModel.m`, `LatAccelResponseTimePlotter.m`,
`YawResponseTimePlotter.m`, `toe_compliance.m`, `ODE Formulation/BicycleModelODE.m`,
and `Simulink Models/{BicycleModelSim,BicycleModelSimv2,FullCarModel,NonLin_BicycleModel}.m`.

Camber scripts: `Testcamber.m`, `test.m`, `generateIdealCamberCurves.m` and
`generateIdealCamberCurves_simplified.m` in each of `Full Car Models/CamberCurveGen/`
and `Full Car Models/CamberCurveGen2/`.

## Changes from the earlier package

- Updated the recorded base to the latest master.
- Excluded the unrelated Car.c property addition and `docs/SHARED_INPUTS.md`.
  Editing notes now live beside the configuration and in this package.
- Restored unrelated whitespace to keep the diff focused.
- Finished the partially moved front/rear and embedded study inputs.
- Kept the understeer state-space output matrix out of the tire configuration.
- Preserved MATLAB string conversion at the existing factory boundary.
- Removed unnecessary full-car coefficient loads from settings-only consumers.
- Re-ran the numerical checks and independently rechecked the reference fixture
  against current master. Verified clean patch application and output hashes.

## Validation and scope

All 40 changed/new MATLAB source files pass MISS_HIT. GNU Octave 8.4.0 passes
checks for original-default equivalence, 185 historical inputs, edit propagation,
375 tire states, suspension inputs, nine vehicle-equation states and a two-second
response integration. Fx/Fy match upstream exactly. See `validation/RESULTS.json`
and the actual logs.

MATLAB production factories, app/UI paths, Simulink and full lap optimization
remain unrun. Existing missing fits and script errors are documented in the
configuration README. These tests check the input refactor, not the accuracy
of the original physical model. No SN6, TUM, original MAT data, solver changes,
Car-class changes, archived rewrites or documentation resurrection are included.
