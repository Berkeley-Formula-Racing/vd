# Shared inputs

Edit `vdVehicleConfig.m` for vehicle, aero, drivetrain, event and suspension
defaults. Edit `vdTireConfig.m` for tire settings and coefficient-file paths.
Both full-car builders now read these files. Their existing call signatures,
DOE overrides and optional tire-model selection are retained.

This is an input-location refactor. Existing numeric values, unit conversions,
tire coefficient files and calculation equations are retained.

## Files

| File | Purpose |
| --- | --- |
| `vdVehicleConfig.m` | Shared full-car defaults; ramp options; named suspension and historical study presets. |
| `vdTireConfig.m` | Default tire settings and asset paths; historical fit choices and stiffness inputs. |
| `vdRampConfig.m` | Adapts shared defaults to the existing ramp setup structure, including its CAD ride-height conversion. |
| `vdApplySuspensionConfig.m` | Applies a study's existing suspension settings. |
| `vdLoadTireParameters.m` | Loads and checks the full-car fit coefficients. |

The full-car defaults remain 162 kg without driver, 62-inch wheelbase,
47-inch track, 11.75-inch CG height and 12 psi tire pressure. `weight_dist`
is the **rear** mass fraction. Geometry in `cfg.car` is in metres; the old
study presets retain the units expected by their source scripts.

The default coefficient files remain:

- `carComponents/Fx_combined_parameters_run38_30.mat`
- `carComponents/Lapsim_Fy_combined_parameters_1965run15.mat`
- `carComponents/camberratiossmoothed.mat`
- `carComponents/camber_models_fast.mat`

Paths resolve relative to the repository. Use `vdTireConfig([],[],false)`
when only settings or paths are needed; this avoids reading unrelated
full-car coefficients. The normal call loads the coefficients as before.
The engine torque function is stored as a handle and evaluated by the two
builders, so geometry-only studies need no engine dependency.

## Historical studies

`cfg.studies` and `cfg.suspension` preserve the different setups already used
by the bicycle, camber, roll, understeer and suspension scripts. These are
not all the same car or the same tire-fit format. Edit the appropriate preset
in the shared files; do not assume editing `cfg.car.mass` changes every old
study's mass. Fields explicitly referencing another shared field follow it.
The toe-compliance study's 520 lb vehicle mass and 535 lb cornering-load input
are both retained; resolving that mismatch belongs in a separate audit.

Full-car consumers already calling `carConfig` inherit its defaults. Functions
accepting a car or parameter structure continue to use their caller's inputs.
Archived/obsolete implementations, tire-fitting tools, measurement processing,
saved car/model data and class compatibility defaults are outside this change.
No SN6 substitution or TUM integration is included.

## Checks and existing limitations

See `tests/shared_inputs/README.md` at the repository root for the regression
command and tested coverage. GNU Octave checks cover input equivalence, tire
forces, suspension inputs and the vehicle numerical core. Production MATLAB
table/aeromap entrypoints, apps and Simulink still need a MATLAB run.

Existing issues remain separate: two historical bicycle scripts refer to the
missing `Fy_pure_parameters_run24_17.mat`; `BicycleModelODE.m` expects `delta`
from its caller; the decoupled modal plot refers to undefined `zeta`; and
`FullCarVibrationModel.m` assigns `car.c`, which the current `Car` class does
not declare. This refactor does not repair those scripts or establish the
engineering accuracy of their original equations.
