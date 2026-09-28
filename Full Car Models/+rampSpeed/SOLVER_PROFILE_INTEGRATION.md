# Ramp-Speed Solver Profiles and Fixed Speed Presets

Both settings pass through the existing simulation entry-point settings struct,
directly or as request.settings to rampSpeed.runStudy.

## Solver profiles

Set a named profile and optional validated overrides:

~~~matlab
request.settings.solverProfile = "fastPreview";
request.settings.solverOptions.maxFunctionEvaluations = 800;
~~~

Supported IDs are accurate (default), fastPreview, and
approximateAeroPreview. runStudy records the resolved profile at study level;
runLongitudinalRamp and runLateralRamp record it per run and apply it to a
local value-copy of Car. max_long_accel consumes the profile optimizer
options. The lateral solver retains its existing optimizer settings.

## Pure-longitudinal speed presets

The Ramp Speed runner uses one deterministic fixed grid per run. Select one of
the three presets:

~~~matlab
[speeds,meta] = rampSpeed.fixedSpeedGrid([5 30],"preview");      % 5 points
[speeds,meta] = rampSpeed.fixedSpeedGrid([5 30],"accurate");     % 11 points
[speeds,meta] = rampSpeed.fixedSpeedGrid([5 30],"highAccuracy"); % 21 points
~~~

The first and last values define the speed range. Every preset includes both
endpoints and uses uniform spacing, so a high-accuracy run always evaluates
all 21 planned speeds rather than stopping after a seed scan. The resolved
grid and per-speed provenance are recorded in runMeta.speedGrid.

## Aero hot paths

Aero.coefficientsNumeric and AeroMap.evaluateNumeric provide batched numeric
coefficient evaluation without per-query structs. Complete planar grids use
gridded interpolation; other maps retain scattered interpolation.
benchmark_aeroMapEvaluation measures the batch and scalar paths.

Car.solveRideHeightAero accepts an explicit context scoped to vehicle and
ride-height settings, AeroMap identity/corrections, and operating point.
max_long_accel carries that context between objective and constraint
evaluations only when coupled ride-height aero is active. It is only a
starting guess: the coupled residual is still solved to the original
tolerance.
