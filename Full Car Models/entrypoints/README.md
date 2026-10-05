# Full Car Models entrypoints

## Build or edit a QSS car

`carConfig()` keeps returning the legacy N-by-2 cell array: lap cars are in
column 1 and acceleration cars are in column 2. The editable defaults now come
from `vehicleConfig.defaultSetup`; make a local setup edit and rebuild it with
`vehicleConfig.buildCar` instead of changing values inside `carConfig.m`:

```matlab
setup = vehicleConfig.defaultSetup("legacy","nominal");
setup.parameters.mass = 160; % kg
[lapCar,derived] = vehicleConfig.buildCar(setup,"lap");
```

Use `carConfig("Explicit",overrides,"legacy","nominal","cases")` when a
table of parameter overrides should produce named cases. The named form
includes each setup plus its lap and acceleration cars. See
[README_car_configuration.md](../README_car_configuration.md) for the setup
fields, units, value-class behavior, and MAT-file guidance. Ramp Speed uses the
independent frozen `carConfigBaseline()` configuration.
The runnable MATLAB scripts live in this folder. Each script bootstraps the
full model path and changes MATLAB's working folder to the repository root, so
it can be launched from any current folder.

To run one scalar baseline QSS autocross simulation and export telemetry:

```matlab
run('C:\VD\.worktrees\telemetry-viewer\Full Car Models\entrypoints\run_single_baseline_qss.m')
```

The output is `baseline.qss.h5` in the repository/worktree root (currently
`C:\VD\.worktrees\telemetry-viewer\baseline.qss.h5`), or a numbered sibling if
that file already exists. To open it in the Python viewer:

```powershell
cd C:\VD\.worktrees\telemetry-viewer\telemetry_viewer
uv run qss-telemetry-viewer C:\VD\.worktrees\telemetry-viewer\baseline.qss.h5
```
