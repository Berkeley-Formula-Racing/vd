# Ramp Speed App

RampSpeedApp is the schema-v1 MATLAB App Designer front end for lateral and longitudinal ramp studies. The project root is this `Full Car Models` folder; `RampSpeedApp.prj` uses `.` to resolve paths relative to its containing folder.

## Open and initialize

1. Open `RampSpeedApp.prj` from MATLAB.
2. Use the project startup shortcut `setup_paths` (the project startup file is `setup_paths.m`). It adds the app, `+rampSpeed`, tests, and model folders to the MATLAB path.
3. Open `RampSpeedApp.mlapp` and press Run. The Setup table supplies the available car/case definitions.

Choose the car role deliberately. `lap` is the normal lateral/cornering role; `acceleration` is the longitudinal role. The role is stored with the case and is part of the study provenance.

## Study setup and interpretation

- Lateral runs support `coast` and `balanced` modes. Coast holds the coast-style longitudinal condition; balanced applies the requested balance condition while the ramp is solved.
- Longitudinal runs use the pure-`Ay=0` convention. The exported `pure_ay0`, `steer_zero`, `lat_velocity_zero`, and `yaw_rate_zero` fields make those constraints explicit instead of treating unavailable lateral metrics as zero.
- The default display is SI. A unit selection may change plot labels and displayed values, but the study and exports remain canonical SI: speed in `m/s`, acceleration in `m/s^2`, force in `N`, length in `m`, angle in `rad`, and angular rate in `rad/s`.
- Understeer is signed. A positive `K_linear_rad_per_mps2` means steering demand grows with lateral acceleration (understeer); a negative value is oversteer. Preserve the sign when comparing cases or exporting deltas.
- `downforce_N` and `drag_N` are separate aerodynamic quantities. Downforce increases normal load; drag is the longitudinal resistance. Use `ClA_m2`, `CdA_m2`, `LoD`, axle loads, and aero-balance fields together when checking an aero result.

Every per-speed row and detailed point carries `valid` and `status`; per-speed rows also carry `reason`. Do not replace invalid rows with interpolated values. Review the warning/status fields before comparing curves, especially:

- `truncated` and `ramp_complete_fraction` for a ramp that did not reach its requested limit;
- `power_limited`, `traction_limited`, and `rear_slip_upper_active` for longitudinal or power-constrained points;
- `wheel_lift` and the minimum normal-load fields for unloaded or lifted tires;
- `aero_outside_map` and `aero_residual_m` for points outside the aero map or with a map residual.

The Raw Ramp and Inspector/Data tabs retain the detailed solver rows. Balance and Aero & Loads plots should be read with validity gaps and warning markers intact.

## Save, cache, and export

Use the app Save/Load controls for schema-v1 `.mat` studies. A runner checkpoint is written only when the request supplies a checkpoint path; it is not a hidden repository cache. Keep temporary checkpoints in a disposable working folder and retain released studies in a named results folder.

`rampSpeed.exportStudy(study,outputDirectory,options)` validates the current schema and writes absolute-path results:

- `<base>_per_speed.csv` with run context and canonical per-speed columns;
- `<base>_points.csv` with run context, detailed points, and point reasons;
- `<base>_study.mat`, loadable by `rampSpeed.loadStudy`;
- `<base>_metadata.json` with schema/app information and canonical units;
- `<base>_terminal_log.csv` plus terminal rows for complete, failed, and cancelled runs;
- optional PNGs for requested figure handles using `options.visibleFigures` (or `options.figures`).

Example:

```matlab
root = fileparts(which('RampSpeedApp'));
study = rampSpeed.loadStudy(fullfile(root,'results','study.mat'),'RampSpeedApp');
options = struct('baseName','cornering_release', ...
    'visibleFigures',struct('figure',gcf,'name','capability'));
exported = rampSpeed.exportStudy(study,fullfile(root,'exports'),options);
```

For a release record, call `manifest = rampSpeed.releaseManifest(root)`. It records the MATLAB release, schema/app version, project/package/test/model inventory, Git commit and dirty state when available, toolbox availability, and creation time. It does not require optional toolboxes or a working Git checkout to return a manifest.

## MATLAB and toolbox requirements

Use MATLAB with App Designer support. The package uses standard MATLAB tables, `jsonencode`, and `exportgraphics` for requested figure files. The solver path may use Optimization Toolbox and project-specific vehicle-model dependencies. `releaseManifest` reports optional license/toolbox availability rather than failing when one is absent.

If parallel execution is unavailable, select the serial execution mode and set workers to one. Serial mode is the supported fallback for machines without Parallel Computing Toolbox or for deterministic debugging. A missing Optimization Toolbox prevents the optimization-based solver path; it is not silently replaced by a different physical model.

## Focused and full verification

From the project root, run the focused export/manifest suites:

```matlab
root = pwd;
addpath(genpath(root));
exportResults = runtests(fullfile(root,'tests','test_rampSpeedExport.m'));
manifestResults = runtests(fullfile(root,'tests','test_rampSpeedReleaseManifest.m'));
assert(all([exportResults.Passed]));
assert(all([manifestResults.Passed]));
```

Run the full repository test folder with:

```matlab
root = pwd;
addpath(genpath(root));
allResults = runtests(fullfile(root,'tests'));
```

The full folder includes legacy model/plot tests in addition to the app contract tests. Class-shadowing warnings from duplicate historical model folders and any pre-existing legacy plot baseline failures should be recorded with the test result rather than “fixed” by changing this app package.
