# Ramp Speed App

RampSpeedApp is the schema-v1 MATLAB App Designer front end for lateral and longitudinal ramp studies. The project root is this `Full Car Models` folder; `RampSpeedApp.prj` uses `.` to resolve paths relative to its containing folder.

## Open and initialize

1. Open `RampSpeedApp.prj` from MATLAB.
2. Use the project startup shortcut `setup_paths` (the project startup file is `setup_paths.m`). It adds the app, `+rampSpeed`, tests, and model folders to the MATLAB path. Keep the generated `resources/project` sidecar directory with `RampSpeedApp.prj`; MATLAB uses it for the project file inventory, name, startup action, and shortcut.
3. Open `RampSpeedApp.mlapp` and press Run. The app opens on the analysis view with a `Setups / Run` sidebar. Use its checkboxes to choose the setups shown in the plots; the shared run configuration, units, Run/Cancel/Clear, save/load/export, and progress controls are in the same sidebar. The `Hide sidebar` button collapses the sidebar when more plot width is useful.

Ramp Speed uses one protected, solver-ready baseline setup for both ramp types. The baseline is built by carConfigBaseline.m, not by the frequently changing carConfig.m. Its physical datum is the axle-plane extrapolation of the CAD ride-height line: 3.865589 in front and 5.769078 in rear. The selected b26 aero map is referenced at those same axle heights, so the baseline has zero ride-height map offset. Open the `Setup` tab to duplicate the baseline or another setup, then edit rear ARB stiffness, front/rear spring selections, front/rear ride height (always entered in inches), driver weight, rear weight distribution, and aero-map ID. The ramp type selects only the solver; it does not select a second car. Each setup is represented by one Car in the N-by-1 setup catalog. The `Setup` tab is intentionally setup-only; run configuration is shared and remains in the sidebar, while the compact overlay legend provides selection control during analysis. The plot axes reset to automatic limits on each redraw.

## Study setup and interpretation

- The speed-grid selector has three deterministic fixed presets: Preview (5 points), Accurate (11 points), and High accuracy (21 points). The start/stop range supplies the endpoints for the selected preset; an explicit speed vector remains available for custom fixed studies.

- Lateral runs support `coast` and `balanced` modes. Coast holds the coast-style longitudinal condition; balanced applies the requested balance condition while the ramp is solved.
- Longitudinal runs use the pure-`Ay=0` convention. The exported `pure_ay0`, `steer_zero`, `lat_velocity_zero`, and `yaw_rate_zero` fields make those constraints explicit instead of treating unavailable lateral metrics as zero.
- Longitudinal powertrain evaluation uses the continuous-envelope model. Each setup builds one serializable `rampSpeedLite` snapshot and caches the bounded ratio envelope for the fixed speed plan; the existing full `Car` equations remain the tire/aero force boundary while the lightweight model is expanded.
- The default display is SI. A unit selection may change plot labels and displayed values, but the study and exports remain canonical SI: speed in `m/s`, acceleration in `m/s^2`, force in `N`, length in `m`, angle in `rad`, and angular rate in `rad/s`.
- Understeer is signed. A positive `K_linear_rad_per_mps2` means steering demand grows with lateral acceleration (understeer); a negative value is oversteer. Preserve the sign when comparing cases or exporting deltas.
- `downforce_N` and `drag_N` are separate aerodynamic quantities. Downforce increases normal load; drag is the longitudinal resistance. Use `ClA_m2`, `CdA_m2`, `LoD`, axle loads, and aero-balance fields together when checking an aero result.
- Corner tire outputs are split into front-left/front-right camber and rear-left/rear-right camber plots. The same four-corner convention is used for tire slip angles, so a left/right balance change is not hidden by axle averaging.

Every per-speed row and detailed point carries `valid` and `status`; per-speed rows also carry `reason`. Do not replace invalid rows with interpolated values. Review the warning/status fields before comparing curves, especially:

- `truncated` and `ramp_complete_fraction` for a ramp that did not reach its requested limit;
- `power_limited`, `traction_limited`, and `rear_slip_upper_active` for longitudinal or power-constrained points;
- `wheel_lift` and the minimum normal-load fields for unloaded or lifted tires;
- `aero_outside_map` and `aero_residual_m` for points outside the aero map or with a map residual.

The Raw Ramp and Inspector/Data tabs retain the detailed solver rows. Balance and Aero & Loads plots should be read with validity gaps and warning markers intact.

The canonical programmatic entry point is the function form of runRampSpeedStudy; the App uses the same session/executor lifecycle:

```matlab
opts = struct('rampType','longitudinal', ...
    'speeds',[5 10 15 17.5 20 22.5 25], ...
    'makeFigures',true);
[study,runs,figures,events] = runRampSpeedStudy(opts);
```

For application integrations, use rampSpeed.RampSpeedSession and rampSpeed.StudyExecutor. rampSpeed.runStudy is the canonical numerical runner. Serial execution is deterministic and supported everywhere; parallel execution is optional and is one study job over the selected setup catalog. Per-speed statuses use planned, running, converged, near_feasible, infeasible, solver_failed, or cancelled; invalid continuous metrics are NaN, never numeric zero.

## Save, cache, and export

Use the app Save/Load controls for schema-v1 `.mat` studies. Saved setup specifications include the setup fields, aero-map ID, driver weight/distribution, and baseline/configuration version; raw Car objects are rebuilt rather than used as the editable source of truth. Existing studies without setup specifications remain loadable as read-only result data. A runner checkpoint is written only when the request supplies a checkpoint path; it is not a hidden repository cache. Keep temporary checkpoints in a disposable working folder and retain released studies in a named results folder.

Longitudinal run metadata includes the data-only ramp model, its solver-profile
version, setup identity, aero-map provenance, and cached envelope records. A
speed with no feasible continuous ratio remains an explicit infeasible or
solver-failed row with a reason; it is never converted to a plotted zero.

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

The focused Ramp Speed acceptance check also exercises one real baseline setup through both ramp types at representative speeds:

```matlab
root = pwd;
addpath(genpath(root));
realCarResults = runtests(fullfile(root,'tests','test_rampSpeedEndToEndRealCar.m'));
assert(all([realCarResults.Passed]));
```
