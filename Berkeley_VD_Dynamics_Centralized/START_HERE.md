# Shared vehicle and tire inputs — ready for a review-branch commit

Prepared against Berkeley-Formula-Racing/vd master
`144b3548425d05988668fcc93e3b2b5bd0830676` (October 2, 2026), checked October 15.

This package centralizes existing inputs without substituting SN6 data or
changing tire fits. It includes 31 modified source files and 12 new configuration,
documentation and test files. The recent steering-limit, CAD ride-height,
camber caching and solver work is preserved. The October 2 docs deletion is
preserved. The earlier unrelated Car.c property addition is excluded.

## What is included

- `shared-inputs.patch`: complete plain Git diff, including the test fixture.
- `revised_files/`: the same changed/new files, with repository-relative paths.
- `CHANGELOG.md`: changes, scope and suggested commit/PR text.
- `validation/`: actual check results, numerical summary and file hashes.

Use the patch on the existing repository. This small package is not a replacement
for the repository and deliberately does not include its large datasets.

## Apply and commit

Start with a clean checkout of your fork or a local copy of the original repo.
Do not overwrite uncommitted work. Fetch the upstream master and create a review
branch. If `origin` is your fork, use the remote pointing to the club's repo in
the first two commands below instead.

```bash
git fetch origin master
git switch -c centralize-vd-inputs origin/master
git rev-parse HEAD
```

The recorded base is `144b3548425d05988668fcc93e3b2b5bd0830676`. If it has advanced,
review the newer changes before applying; do not force an old patch over them.
Set the patch path below to where you extracted this package:

```bash
git apply --check /absolute/path/to/shared-inputs.patch
git apply --index /absolute/path/to/shared-inputs.patch
git diff --cached --check
git diff --cached --stat
git diff --cached
git commit -m "Centralize existing vehicle and tire inputs"
```

The patch stages only its own 43 files. There is no need for `git add .` or to
copy the package's top-level notes into the repository. Commit using your own
Git identity. Nothing has been pushed or published by this work.

## What passed

- Clean patch application against the current upstream file versions; every
  resulting file matches `revised_files/` exactly.
- MISS_HIT source checks for all 40 changed/new MATLAB source files.
- GNU Octave 8.4.0 regression checks: complete main/ramp defaults, 185 historical
  inputs, temporary master-edit propagation and coefficient-loading validation.
- 375 tire states / 750 force outputs match upstream exactly (maximum difference 0).
- Three suspension Car presets and 7-DOF matrices; vibration inputs checked
  separately on a struct; nine Car equation states and a 101-sample response.
- The stored baseline fixture was independently reproduced from current master.

## Remaining validation before merging into shared master

Run the supplied checks in MATLAB, then the club's usual MATLAB/app checks:

```matlab
addpath(fullfile(pwd,'tests','shared_inputs'));
results = runSharedInputChecks();
```

This MATLAB branch builds both production factories and exercises the table/aero
map path. It has not been run here because MATLAB is unavailable. Apps, Simulink
and full lap optimization have not been run. The Octave run checks the numerical
core with one isolated interpolation adapter; it is not a MATLAB certification.

Several older scripts already have missing inputs or errors, including the
missing `Fy_pure_parameters_run24_17.mat`, undefined `delta`/`zeta` uses, and the
vibration script's undeclared `Car.c`. They are documented in the configuration
README and intentionally not repaired in an input-location refactor. The package
is ready to commit for review; it does not claim every old script is runnable or
that the original engineering equations have been audited.

## Where future edits go

- `Full Car Models/config/vdVehicleConfig.m`
- `Full Car Models/config/vdTireConfig.m`

Historical studies retain named presets in those files so their different
calibrations and native units are not silently overwritten. See
`Full Car Models/config/README.md` for the mapping and editing conventions.
