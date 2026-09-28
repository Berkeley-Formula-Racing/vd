# Merge All Worktrees into Master Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to execute this plan task-by-task. Keep every checkpoint and validation gate; do not resolve overlapping files by blindly taking one side.

**Goal:** Integrate every feature-bearing committed and uncommitted change from the active VD worktrees into `master` without losing source files, tests, required model artifacts, local master commits, or the behavior represented by the worktree-specific changes.

**Architecture:** Treat `master` as the integration target and first preserve each dirty worktree as a named, reversible capture. Fast-forward the committed ramp-speed history because `master` is its direct ancestor, then reapply and reconcile dirty states in dependency order: ramp solver/core vehicle interfaces, master tire/LC0 work, adaptive DOE, and telemetry entrypoint relocation. Keep legacy scripts and generated artifacts until their replacement is proven; make cleanup a separate post-merge decision.

**Tech Stack:** Git linked worktrees; MATLAB vehicle dynamics, ramp-speed, DOE, and LC0 tests; Python 3.12–3.14 QSS telemetry viewer with `pytest`/`pytest-qt` and `uv`.

**Spec:** Repository state inspected on 2026-09-27 in `C:\VD`; no separate merge specification was supplied.

## Global Constraints

- Target branch is `master` at `61ce55d`; it is already 8 commits ahead of `origin/master`, and those local commits must remain in history.
- `codex/ramp-speed-app` is 34 commits ahead of `master` and 0 commits behind it; its committed history is eligible for a fast-forward after dirty-state capture.
- `adaptive-doe` and `codex/telemetry-viewer` have no commits unique to their branch tips relative to `master`; their dirty working trees still contain changes that must be captured and reconciled.
- The detached `C:\Users\johny\.codex\worktrees\4a35\VD` checkout is clean at `4ca2627`, a prefix of the ramp-speed history; do not apply it a second time.
- Do not delete, reset, clean, archive, or drop a worktree/stash until the final feature inventory and validation gates pass.
- Do not commit Python bytecode, pytest caches, MATLAB runtime caches, or an unreviewed DOE checkpoint merely because `git stash --include-untracked` captured them.
- Preserve the old `Full Car Models` launcher paths while the new `Full Car Models/entrypoints` wrappers are being validated; the telemetry relocation’s deletions are not an automatic license to remove feature scripts.
- Use `git -c safe.directory=C:/VD -C C:\VD ...` for repository commands in this checkout.

## Review Focus

- **Dirty-state loss:** every tracked, untracked, and deleted path must be present in a named capture before any fast-forward; verify with per-worktree status and capture statistics in Task 1.
- **Vehicle-interface collision:** the ramp branch adds continuous-ratio/gear overrides, ride-height warm starts, and cached camber interpolation while master adds nonpositive-load handling and LC0 model selection; verify with the ramp contract tests plus `test_tireNonpositiveLoad`, `test_buildTireModel`, and `test_nondimensionalTire` in Tasks 3 and 5.
- **Divergent LC0 copies:** adaptive-doe contains older/different copies of overlapping LC0 functions while master contains target-scaling, plotting, uncertainty, donor-validation, and comparison artifacts; use a three-way review and run the complete LC0 test directory in Task 5.
- **Entrypoint relocation deletes:** telemetry-viewer deletes root DOE/ramp/simulation scripts while adding wrappers; retain both locations until `test_entrypointBootstrap`, adaptive DOE tests, ramp tests, and reference searches pass in Task 6.
- **Generated/project artifacts:** the ramp branch carries 722 MATLAB project resource files and the worktrees contain figures, MAT/CSV outputs, `uv.lock`, and bytecode; classify them against actual runtime/test references and finish with `git diff --check`, clean-status checks, and visual QA in Task 7.

### Task 1: Freeze and capture every dirty worktree

**Files:**
- Read-only inventory: `C:\VD\.git\worktrees\`, all five paths from `git worktree list --porcelain`
- External capture directory: `%TEMP%\vd-merge-capture-20260927` (outside the repository)
- No source files are modified in this task.

**Interfaces:**
- Produces one capture identifier per dirty checkout: `master`, `adaptive-doe`, `ramp-speed-app`, and `telemetry-viewer`.
- Preserves the existing `codex-pre-adaptive-doe-merge` stash; it must not be overwritten or dropped.

- [ ] **Step 1: Record the starting topology and existing stash list**

```powershell
$capture = Join-Path ([System.IO.Path]::GetTempPath()) 'vd-merge-capture-20260927'
New-Item -ItemType Directory -Path $capture -Force | Out-Null
git -c safe.directory=C:/VD -C C:\VD worktree list --porcelain | Set-Content -LiteralPath (Join-Path $capture 'worktrees.txt')
git -c safe.directory=C:/VD -C C:\VD status --short --branch --untracked-files=all | Set-Content -LiteralPath (Join-Path $capture 'master-status.txt')
git -c safe.directory=C:/VD -C C:\VD stash list --date=local | Set-Content -LiteralPath (Join-Path $capture 'stashes-before.txt')
```

Expected: the inventory lists `C:\VD`, `.worktrees\adaptive-doe`, `.worktrees\ramp-speed-app`, `.worktrees\telemetry-viewer`, and the detached `4a35` checkout; the pre-existing stash remains visible.

- [ ] **Step 2: Capture each dirty checkout with a named stash**

```powershell
$captures = [ordered]@{
    'master' = 'C:\VD'
    'adaptive-doe' = 'C:\VD\.worktrees\adaptive-doe'
    'ramp-speed-app' = 'C:\VD\.worktrees\ramp-speed-app'
    'telemetry-viewer' = 'C:\VD\.worktrees\telemetry-viewer'
}
foreach ($item in $captures.GetEnumerator()) {
    git -c safe.directory=C:/VD -C $item.Value stash push --include-untracked --message "codex-merge-capture-$($item.Key)-20260927"
    git -c safe.directory=C:/VD -C $item.Value status --short --branch --untracked-files=all | Set-Content -LiteralPath (Join-Path $capture "$($item.Key)-after-status.txt")
    git -c safe.directory=C:/VD -C $item.Value stash list -n 1 --format='%H %s' | Set-Content -LiteralPath (Join-Path $capture "$($item.Key)-stash.txt")
}
```

Expected: each dirty worktree is clean after capture; the detached worktree is unchanged and produces no new stash. If a worktree reports no changes, record `CLEAN` rather than inventing a stash id.

- [ ] **Step 3: Verify every capture before proceeding**

```powershell
foreach ($item in $captures.GetEnumerator()) {
    $stashLine = Get-Content -LiteralPath (Join-Path $capture "$($item.Key)-stash.txt") -ErrorAction SilentlyContinue
    if ($stashLine) {
        $stashId = ($stashLine -split ' ')[0]
        git -c safe.directory=C:/VD -C $item.Value stash show --stat --include-untracked $stashId
    }
    git -c safe.directory=C:/VD -C $item.Value status --short --branch --untracked-files=all
}
```

Expected: each stash statistic accounts for the status captured in Step 1, and no dirty source change has disappeared. Stop if counts or named paths do not match.

### Task 2: Establish the integration base and preserve committed history

**Files:**
- Git refs only: `master`, `codex/ramp-speed-app`, `adaptive-doe`, `codex/telemetry-viewer`
- Safety ref: `codex/merge-base-20260927`

**Interfaces:**
- Consumes the clean post-capture worktrees and the capture ids from Task 1.
- Produces a `master` fast-forward containing all 34 committed ramp-speed commits while retaining the existing 8 local master commits.

- [ ] **Step 1: Create a recoverable pre-merge ref**

```powershell
git -c safe.directory=C:/VD -C C:\VD branch codex/merge-base-20260927 master
git -c safe.directory=C:/VD -C C:\VD rev-parse master
git -c safe.directory=C:/VD -C C:\VD rev-list --left-right --count origin/master...master
```

Expected: the safety ref points to `61ce55d`, and the local/remote count remains `0 8`.

- [ ] **Step 2: Reconfirm source relationships**

```powershell
git -c safe.directory=C:/VD -C C:\VD rev-list --left-right --count master...codex/ramp-speed-app
git -c safe.directory=C:/VD -C C:\VD rev-list --left-right --count master...adaptive-doe
git -c safe.directory=C:/VD -C C:\VD rev-list --left-right --count master...codex/telemetry-viewer
git -c safe.directory=C:/VD -C C:\VD merge-base --is-ancestor 4ca2627aa3483dbfbfdf5b6d1c6cc583375e9c8b codex/ramp-speed-app
```

Expected: `0 34`, `11 0`, `3 0`, and a successful ancestry check respectively. The adaptive and telemetry committed histories require no second merge; their captured dirty states still do.

- [ ] **Step 3: Fast-forward master to the ramp branch**

```powershell
git -c safe.directory=C:/VD -C C:\VD merge --ff-only codex/ramp-speed-app
git -c safe.directory=C:/VD -C C:\VD log --oneline --decorate -n 5 master
git -c safe.directory=C:/VD -C C:\VD status --short --branch --untracked-files=all
```

Expected: `master` advances to `25ebb37` with the ramp history intact and remains clean. Do not use a squash merge or rebase; the ordered ramp commits are the audit trail for the solver/app feature.

### Task 3: Reapply the ramp-speed working tree and reconcile vehicle interfaces

**Files:**
- Apply the `ramp-speed-app` capture to `Full Car Models/+rampSpeed/`, `Full Car Models/RampSpeedApp.mlapp`, `Full Car Models/runRampSpeedStudy.m`, `Full Car Models/events/max_long_accel.m`, and their tests/plans/specs.
- Reconcile `Full Car Models/carComponents/Car.m`, `Full Car Models/carComponents/Tire2.m`, `Full Car Models/utilities/parameters_loop.m`, and `Full Car Models/utilities/buildCarFromParameterSet.m`.
- Preserve `Full Car Models/carComponents/Aero.m`, `Full Car Models/carComponents/AeroMap.m`, `Full Car Models/carConfigBaseline.m`, and the ramp support files.

**Interfaces:**
- Consumes `codex/ramp-speed-app` at `25ebb37` and its capture from Task 1.
- Produces ramp code that retains explicit-gear and continuous-envelope evaluation, status-aware results, ride-height context, aero-map evaluation, cancellation/validation, and the current dirty execution rewrite.

- [ ] **Step 1: Apply the ramp capture at its original base**

```powershell
$rampStash = (Get-Content -LiteralPath (Join-Path $capture 'ramp-speed-app-stash.txt') -ErrorAction Stop) -split ' ' | Select-Object -First 1
git -c safe.directory=C:/VD -C C:\VD stash apply $rampStash
```

Expected: ramp-specific files and current untracked ramp tests/specs reappear; any conflict is limited to files touched by later captures, not silently discarded.

- [ ] **Step 2: Resolve the vehicle-interface files by behavior**

In `Car.m`, keep the ramp branch’s optional `gearOverride`/`continuousRatio`, drivetrain reduction override, ride-height context/warm-start propagation, and status diagnostics. Reapply master’s tighter ride-height convergence (`toleranceIn = 1e-12`, `maxIterations = 20`) on the same solver path. In `Tire2.m`, keep the ramp branch’s injectable camber data/interpolator and warning behavior, then retain master’s `invalidLoad` clamp and zero-force behavior for zero and negative normal loads. In `parameters_loop.m`, combine the ramp helper construction with master’s `buildTireModel` selection so the ramp path can construct either `Tire2` or `NondimensionalTire` without bypassing ride-height/camber configuration.

```powershell
git -c safe.directory=C:/VD -C C:\VD diff --check
git -c safe.directory=C:/VD -C C:\VD status --short --untracked-files=all
```

Expected: no conflict markers remain; both `buildCarFromParameterSet` and `buildTireModel` are reachable from the parameter loop, and no `Tire2` call reintroduces `abs(F_z)` for contact-loss loads.

- [ ] **Step 3: Run the focused ramp and vehicle tests before adding other captures**

```matlab
cd('C:\VD');
addpath('Full Car Models');
setup_paths;
results = runtests({'Full Car Models/tests/test_rampSpeedRequestContract.m', ...
    'Full Car Models/tests/test_rampSpeedExecutionContract.m', ...
    'Full Car Models/tests/test_rampSpeedContinuousEnvelope.m', ...
    'Full Car Models/tests/test_rampSpeedContinuousEnvelopeSolver.m', ...
    'Full Car Models/tests/test_rampSpeedSession.m', ...
    'Full Car Models/tests/test_rampStudyRunner.m', ...
    'Full Car Models/tests/test_tireNonpositiveLoad.m'});
assertSuccess(results);
```

Expected: all focused tests pass on the combined ramp/vehicle interfaces. Do not proceed on a partial pass; repair the interface or update the owning test with an explicit compatibility reason.

### Task 4: Reapply master’s own dirty feature set

**Files:**
- `Full Car Models/carConfig.m`
- `Full Car Models/carComponents/Car.m`
- `Full Car Models/carComponents/Tire2.m`
- `Full Car Models/utilities/parameters_loop.m`
- `Full Car Models/carComponents/NondimensionalTire.m`
- `Full Car Models/utilities/buildTireModel.m`
- `Full Car Models/tests/test_buildTireModel.m`
- `Full Car Models/tests/test_nondimensionalTire.m`
- `Full Car Models/tests/test_tireNonpositiveLoad.m`
- `Magic Formula/experimental/lc0_nd_tire/` current master files, tests, and comparison results

**Interfaces:**
- Consumes the ramp-compatible vehicle interfaces from Task 3 and the `master` capture.
- Produces switchable legacy/LC0 tire construction, nonpositive-load safety, and the master’s current LC0 evaluation/plot/uncertainty additions.

- [ ] **Step 1: Apply the master capture after the ramp fast-forward**

```powershell
$masterStash = (Get-Content -LiteralPath (Join-Path $capture 'master-stash.txt') -ErrorAction Stop) -split ' ' | Select-Object -First 1
git -c safe.directory=C:/VD -C C:\VD stash apply $masterStash
```

Expected: master-only Car/Tire/LC0 changes and required untracked source/tests reappear; generated `__pycache__` files may also reappear but remain excluded from commits.

- [ ] **Step 2: Resolve the model-selection contract**

Keep `carConfig`’s default `legacy` behavior, its explicit `tireModel`/`uncertaintyMode` selectors, and its `model_artifact` path. Ensure `parameters_loop` passes those fields into `buildTireModel` and then into the ramp-compatible car constructor. Ensure `NondimensionalTire` keeps zero-load behavior, finite clamping diagnostics, pressure/camber diagnostics, and uncertainty-scenario validation.

- [ ] **Step 3: Commit the reconciled master/tire slice before the next subsystem**

```powershell
git -c safe.directory=C:/VD -C C:\VD add -- 'Full Car Models/carConfig.m' 'Full Car Models/carComponents/Car.m' 'Full Car Models/carComponents/Tire2.m' 'Full Car Models/carComponents/NondimensionalTire.m' 'Full Car Models/utilities/buildTireModel.m' 'Full Car Models/utilities/parameters_loop.m' 'Full Car Models/tests/test_buildTireModel.m' 'Full Car Models/tests/test_nondimensionalTire.m' 'Full Car Models/tests/test_tireNonpositiveLoad.m'
git -c safe.directory=C:/VD -C C:\VD commit -m "feat: integrate switchable tire models with ramp interfaces"
```

Expected: the commit contains source and tests only for this reconciled slice; LC0 result artifacts are staged in Task 5 after their reference policy is checked.

### Task 5: Reconcile adaptive-doe’s LC0 and DOE working tree

**Files:**
- LC0 overlap: `Magic Formula/experimental/lc0_nd_tire/README.md`, `lc0NDBuildModel.m`, `lc0NDCompareForceModels.m`, `lc0NDConfig.m`, `lc0NDEvaluate.m`, three `run_lc0_nd_*.m` scripts, and their tests.
- Master-only LC0 extensions: `lc0NDFitTargetScaling.m`, `lc0NDPlotModelComparisons.m`, `lc0NDUncertaintyScenarios.m`, `lc0NDValidateDonorChoice.m`, and their tests/results.
- Adaptive DOE: `Full Car Models/DOEStudyConfig.m`, `Full Car Models/DOE_Fitting.m`, `Full Car Models/SteadyStateLapsim.m`, `Full Car Models/runAdaptiveDOE.m`, `Full Car Models/utilities/doe*.m`, `Full Car Models/tests/test_doe_*.m`, `ADAPTIVE_DOE_WORKFLOW.md`, and the adaptive DOE design documents.

**Interfaces:**
- Consumes the master canonical LC0 path plus the adaptive capture.
- Produces one coherent LC0 API consumed by `NondimensionalTire` and one DOE analysis API consumed by `DOE_Fitting`, `doeSensitivityViewer`, and `doeInteractionSurface`.

- [ ] **Step 1: Apply the adaptive capture and identify true source conflicts**

```powershell
$adaptiveStash = (Get-Content -LiteralPath (Join-Path $capture 'adaptive-doe-stash.txt') -ErrorAction Stop) -split ' ' | Select-Object -First 1
git -c safe.directory=C:/VD -C C:\VD stash apply $adaptiveStash
git -c safe.directory=C:/VD -C C:\VD status --short --untracked-files=all
```

Expected: adaptive DOE utilities/tests are added; overlapping LC0 files report conflicts or divergent content instead of silently replacing the master versions.

- [ ] **Step 2: Resolve LC0 overlaps from the master feature superset**

Use a three-way diff for every overlapping LC0 file. Start from the current master version because it contains the target-scaling, comparison plotting, uncertainty scenarios, donor provenance, and vehicle-facing artifact behavior. Transplant any adaptive-only behavior only when it has a distinct documented contract and a test; preserve both test families when they assert different behavior. Do not replace the master comparison artifact with the older adaptive copy.

```powershell
git -c safe.directory=C:/VD -C C:\VD diff --cc -- 'Magic Formula/experimental/lc0_nd_tire'
git -c safe.directory=C:/VD -C C:\VD diff --check
```

- [ ] **Step 3: Keep the DOE feature set and make paths portable**

Retain resumable DOE state, parallel/serial guards, response fitting, Sobol sensitivity, interaction surfaces, the sensitivity viewer, and the explicit small-study configuration. Ensure scripts use `setup_paths` or the entrypoint bootstrap consistently and that `doeFindResultPath` resolves the same result locations used by tests. Keep `DOE_checkpoint.mat` in the capture for recovery, but stage it only if a test or documented fixture requires that exact state; otherwise leave it outside the source commit.

- [ ] **Step 4: Run the LC0 and DOE test slices**

```matlab
cd('C:\VD');
addpath('Full Car Models');
setup_paths;
lc0 = runtests('Magic Formula/experimental/lc0_nd_tire/tests');
doe = runtests({'Full Car Models/tests/test_adaptive_doe_resume.m', ...
    'Full Car Models/tests/test_doe_adaptive_sampling.m', ...
    'Full Car Models/tests/test_doe_fitting_helpers.m', ...
    'Full Car Models/tests/test_doe_interaction_surface.m', ...
    'Full Car Models/tests/test_doe_metrics.m', ...
    'Full Car Models/tests/test_doe_sensitivity_viewer.m', ...
    'Full Car Models/tests/test_doe_sobol.m', ...
    'Full Car Models/tests/test_nondimensionalTire.m', ...
    'Full Car Models/tests/test_buildTireModel.m'});
assertSuccess([lc0; doe]);
```

Expected: all LC0 and DOE tests pass, including both nominal and uncertainty paths, without requiring the generated checkpoint to be present.

### Task 6: Reapply telemetry-viewer’s entrypoint and viewer changes without dropping DOE/ramp features

**Files:**
- `Full Car Models/entrypoints/` new bootstrap and launcher scripts
- Existing root launchers listed as deleted in the telemetry capture, including `DOE_Fitting.m`, `SteadyStateLapsim.m`, and `runRampSpeedStudy.m`
- `Full Car Models/tests/test_entrypointBootstrap.m`, `test_adaptive_grip_sweep.m`
- `.gitignore`, `ADAPTIVE_DOE_WORKFLOW.md`, `telemetry_viewer/qss_telemetry/app.py`, `telemetry_viewer/uv.lock`

**Interfaces:**
- Consumes the reconciled ramp and DOE sources from Tasks 3–5.
- Produces portable MATLAB entrypoints and the already-merged QSS viewer interaction behavior without breaking direct legacy launch paths.

- [ ] **Step 1: Apply the telemetry capture**

```powershell
$telemetryStash = (Get-Content -LiteralPath (Join-Path $capture 'telemetry-viewer-stash.txt') -ErrorAction Stop) -split ' ' | Select-Object -First 1
git -c safe.directory=C:/VD -C C:\VD stash apply $telemetryStash
```

- [ ] **Step 2: Resolve relocation conflicts conservatively**

Keep the new `entrypoints/bootstrap.m` and launcher copies, including the entrypoint tests. Keep the original root scripts while references and tests are being checked; do not accept the telemetry stash’s root-file deletions if they would remove the adaptive DOE or ramp behavior. Combine `ADAPTIVE_DOE_WORKFLOW.md` rather than choosing the shorter copy. Union the `.gitignore` rules, including `baseline*.qss.h5`, `.venv/`, `__pycache__/`, `*.py[cod]`, `.pytest_cache/`, and `*.egg-info/`.

```powershell
rg -n --glob '*.m' --glob '*.md' 'DOE_Fitting|SteadyStateLapsim|runRampSpeedStudy|entrypoints/bootstrap' 'C:\VD'
git -c safe.directory=C:/VD -C C:\VD diff --check
```

Expected: every old or new launcher reference resolves to an existing file, and the root scripts are not removed until a separate cleanup review approves the move.

- [ ] **Step 3: Run MATLAB entrypoint and Python viewer tests**

```matlab
cd('C:\VD');
addpath('Full Car Models');
setup_paths;
results = runtests({'Full Car Models/tests/test_entrypointBootstrap.m', ...
    'Full Car Models/tests/test_adaptive_grip_sweep.m', ...
    'Full Car Models/tests/test_qssTelemetryExport.m'});
assertSuccess(results);
```

```powershell
Set-Location 'C:\VD\telemetry_viewer'
uv sync --extra dev
uv run pytest -q
```

Expected: MATLAB entrypoint/bootstrap/export tests pass and the Python telemetry suite passes without committing bytecode. If `uv` is unavailable, use the project’s configured Python 3.12–3.14 environment with `python -m pytest -q` and record the exact interpreter.

### Task 7: Full validation, visual QA, and finalize master

**Files:**
- All merged source and tests
- `docs/superpowers/plans/2026-09-27-merge-all-worktrees-into-master.md`
- `Full Car Models/RampSpeedApp.prj`, `Full Car Models/resources/project/`, and ramp-generated project metadata
- `.gitignore` and any explicitly retained result artifacts

**Interfaces:**
- Consumes the fully reconciled integration tree and all focused-test results.
- Produces a validated `master` with no unresolved paths and a documented list of intentionally excluded generated files.

- [ ] **Step 1: Run the full MATLAB suite**

```matlab
cd('C:\VD');
addpath('Full Car Models');
setup_paths;
results = runtests('Full Car Models/tests');
assertSuccess(results);
```

Expected: the full test suite passes. Record failures by test name and fix the owning interface before finalizing; do not mask failures by excluding tests.

- [ ] **Step 2: Validate the MATLAB project and visual outputs**

Open `Full Car Models/RampSpeedApp.prj`, run the ramp app smoke path, and inspect explicit-gear versus continuous-envelope results, cancellation/status display, aero-map/ride-height traces, and exported figures. Open the LC0 comparison `.fig`/`.png` artifacts and verify legends, uncertainty cases, donor provenance, and no helper/reference lines are presented as measured model curves. Run the QSS exporter and Python viewer, then inspect absolute and Delta (active minus datum) modes, cursor ownership, units/warnings, and invalid-gap handling.

- [ ] **Step 3: Audit the final file set**

```powershell
git -c safe.directory=C:/VD -C C:\VD diff --check
git -c safe.directory=C:/VD -C C:\VD status --short --untracked-files=all
git -c safe.directory=C:/VD -C C:\VD status --porcelain=v1 | Select-String '^(DD|AU|UD|UA|DU|AA|UU) '
git -c safe.directory=C:/VD -C C:\VD ls-files --others --exclude-standard
```

Expected: no unmerged entries; only intentionally retained source/tests/docs/required artifacts are tracked. Python bytecode, caches, and an unneeded DOE checkpoint are either ignored or remain in the external capture, not in the final commit.

- [ ] **Step 4: Commit logical integration slices and preserve the ramp history**

Use separate commits for the reconciled vehicle/tire interfaces, LC0/DOE work, telemetry entrypoints, and generated-ignore/project metadata. Do not squash the 34 ramp commits or the existing local master commits. The final commit messages should identify the feature slice and its test command, for example:

```powershell
git -c safe.directory=C:/VD -C C:\VD commit -m "feat: integrate adaptive DOE and LC0 workflows"
git -c safe.directory=C:/VD -C C:\VD commit -m "feat: add portable MATLAB and telemetry entrypoints"
git -c safe.directory=C:/VD -C C:\VD log --oneline --decorate -n 12 master
```

- [ ] **Step 5: Keep rollback refs until the user approves cleanup**

Retain `codex/merge-base-20260927`, the four named capture stashes, the existing `codex-pre-adaptive-doe-merge` stash, and all active worktrees until the user confirms the merge is complete. Only then may the user choose whether to delete obsolete root launchers, drop captures, archive worktrees, or push `master`.

## Self-review checklist

- [x] Committed topology distinguishes already-merged adaptive/telemetry histories from the only ahead branch, ramp-speed.
- [x] Dirty working trees are captured before fast-forward and re-applied in an order that minimizes base mismatches.
- [x] The high-risk overlap files and their required behavior are named explicitly.
- [x] Generated files are preserved for recovery but excluded from automatic feature commits.
- [x] MATLAB, Python, focused, full-suite, and visual validation gates are specified.
