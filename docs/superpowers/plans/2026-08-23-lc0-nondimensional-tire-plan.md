# LC0 Nondimensional Tire Prototype Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build an untracked, testable LC0 tire-model prototype from the supplied Round 8 TTC free-rolling data.

**Architecture:** A new experimental package loads and summarizes the measured LC0 lateral data separately from the existing car tire model. A bounded normalized combined-slip function accepts pure-force inputs and donor-derived calibration parameters but makes no unsupported claim that target combined-slip data exists.

**Tech Stack:** MATLAB R2026a, base MATLAB tables/MAT files, Statistics and Machine Learning Toolbox only when already present.

**Spec:** `docs/superpowers/specs/2026-08-23-lc0-nondimensional-tire-design.md`

## Global Constraints

- Create files only under `Magic Formula/experimental/lc0_nd_tire/` and its `tests/` subfolder.
- Do not modify or stage existing tracked model, configuration, parser, or data files.
- Target runs are Round 8 `A1965run15`, `16`, `18`, and `19`.
- Keep force/load units explicit at every package boundary.
- Use test-first development for every production function.

---

### Task 1: Target-data configuration and loader

**Files:**
- Create: `Magic Formula/experimental/lc0_nd_tire/lc0NDConfig.m`
- Create: `Magic Formula/experimental/lc0_nd_tire/lc0NDLoadFreeRolling.m`
- Test: `Magic Formula/experimental/lc0_nd_tire/tests/test_lc0NDLoadFreeRolling.m`

**Interfaces:**
- Consumes: Round 8 MAT files containing `SA`, `SL`, `FY`, `FZ`, `IA`, `P`, and `V` in USCS units.
- Produces: `data` struct with SI vectors and `metadata` describing the source runs.

- [ ] Write a failing test that uses a temporary TTC-like MAT file with zero slip and asserts N/deg/psi conversion and run metadata.
- [ ] Run the test and confirm `lc0NDLoadFreeRolling` is undefined.
- [ ] Implement configuration and loader validation: required channels, finite samples, target run match, and maximum absolute `SL`.
- [ ] Re-run the test and confirm it passes.

### Task 2: Measured lateral summary

**Files:**
- Create: `Magic Formula/experimental/lc0_nd_tire/lc0NDLateralSummary.m`
- Test: `Magic Formula/experimental/lc0_nd_tire/tests/test_lc0NDLateralSummary.m`

**Interfaces:**
- Consumes: SI target `data` from `lc0NDLoadFreeRolling`.
- Produces: one row per requested load/pressure/camber bin with peak friction and near-zero-slip lateral stiffness.

- [ ] Write a failing synthetic-data test that expects a known peak `abs(Fy)/Fz` and initial slope.
- [ ] Run the test and confirm `lc0NDLateralSummary` is undefined.
- [ ] Implement bin selection and linear near-zero-slip stiffness fit, retaining sample count and actual bin means.
- [ ] Re-run the test and confirm it passes.

### Task 3: Bounded normalized combined-force law

**Files:**
- Create: `Magic Formula/experimental/lc0_nd_tire/lc0NDCombinedForce.m`
- Test: `Magic Formula/experimental/lc0_nd_tire/tests/test_lc0NDCombinedForce.m`

**Interfaces:**
- Consumes: pure `Fx0`, `Fy0`, capacities `muX*Fz`, `muY*Fz`, and coupling exponent.
- Produces: coupled `Fx`, `Fy`, and a utilization value.

- [ ] Write failing tests for pure-axis invariance and anisotropic-limit compliance.
- [ ] Run the tests and confirm `lc0NDCombinedForce` is undefined.
- [ ] Implement normalized-demand coupling with an explicit smoothness exponent and numerical zero-demand guard.
- [ ] Re-run the tests and confirm they pass.

### Task 4: Data-backed prototype runner and diagnostics

**Files:**
- Create: `Magic Formula/experimental/lc0_nd_tire/run_lc0_nd_prototype.m`
- Create: `Magic Formula/experimental/lc0_nd_tire/README.md`
- Test: `Magic Formula/experimental/lc0_nd_tire/tests/test_lc0NDPrototypePaths.m`

**Interfaces:**
- Consumes: `lc0NDConfig`, target TTC files, and tasks 1–3.
- Produces: untracked `results/lc0_nd_target_summary.mat` and PNG diagnostic figures.

- [ ] Write a failing path test that requires results to live only under the experimental package.
- [ ] Run the test and confirm the runner is undefined.
- [ ] Implement the runner, result save path, source/run manifest, lateral plots, and coupling-envelope plot.
- [ ] Run the package tests and the runner on Round 8 data.
- [ ] Inspect output paths with `git status --short` and confirm no tracked source file changed.
