# LC0 Experimental Evaluator and Pacejka Comparison Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build an untracked, switchable LC0 force evaluator from measured target lateral behavior plus normalized provisional donor longitudinal/coupling behavior, and compare it fairly against the current Pacejka model.

**Architecture:** The evaluator will consume saved experimental target and donor artifacts and return `Fx`, `Fy`, pure-force demands, and coupled-force utilization in SI units. A separate comparison runner will evaluate both tyre models at matched operating points and write only untracked MAT/PNG/FIG outputs.

**Tech Stack:** MATLAB R2026a, existing experimental TTC loaders and MATLAB tables.

**Spec:** `docs/superpowers/specs/2026-08-23-lc0-nondimensional-tire-design.md`

## Global Constraints

- Create files only under `Magic Formula/experimental/lc0_nd_tire/` and `docs/superpowers/`; do not modify tracked car-model code.
- The experimental evaluator is not connected to `Tire2` or vehicle simulation.
- Target LC0 lateral force comes only from Round 8 runs 15, 16, 18, and 19.
- Donor force level is controlled only through explicit `rhoMu` and `rhoStiff` fields; no implicit matching is allowed.
- Use test-first development and save exact input configuration with every result.

---

### Task 1: Build a target-lateral / donor-longitudinal evaluator

**Files:**
- Create: `Magic Formula/experimental/lc0_nd_tire/lc0NDEvaluate.m`
- Test: `Magic Formula/experimental/lc0_nd_tire/tests/test_lc0NDEvaluate.m`

**Interfaces:**
- Consumes: target-lateral curve, donor longitudinal curve, explicit capacities, `rhoMu`, `rhoStiff`, and coupling exponent.
- Produces: signed `Fx`, `Fy`, pure `Fx0`, `Fy0`, and normalized utilization.

- [ ] Write a failing synthetic test proving pure-axis force invariance and explicit scaling behavior.
- [ ] Run the test and confirm the evaluator is undefined.
- [ ] Implement bounded interpolation, extrapolation guards, and `lc0NDCombinedForce` coupling.
- [ ] Run the evaluator test and confirm it passes.

### Task 2: Create reproducible target-reference curves

**Files:**
- Create: `Magic Formula/experimental/lc0_nd_tire/lc0NDBuildTargetReference.m`
- Test: `Magic Formula/experimental/lc0_nd_tire/tests/test_lc0NDBuildTargetReference.m`

**Interfaces:**
- Consumes: Round 8 target samples and an explicit load/pressure/camber condition.
- Produces: lateral `mu_y(alpha)` curve with sample counts, coverage qualification, peak capacity, and slope magnitude.

- [ ] Write a failing synthetic test for load-normalized lateral curve construction.
- [ ] Run the test and confirm the builder is undefined.
- [ ] Implement the binned TTC-sign-preserving reference-curve builder.
- [ ] Run the test and confirm it passes.

### Task 3: Compare experimental evaluator with legacy Pacejka

**Files:**
- Create: `Magic Formula/experimental/lc0_nd_tire/run_lc0_nd_pacejka_comparison.m`
- Test: `Magic Formula/experimental/lc0_nd_tire/tests/test_run_lc0_nd_pacejka_comparison.m`
- Modify: `Magic Formula/experimental/lc0_nd_tire/README.md`

**Interfaces:**
- Consumes: target and donor artifacts plus the legacy `Tire2` interface through an adapter, at matched SI operating points.
- Produces: a results MAT file and pure lateral, pure longitudinal, and combined-slip overlays; labels baseline uncertainty and donor provenance.

- [ ] Write a failing synthetic comparison-runner test requiring an untracked result artifact and matched case table.
- [ ] Run the test and confirm the runner is undefined.
- [ ] Implement a legacy adapter only after inspecting `Tire2` force-interface inputs and sign convention; do not alter it.
- [ ] Plot experimental versus legacy forces on the same load/pressure/camber grid, and record unsupported operating points as missing rather than extrapolating.
- [ ] Run all package tests, run the comparison with the actual LC0 data, and inspect saved figures.
