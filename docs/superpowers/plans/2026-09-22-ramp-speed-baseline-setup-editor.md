# Ramp-Speed Baseline and Setup Editor Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task.

**Goal:** Replace placeholder ramp-speed cars with a stable single-car baseline and add reproducible setup duplication/editing for suspension, ride height, and aero-map variants.

**Architecture:** `carConfigBaseline.m` owns frozen ramp-speed defaults and editable option lists. A single-car builder creates one solver-ready `Car` from a serializable setup specification. `RampSpeedApp.mlapp` stores one car per setup, rebuilds duplicated setups from specs, and persists specs with study cases.

**Tech Stack:** MATLAB R2026a, App Designer `.mlapp`, MATLAB unit tests, existing `Car`, `Aero`, and `AeroMap` classes.

**Spec:** User-provided Ramp-Speed Baseline and Setup Editor Plan in the conversation.

## Global Constraints

- `carConfigBaseline.m` must not call `carConfig()`.
- Each setup factory/build operation creates exactly one solver-ready `Car` object.
- The first baseline uses acceleration-oriented car parameters for the pure-longitudinal priority.
- ARB and spring inputs are discrete baseline-defined choices; ride heights are numeric and stored internally in inches.
- Rear ARB and spring selections derive `R_sf` from front/rear total roll stiffness; full transient suspension solving is out of scope.
- Aero maps use an explicit catalog and the existing `AeroMap` CSV contract.
- Existing studies without setup specifications remain loadable as read-only result data.
- No new Codex project chats; use subagents only.

## Review Focus

- The baseline must construct a real `Car`, not a name/source placeholder.
- An N-by-1 car catalog must select the same setup for lateral and longitudinal modes.
- Changing rear ARB/springs must change derived `R_sf` in the expected direction.
- Map IDs must remain portable and resolve through the explicit catalog.
- Save/load must preserve setup specs without serializing raw `Car` objects.

### Task 1: Stable baseline and single-car construction

Create `carConfigBaseline.m` and extract reusable single-car construction from `utilities/parameters_loop.m`. Add tests proving the baseline is independent of `carConfig`, returns one valid acceleration-oriented `Car`, and exposes serializable setup metadata/options.

### Task 2: Setup builder, derived roll split, and aero catalog

Create setup-spec validation/building, roll-stiffness derivation, and an explicit aero-map catalog. Add tests for duplicate specs, units, invalid choices, ride-height propagation, map validation, and the `R_sf` formula.

### Task 3: One-car study selection and persistence

Update case selection and study metadata so setups are an N-by-1 car catalog shared by both ramp types. Persist setup specs/map IDs/config version, rebuild editable cars on load, and preserve read-only loading for legacy studies. Add runner/schema regression tests.

### Task 4: RampSpeedApp setup editor

Update `RampSpeedApp.mlapp` with protected baseline, duplicate/edit/delete workflow, ARB/spring dropdowns, numeric ride heights, aero-map dropdown, setup table/provenance, and overlay selection. Add app smoke/request tests covering default construction, both ramp types, duplication, validation, and persistence.

### Task 5: Integrated verification

Run focused new tests, app smoke/request tests, the broader ramp regression suite, archive/XML validation for the `.mlapp`, and `git diff --check`. Review the complete branch before reporting completion.

## As-built integration notes (2026-09-27)

The baseline/setup editor is integrated with the execution rewrite. RampSpeedApp.mlapp delegates setup duplication, editing, deletion, save/load, run, cancel, and result clearing to rampSpeed.RampSpeedSession; the App is not a second numerical runner. rampSpeed.StudyExecutor and rampSpeed.runStudy consume the same one-Car-per-setup catalog for both ramp types.

Driver weight and rear weight distribution are serializable setup fields alongside rear ARB stiffness, spring selections, ride heights, aero-map ID, and baseline/configuration version. Saved legacy studies without setup specifications load as read-only result data. The canonical function wrapper is runRampSpeedStudy(options), and tests/test_rampSpeedEndToEndRealCar.m verifies a real baseline setup at representative speeds in both modes.
