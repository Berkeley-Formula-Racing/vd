# Adaptive DOE and Parallel Ramp-Metrics Design

## Purpose

The vehicle DOE workflow must support two different phases of vehicle design without discarding earlier simulations:

1. broad sensitivity work, where the objective is to learn which parameters and interactions control the vehicle responses across the full design space; and
2. design optimization, where the objective is to find configurations that maximize modeled FSAE dynamic-event points subject to optional engineering penalties.

The same runner will support `sensitivity`, `optimization`, and `hybrid` acquisition modes. A saved sensitivity study can later resume in optimization mode. Parameter ranges will no longer be hard-coded in `carConfig.m`; they will be declared by name in a dedicated study configuration.

## Scope

This change will:

- move DOE parameter selection and ranges into a user-editable study definition;
- preserve `carConfig()` as the source of the baseline vehicle;
- generate an initial space-filling design and subsequent adaptive batches;
- fit Gaussian-process surrogates for adaptive point selection;
- maximize modeled dynamic-event points with configurable penalties;
- execute optional ramp metrics in parallel with each car case;
- checkpoint every completed batch and support resume and mode changes;
- cache per-case metrics so analysis does not repeatedly traverse every saved car; and
- preserve a final `DOE_results.mat` interface for the existing analysis workflow.

This change will not add an efficiency-event model, replace the QSS vehicle equations, or perform multi-objective Pareto optimization. Efficiency points will remain explicitly excluded until a validated efficiency model exists.

## User-facing workflow

`SteadyStateLapsim.m` remains the entry point, but becomes a thin wrapper:

```matlab
study = DOEStudyConfig();
runAdaptiveDOE(study);
```

The user edits only `DOEStudyConfig.m` to choose parameters, ranges, events, adaptive mode, penalties, worker count, ramp fidelity, output location, and resume behavior.

A typical staged workflow is:

1. set `study.mode = "sensitivity"` and `study.maxCases = 256`;
2. run the initial LHS and sensitivity-focused adaptive batches;
3. review sensitivities and narrow the design ranges if the physical design changes;
4. for the same unchanged design space, set `study.mode = "optimization"`, increase `study.maxCases`, and resume the checkpoint; and
5. directly simulate the best predicted configurations for final validation.

Changing the parameter names, baseline values, or parameter bounds creates a different design space and therefore requires a new checkpoint. Changing mode, objective penalties, workers, batch size, ramp storage detail, or maximum case count is allowed when resuming.

## Study configuration

The study configuration is a scalar structure. Its principal fields are:

```matlab
study.name         = "2026_design_study";
study.mode         = "sensitivity"; % sensitivity | optimization | hybrid
study.randomSeed   = 1;
study.initialCases = 128;
study.batchSize    = 64;
study.maxCases     = 512;
study.numWorkers   = 16;
study.resume       = true;

study.parameters = table( ...
    ["mass"; "cg_height"; "R_sf"; "cla"; "distribution"], ...
    [-5; -10; -15; -10; -5], ...
    [ 5;  10;  15;  10;  5], ...
    ["percent"; "percent"; "percent"; "percent"; "percent"], ...
    'VariableNames',{'name','lower','upper','rangeType'});
```

`rangeType="percent"` interprets the bounds as percentage deltas about the scalar baseline. It is allowed only for a positive, nonzero baseline. Parameters whose baseline is zero, negative, or intentionally crosses zero must use `rangeType="absolute"`. The resolved absolute lower and upper bounds are stored in the checkpoint and final design table.

Every parameter name is validated against the parameter registry already used to construct `Car`, `Aero`, `Powertrain`, and `Tire2`. Unknown names, duplicate names, nonfinite bounds, reversed bounds, and invalid percentage baselines fail before the parallel pool starts.

The default study will mirror the current DOE variables: `mass`, `wheelbase`, `weight_dist`, `track_width`, `cg_height`, `R_sf`, `accel_cda`, `accel_cla`, `final_drive`, `gamma_f`, `gamma_r`, and `p_i`. Camber ranges use absolute units because their baselines are negative. Users can add full-car `cda`, `cla`, `distribution`, grip scales, or any other registered parameter without editing `carConfig.m`. No DOE ranges remain embedded in `carConfig.m`.

## Baseline-car and explicit-case construction

`carConfig()` will continue to return one calibrated baseline lap car, one acceleration car, and `eventParams` for existing tools. Baseline physical values remain defined there as scalars.

An explicit-case construction path will be added so the adaptive runner can pass a table of absolute sampled values. The builder will start from the baseline parameter structures, override only the named fields, and use the existing `parameters_loop` construction logic. Derived quantities will therefore still be established by the normal constructors rather than by mutating completed `Car` objects.

Legacy no-argument calls remain compatible. For transition compatibility, `carConfig("LHS",N)` and `carConfig("Random",N)` will remain available, issue a deprecation warning, and construct one non-adaptive sample set from the ranges in `DOEStudyConfig.m`. New studies use `runAdaptiveDOE` because the legacy calls cannot describe adaptive batches or resume state.

## Adaptive sampling process

All predictors are represented internally in normalized coordinates on `[0,1]^p`. The stored `designTable` contains physical units, while the normalized matrix is retained for fitting, distance calculations, and reproducibility.

### Initial batch

The first batch is a repeatable maximin-style Latin-hypercube design using `study.randomSeed`. The baseline point is included explicitly if it is not already present. Duplicate normalized points are rejected.

### Surrogate models

After the initial batch, the runner fits one standardized ARD Gaussian-process regression model per configured adaptive response. With at most 512 cases, exact GPR is appropriate. An ARD Matérn 3/2 kernel is the default because it tolerates less-than-perfectly-smooth simulation responses better than a squared-exponential kernel.

Only finite, valid cases train a continuous-response model. A feasibility model is fitted when the results contain both valid and invalid cases; otherwise predicted feasibility is one everywhere. Surrogate fitting occurs on the MATLAB client after each batch and is not nested inside `parfor`.

### Sensitivity acquisition

Sensitivity mode ranks a large repeatable candidate pool by the root-mean-square normalized predictive standard deviation across the configured responses. Normalizing by each response's observed spread prevents a response with large units from dominating the acquisition.

This mode selects points that reduce uncertainty across the complete parameter space. It is intended to support response surfaces, interactions, and surrogate-based Sobol indices.

### Optimization acquisition

Optimization mode fits a GP to `objective_score` and ranks candidates by expected improvement for a maximization objective. Expected improvement is multiplied by predicted feasibility when a feasibility model is available. This balances predicted points and model uncertainty without treating an obviously invalid design as attractive merely because its prediction is uncertain.

### Hybrid acquisition

Hybrid mode selects a configurable fraction of each batch using sensitivity acquisition and the remainder using optimization acquisition. The default is 70% sensitivity and 30% optimization. The two selections share one exclusion set, so the same candidate cannot be selected twice.

### Batch diversity

Candidates are selected greedily from their acquisition ranking. After each selection, candidates closer than a configurable normalized distance to any completed or newly selected point are excluded. If the threshold prevents filling a batch, it is reduced deterministically until the batch is full. This avoids spending a parallel batch on nearly identical cars.

## Event points and penalties

Optimization maximizes:

```text
objective_score = modeled_dynamic_points - total_penalty_points
```

The raw point columns are retained separately:

- skidpad points;
- acceleration points;
- autocross points;
- endurance points; and
- modeled dynamic-event total.

The existing `Events2.computePoints()` formulas and `eventParams.winning_time` values remain the single source of scoring truth. The optimization runner must execute skidpad, acceleration, autocross, and endurance before calculating the full modeled dynamic-event total. The current efficiency-event comment is not treated as an implemented score: efficiency points are excluded and the output is named `modeled_dynamic_points` so it cannot be mistaken for the complete competition score.

Penalties are configured as named rules with an enabled flag, threshold, direction, and points-per-unit weight. Initial supported rules are:

- invalid or failed case;
- g-g or ramp solve-failure fraction;
- wheel lift or negative minimum wheel load;
- energy above a configured allowance;
- understeer-gradient bounds; and
- rebalance-speed error from a configured target.

Hard-invalid cases do not train continuous GPs and receive no optimization benefit. Other penalties remain soft and are stored individually so the vehicle-dynamics group can audit why a design's objective differs from its raw event points. Penalty weights are never hidden defaults: disabled rules have zero weight in `DOEStudyConfig.m`.

Changing penalty weights on resume recomputes `objective_score` from cached raw metrics without rerunning the cars.

## Parallel execution and ramp metrics

The outer case loop remains the only active parallel loop. Each worker:

1. receives one explicit car case;
2. builds its g-g diagram;
3. runs the configured events;
4. calculates event points;
5. optionally runs `rampSweep` on the same solved car;
6. reduces the result to one metric row; and
7. returns the car, ramp result, metrics, timing, and any error information.

Running ramps within the same outer `parfor` avoids the current serial `doeMetrics` ramp loop and avoids retransmitting every car to the pool. MATLAB nested parallelism is not introduced; inner g-g work executes serially when called from a worker.

Ramp configuration includes enabled state, speeds, ramp points, bisection count, mode, and whether full point-by-point ramp tables are retained. Per-speed summaries and the scalar DOE metrics are always saved when ramps are enabled. Full ramp points default to disabled for a 512-car study to control memory and file size.

If ramps were disabled during an earlier run and are enabled on resume, the runner performs a parallel backfill for completed cars before selecting the next adaptive batch.

## Metrics and analysis efficiency

Metric extraction will be split into a per-case reducer and an all-cases wrapper. The worker reduces each newly completed car immediately and the client appends those rows to the checkpoint. Adaptive iterations therefore do not rebuild `ggMetrics` for every old case.

`doeAnalyze` will consume cached metrics when present and retain its current extraction path for older result files. Final model fitting and plot generation occur once after the requested case count is reached, or explicitly when the user runs `DOE_Fitting.m` against a checkpoint.

The final analysis will support quadratic response surfaces and GPs. Cross-validation statistics determine which model is preferred per response. ARD length scales are diagnostic only; final global sensitivity rankings use surrogate-based Sobol indices rather than treating inverse GP length scale as a Sobol measure.

## Persistence and resume behavior

The output directory contains:

- `DOE_checkpoint.mat`: authoritative resumable state after every batch;
- `DOE_results.mat`: compatibility-oriented final results;
- `DOE_analysis.mat`: fitted models, sensitivity results, and plot inputs; and
- the configured figure directory.

The checkpoint stores:

- resolved study configuration and a design-space signature;
- baseline values and absolute parameter bounds;
- normalized and physical design matrices;
- batch number, acquisition mode, acquisition values, and random state;
- completed cars and event outputs;
- raw event points, penalties, objective score, and cached metrics;
- ramp summaries and optional full ramp data;
- per-case runtime, worker count, and failure diagnostics; and
- total elapsed time.

Checkpoint writes happen on the client after a batch completes. The file is written to a temporary MAT file and then replaced, preventing a partially written checkpoint from destroying the last good state. MAT v7.3 remains appropriate because the checkpoint can exceed 2 GB and benefits from structured partial loading.

On resume, the runner validates the design-space signature. Mode, maximum cases, workers, objective penalties, and storage detail may change. Baseline values, parameter names, range types, and resolved bounds may not change within the same checkpoint.

## Failure handling

A failed case records its identifier, sampled parameters, exception identifier and message, elapsed time, and solve-health metrics available before the failure. One failed case does not terminate a batch.

The adaptive engine requires a minimum number of valid cases greater than the number of varying predictors before fitting. If that condition is not met, it adds another space-filling batch rather than attempting an unreliable GP. If no response has enough finite variation to fit, the runner stops with an actionable diagnostic while preserving the checkpoint.

Missing required toolboxes are detected before the initial batch. Adaptive GP fitting requires Statistics and Machine Learning Toolbox; the parallel runner falls back to serial only when explicitly allowed by configuration, otherwise an unavailable Parallel Computing Toolbox is reported before simulation starts.

## Verification strategy

Implementation will follow test-driven development. Required tests are:

1. baseline regression: `carConfig()` constructs the same scalar baseline vehicle as before the refactor;
2. parameter resolution: percent and absolute ranges map normalized samples to expected physical values and reject invalid definitions;
3. explicit construction: sampled table values reach the corresponding `Car`, `Aero`, `Powertrain`, and tyre/grip properties;
4. deterministic sampling: the same seed and study definition produce identical initial and adaptive candidates;
5. acquisition behavior: sensitivity, optimization, and hybrid modes select the intended batch allocation without duplicates;
6. objective accounting: raw points, each penalty, and final objective reconcile exactly;
7. checkpoint resume: a mode change resumes without repeating completed cases, while a range change is rejected;
8. ramp equivalence: parallel per-case ramp summaries match the current serial calculation within numerical tolerance;
9. failure isolation: one intentionally invalid case is recorded and does not terminate the batch; and
10. end-to-end smoke test: a small initial batch plus one adaptive batch runs with two workers, saves, resumes, and produces analyzable results.

Before a production 512-car run, a small benchmark will compare serial and parallel ramp runtime and verify that saved metric rows and modeled points are complete.

## Migration

The current uncommitted DOE-analysis work remains the starting point. Migration proceeds without altering unrelated user changes:

1. add `DOEStudyConfig.m` and external parameter validation;
2. remove only the hard-coded DOE range block from `carConfig.m` and add explicit-case construction;
3. add adaptive sampling and checkpoint functions;
4. move per-case metric and ramp work into the outer parallel case execution;
5. update `SteadyStateLapsim.m` to call the new runner;
6. update `doeAnalyze` and `DOE_Fitting.m` to prefer cached metrics and add GP/Sobol outputs; and
7. preserve compatibility with existing saved `DOE_results.mat` files.

No existing user calibration values or unrelated working-tree edits will be overwritten during this refactor.
