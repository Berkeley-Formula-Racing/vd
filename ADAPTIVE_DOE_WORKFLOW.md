# Adaptive DOE workflow

This branch provides a checkpointed, adaptive design-of-experiments (DOE)
workflow for the quasi-steady-state vehicle model. It builds a
space-filling initial design, runs lap and optional ramp metrics, then fits
Gaussian-process surrogate models to choose useful subsequent designs.

## Before running

Use this workflow from the `adaptive-doe` branch and start MATLAB in the
`Full Car Models` folder:

```matlab
cd("Full Car Models")
setup_paths
```

The adaptive portion requires the Statistics and Machine Learning Toolbox
(`fitrgp`), and parallel runs require Parallel Computing Toolbox. A local
worker pool must contain at least as many workers as the requested study
uses. For example:

```matlab
delete(gcp('nocreate'))
parpool('local',12)
```

`study.numWorkers` controls the outer `parfor` loop. It cannot create more
physical workers than the active pool or the local cluster profile allows.
Use `study.numWorkers = 0` for a serial run.

## Configure a study

Open `DOEStudyConfig.m` and save a copy of its settings in a launch script
for each named study. The default configuration is a 12-parameter
sensitivity study with 128 initial designs, 12 workers, and a maximum of
512 designs.

Create a separate output directory for every meaningfully different study:

```matlab
study = DOEStudyConfig();
study.name = "design_screen_v1";
study.mode = "sensitivity";
study.initialCases = 128;
study.batchSize = 16;
study.maxCases = 384;
study.numWorkers = 12;
study.output.directory = fullfile(fileparts(which('DOEStudyConfig')), ...
    "DOE_output_design_screen_v1");
```

### Parameter ranges

`study.parameters` is the complete external definition of the design
space. It is a table with these columns:

| Column | Meaning |
| --- | --- |
| `name` | A parameter registered by `parameters_loop.m` |
| `lower`, `upper` | Lower and upper variation bounds |
| `rangeType` | `"percent"` relative to nominal `carConfig`, or `"absolute"` |

The supplied starting set is:

`mass`, `wheelbase`, `weight_dist`, `track_width`, `cg_height`, `R_sf`,
`accel_cda`, `accel_cla`, `final_drive`, `gamma_f`, `gamma_r`, and `p_i`.

Only include parameters with independent and physically plausible ranges.
For example, do not independently sweep a mass distribution and component
masses if they are meant to describe the same packaging change.

### Modes

- `"sensitivity"` — recommended first. Prioritizes designs that improve
  response-surface coverage for the configured responses.
- `"optimization"` — prioritizes the score from `study.objective`.
- `"hybrid"` — combines both behaviours. Use after the important variables
  and credible constraints are known.

For a design screen, begin with `sensitivity`. Do not tune objective
penalties until the model has been calibrated sufficiently to trust their
absolute values.

### Events and ramp metrics

All four events are needed when dynamic points are an analysis response:

```matlab
study.events = ["skidpad","accel","autocross","endurance"];
```

Ramp metrics provide understeer gradients and rebalance speed. The normal
handling screen is:

```matlab
study.ramps.enabled = true;
study.ramps.speeds = [10 25];
study.ramps.nRamp = 8;
study.ramps.nBisect = 3;
study.ramps.mode = "balanced";
study.ramps.saveFullPoints = false;
```

For a cheaper broad screen, set `study.ramps.enabled = false` and run a
separate, smaller ramp-focused study for the selected parameters. Enabling
ramps during a resume intentionally backfills compatible cached cases; it
does not only calculate new cases.

## Run the study

Run directly:

```matlab
state = runAdaptiveDOE(study);
```

Or run the thin convenience script after editing `DOEStudyConfig.m`:

```matlab
SteadyStateLapsim
```

The runner writes a checkpoint after the initial design, before/after
selection, and after every completed batch. A failed MATLAB session can
therefore be resumed safely:

```matlab
study.resume = true;
state = runAdaptiveDOE(study);
```

When resuming, these settings may change: `mode`, `maxCases`,
`numWorkers`, `batchSize`, `objective`, and `ramps`. Do **not** change the
parameter list/bounds, event list, or adaptive-selection settings for an
existing output directory; create a new study directory instead.

## Files produced

The output directory contains:

| File | Purpose |
| --- | --- |
| `DOE_checkpoint.mat` | Restartable state, including pending designs and cached cars |
| `DOE_results.mat` | Final/exported design table, metrics, point data, and selection history |
| `DOE_analysis.mat` | Offline analysis output, created by `DOE_Fitting` |
| `figures/doe/*.png` | Saved analysis figures |

The global run ledger is `Full Car Models/sim_log.csv`. It records the DOE
run name, start/finish timestamps, events, requested workers, case count,
status, output paths, and elapsed time.

## Analyze without rerunning simulations

Open and run `DOE_Fitting.m`. It automatically selects the configured DOE
output, or the newest `DOE_output*` result folder when the configured one
does not exist. It does not launch new event simulations when
`runRampMetrics` is false.

```matlab
DOE_Fitting
```

One press of MATLAB's **Run** button generates and saves the configured
static plots plus one tiled, rotatable interaction figure. Each 3-D tile
uses a fitted response and its two largest total-order Sobol inputs. It also
opens one interactive DOE sensitivity viewer: select the output, x input,
and y input from dropdowns to see the Sobol ranking, main effect, and
interaction surface. Set `resultSource` near the top of `DOE_Fitting.m` to a
specific output folder or MAT-file when several studies are present.

Useful default responses include dynamic points, Autocross/Accel/Skidpad
time, total energy/work, g-g metrics, understeer gradients at 10 and
25 m/s, and rebalance speed. The script saves the analysis MAT-file and
all enabled figures.

Prioritize these plots when reviewing a design screen:

1. **Quality** — invalid-case count and g-g coverage; resolve model failures
   before treating a sensitivity as real.
2. **Sobol total-order sensitivity** — importance including interactions.
3. **Main effects** — direction and approximate size of each parameter's
   effect.
4. **Pareto event/energy plots** — lap-time or points tradeoffs against work.
5. **Balance plots** — low/high-speed understeer gradient and rebalance speed.
6. **Validation** — surrogate prediction quality. Do not rank small effects
   more finely than the model validation error supports.

For an interactive rotatable two-parameter surface, use the returned
analysis structure and any two predictors in the study:

```matlab
f = doeInteractionSurface(analysis,'t_autox','mass','wheelbase');
```

All remaining predictors are held at their valid-design mean. The surface
uses the preferred validated GP or quadratic surrogate; drag the normal
MATLAB 3-D figure to rotate it.

## Recommended study sequence

1. Run 16 cases with the final event/ramp settings as a timing and stability
   pilot. Inspect `sim_log.csv` and invalid-case metrics.
2. Run a 128–192 case broad sensitivity screen. If serial, disable ramps or
   use a smaller study; g-g construction is the main time cost.
3. Review Sobol and validation figures. Retain roughly 6–8 influential,
   independently controllable factors.
4. Run a 128–192 case focused study with full ramp metrics and tighter,
   design-realistic bounds.
5. Only then use `hybrid` or `optimization` mode with handling, energy, or
   rebalance penalties that reflect known targets.

## Common issues

- **`fitrgp` is unavailable:** install/enable Statistics and Machine Learning
  Toolbox; adaptive DOE cannot select later batches without it.
- **Requested worker count has no effect:** create a matching local pool and
  check the local cluster's worker limit. `study.numWorkers` is not a pool
  creation command.
- **Resume mismatch error:** a protected study definition changed. Restore
  the original definition or use a new output directory.
- **Ramp values are NaN:** inspect case validity and
  `state.rampBackfillDiagnostics`; a failed ramp backfill does not overwrite
  otherwise successful lap results.
- **Long serial runtime:** lower `maxCases`, disable ramps for screening, or
  run a separate focused study. Do not silently loosen solver tolerances just
  to shorten a design decision.
