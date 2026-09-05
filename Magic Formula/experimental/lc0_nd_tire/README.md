# 16x7.5-10 LC0 nondimensional tire prototype

This is an isolated calibration workspace for the Hoosier 43075
16x7.5-10 LC0 on the 8-inch rim. It intentionally does **not** replace
`Tire2`, `carConfig`, or the current lap simulator.

## What is measured

The default configuration reads Round 8 TTC free-rolling runs 15, 16, 18,
and 19. The loader converts USCS data to SI units and rejects nonzero slip
ratio, so these runs are only used to establish the LC0 pure-lateral surface:

- peak lateral friction versus load, pressure, and camber;
- small-slip cornering-stiffness magnitude versus the same conditions;
- coverage diagnostics, which exclude bins that do not contain a meaningful
  slip-angle sweep.

## What is not measured

These runs do not identify LC0 pure longitudinal force or LC0 combined-slip
coupling. `lc0NDCombinedForce` is therefore only a transparent, bounded
friction-ellipse baseline for later comparison with a carefully selected TTC
donor tyre. It must not be calibrated or connected to vehicle predictions as
though it were LC0 combined-slip data.

## Run

From this folder in MATLAB:

```matlab
result = run_lc0_nd_prototype;
```

It writes the exact input configuration, target data, run manifest, summary
table, and figures under `results/`:

- `lc0_nd_target_summary.mat`
- `lc0_nd_peak_mu_vs_load.png` / `.fig`
- `lc0_nd_stiffness_vs_load.png` / `.fig`
- `lc0_nd_data_coverage.png` / `.fig`

Build the provisional donor references separately:

```matlab
donor = run_lc0_nd_donor;
coupling = run_lc0_nd_coupling;
```

The donor source is explicit in `lc0NDConfig`: paired Round 6 LC0 C2000
runs 46 (free rolling) and 47 (drive/brake/combined), both on a 7-inch rim.
They provide a normalized longitudinal reference and a donor coupling-shape
fit only. Their 18x6-10 size and rim mismatch mean their force levels are
never applied directly to the 16x7.5-10/8-inch target.

The coupling artifact and score figure are saved as:

- `lc0_nd_donor_coupling.mat`
- `lc0_nd_coupling_score.png` / `.fig`

## Experimental model and legacy comparison

Run the matched force-level comparison with:

```matlab
comparison = run_lc0_nd_pacejka_comparison;
```

It builds an isolated evaluator with:

- measured Round 8 target `mu_y(alpha)` at 200 lbf, 12 psi, and 0 deg;
- donor `mu_x(kappa)` scaled by the target/donor lateral-peak ratio and the
  explicit `rhoMu` / `rhoStiff` configuration fields;
- the provisional paired-donor coupling exponent (`p = 1.30`).

The comparison invokes the current configured `Tire2` at the exact same
alpha, kappa, normal load, camber, and 12 psi pressure. It writes:

- `lc0_nd_pacejka_comparison.mat` — model inputs and every matched case;
- `lc0_nd_vs_pacejka.png` / `.fig` — pure lateral, pure longitudinal, and
  combined-slip overlays.

Unsupported experimental interpolation points are left as `NaN`; neither
the measured target curve nor the donor curve is silently extrapolated.

Run all isolated checks with:

```matlab
runtests(fullfile(pwd,'tests'))
```

To inventory potential donor data by recorded channels rather than by its
folder name:

```matlab
round6 = fullfile(fileparts(fileparts(pwd)), ...
    'TTC Documentation','Round 6 (18 in Hoosiers)');
catalog = lc0NDCatalogRuns(round6);
catalog(catalog.has_combined_coverage,:)
```

## Next calibration stage

1. Select a TTC donor only after verifying its construction, rim, pressure,
   load, temperature, and slip-ratio coverage.
2. Fit donor longitudinal force in nondimensional coordinates, then use
   bounded `rhoMu` and `rhoStiff` factors at the LC0 reference condition.
3. Compare several bounded combined-slip exponents/curves as an uncertainty
   band, rather than treating an arbitrary donor curve as known LC0 behavior.
4. Validate against logged acceleration, braking, and combined-corner
   vehicle data before moving any calibrated model into `Tire2`.
