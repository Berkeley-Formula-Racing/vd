# Shared-input regression checks

From the repository root in MATLAB or GNU Octave:

```matlab
addpath(fullfile(pwd,'tests','shared_inputs'));
results = runSharedInputChecks();
```

The fixture was captured from commit
`0b02ea46df58c6f03434c3f0c7021880487804ed` and independently reproduced from
current master `144b3548425d05988668fcc93e3b2b5bd0830676` on 2026-10-15.
Those commits have identical calculation code and tire assets. The capture
normalizes repository paths and represents MATLAB strings as character data.
It includes complete main/ramp defaults and 375 tire states, with both Fx and
Fy outputs from the original Tire2 implementation.

Checks cover:

- Main and complete ramp default equivalence, including the steering limit
  and CAD-to-axle ride-height conversion.
- 185 relocated historical numeric inputs against their original expressions.
- Propagation of temporary master mass/pressure edits; settings-only access
  without the full-car coefficient files. Temporary files and paths are restored.
- Missing, non-finite and wrong-length coefficient validation.
- 375 tire states and an independent nested-interpolation comparison.
- Three Car suspension presets and finite 7-DOF matrices. The vibration
  preset is checked on a struct because upstream Car has no `c` property.
- Nine Car.equations evaluations and a two-second, 101-sample steering response.

GNU Octave 8.4.0 passed these checks. The test-only `octave_compat` adapter
implements the two-dimensional linear griddedInterpolant operation used here.
It is added only for Octave and removed when the check ends; it does not replace
vehicle or tire equations. Do not permanently add the tests tree with genpath.

The MATLAB branch additionally builds both production factories. It has **not**
been executed here. Apps, full lap optimization and Simulink are not covered.
The tests establish refactor regression coverage, not physical model accuracy.
