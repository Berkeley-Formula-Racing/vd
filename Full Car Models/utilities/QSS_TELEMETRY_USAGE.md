# QSS telemetry export

The telemetry adapter is additive. It reads a solved `Events2` object or a
plain MATLAB result struct and writes the version 1 HDF5 interchange described
in `docs/qss-telemetry-format.md`:

```matlab
addpath('Full Car Models/utilities')
result = buildTelemetryResult(events2Object);
path = exportQSSResults(result,'results/qss_run.h5');
```

`exportQSSResults` uses a unique sibling name when the requested file already
exists. Set `struct('overwrite',true)` only when replacing a file is intended.
The writer creates a sibling temporary HDF5 file, validates it, and then moves
it into place. No event method is run by either capture or export.

The supported plain source form is a scalar struct with a `track` (or
`tracks`) and `laps` field. A track has `id`, `distance_m`, and
`curvature_per_m`. Each lap has `id`, `track_id`, `time_s` and either a
`channels` struct array or well-known vectors such as `speed_mps`,
`long_accel_mps2`, `lat_accel_mps2`, `gear`, and `steering_angle_rad`. A channel
has `id`, `values`, `label`, `unit`, `origin`, and `interpolation`; the other
metadata fields default to the v1 vehicle conventions. A channel whose vector
length differs from the native time/profile vector receives its own axis and
is retained in full. `lap.metadata.node_segment_mapping` records each source
field, length, axis, node indices, and segment indices.

Native `Events2` autocross, endurance, and acceleration outputs are captured
when their output structs are populated. Endurance is represented by one
steady lap and records its lap count, raw representative lap time, and the
driver-adjusted event time separately. Lap metadata records the raw time field
and whether time was already cumulative.

Optional detail reconstruction is enabled with
`struct('reconstruct_detail',true)`, or can be called directly:

```matlab
detail = reconstructLapTelemetry(events2Object.car,result.cases(1).laps(1));
```

Detail points include approximately every two metres, start/end, curvature
extrema, and detected gear, propulsion, or braking transitions. The QSS state
is `[steering_rad, signed_control_demand, v_long, v_lat, yaw_rate,
kappa_FL, kappa_FR, kappa_RL, kappa_RR]`. Launch and shift points are marked
`unavailable`; solver failures keep NaN values and a zero validity mask. Every
point carries status, residuals, evaluation count, and failure reason in
`lap.diagnostics` and `/diagnostics_json`.

`validateQSSResults` returns `[ok,report]`; with no output arguments it raises
on invalid content. MATLAB startup in the current headless environment may
fail before tests can be discovered; the integration report records that
result rather than treating it as a passing MATLAB run.
