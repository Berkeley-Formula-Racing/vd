# QSS telemetry HDF5 interchange format, version 1.0

This contract is owned by the telemetry-viewer integration task. MATLAB writes it; Python reads it. All numeric values use SI units except rotational speed (`rpm`) and dimensionless quantities such as gear (`1`). Arrays are `float64`; validity masks are `uint8`; JSON is UTF-8 encoded bytes in a one-dimensional `uint8` dataset.

## File tree

```text
/manifest_json
/tracks/<track_id>/metadata_json
/tracks/<track_id>/distance_m
/tracks/<track_id>/curvature_per_m
/tracks/<track_id>/x_m                 # optional, only with y_m
/tracks/<track_id>/y_m                 # optional, only with x_m
/cases/<case_id>/metadata_json
/cases/<case_id>/setup_json
/cases/<case_id>/envelope/...           # optional self-describing channels
/cases/<case_id>/laps/<lap_id>/metadata_json
/cases/<case_id>/laps/<lap_id>/axes/<axis_id>/time_s
/cases/<case_id>/laps/<lap_id>/axes/<axis_id>/distance_m
/cases/<case_id>/laps/<lap_id>/channels/<channel_id>/values
/cases/<case_id>/laps/<lap_id>/channels/<channel_id>/valid
/cases/<case_id>/laps/<lap_id>/channels/<channel_id>/metadata_json
```

`manifest_json` must include `format: "qss-telemetry"`, `schema_version: "1.0"`, `run_uuid`, `creation_utc`, `source_type`, producer and MATLAB versions, Git revision/dirty state, and a case/track inventory. Writer implementations may add fields.

Track distance, every axis time, and every axis distance are finite and monotonically non-decreasing. `x_m` and `y_m` either both exist or neither exists. When no geometry is present, the viewer may integrate curvature but must label the result schematic.

Each channel metadata JSON object has exactly these required fields (additional fields are allowed):

```json
{
  "id": "speed_mps",
  "label": "Vehicle speed",
  "unit": "m/s",
  "axis_id": "native",
  "origin": "simulation_output",
  "interpolation": "linear",
  "description": "",
  "coordinate_frame": "vehicle",
  "sign_convention": ""
}
```

`origin` is one of `simulation_output`, `derived`, `qss_reconstructed`, or `measured`. `interpolation` is `linear`, `previous`, or `none`. Invalid samples use `valid == 0` and `values = NaN`; zeros never stand for missing data. A channel is absent when it was not produced.

Files are written atomically: write a temporary sibling, validate the completed content, then rename. Exporters default to unique names and require explicit overwrite.
