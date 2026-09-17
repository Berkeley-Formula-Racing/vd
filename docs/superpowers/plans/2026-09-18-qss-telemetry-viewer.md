# QSS telemetry viewer implementation plan

Create a versioned HDF5 interchange between the MATLAB QSS solver and a local PySide6/PyQtGraph desktop viewer. The v1 UI opens one or more simulation result files, compares laps by physical distance or time with manual offsets, shows linked plots and an explicitly schematic map when only curvature is available, playback, setup differences, quality information, CSV export, and screenshots. Future real telemetry must adapt into the same Python `LapData` interface; importing it is outside v1.

The immutable shared interface is `docs/qss-telemetry-format.md` and `telemetry_viewer/qss_telemetry/schema.py`. The demo fixture in `qss_telemetry/fixtures.py` is the cross-team development fixture. Preserve `python/` and `.claude/` completely.

## Work packages

1. MATLAB only: add additive result capture/export and documented QSS reconstruction in `Full Car Models/utilities/` plus MATLAB tests. Reconstruct at approximately 2 m and significant transitions using the QSS state `[steering, signed control demand, longitudinal velocity, lateral velocity, yaw rate, kappaFL, kappaFR, kappaRL, kappaRR]`; preserve solver bounds; independently store residuals/status/eval-count/failure reason; enforce residual limits of 1e-3 m/s^2 acceleration, 1e-3 rad/s^2 yaw, 1e-2 N m torque, virtual load >= -0.1 N, aero residual 1e-7 inches, and max 1500 evaluations. Do not reconstruct unsupported launch/shift states.
2. Python core only: HDF5 reader, format validation, `load_results`, `load_lap`, distance/time alignment, and delta time. Continuous channels are linearly interpolated; discrete channels use previous sample; preserve gaps and only compare overlap; positive delta means comparison slower.
3. UI only: PySide6/PyQtGraph app built against the shared dataclasses/fixture, linked plots, playback, track map cursor, case/lap browser, channel selection, delta display, setup/quality panes, and UI tests.

Integration validates the MATLAB-written artifact in Python, checks corruption/schema rejection and synthetic delta-time cases, runs UI tests where dependencies exist, and reports the known headless MATLAB startup failure rather than claiming MATLAB verification.
