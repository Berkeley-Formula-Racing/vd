# Full Car Models entrypoints

The runnable MATLAB scripts live in this folder. Each script bootstraps the
full model path and changes MATLAB's working folder to the repository root, so
it can be launched from any current folder.

To run one scalar baseline QSS autocross simulation and export telemetry:

```matlab
run('C:\VD\.worktrees\telemetry-viewer\Full Car Models\entrypoints\run_single_baseline_qss.m')
```

The output is `baseline.qss.h5` in the repository/worktree root (currently
`C:\VD\.worktrees\telemetry-viewer\baseline.qss.h5`), or a numbered sibling if
that file already exists. To open it in the Python viewer:

```powershell
cd C:\VD\.worktrees\telemetry-viewer\telemetry_viewer
uv run qss-telemetry-viewer C:\VD\.worktrees\telemetry-viewer\baseline.qss.h5
```
