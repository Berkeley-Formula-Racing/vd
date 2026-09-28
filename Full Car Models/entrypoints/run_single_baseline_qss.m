%% run_single_baseline_qss
% Run one scalar baseline QSS autocross simulation and export its telemetry.
%
% From MATLAB, run:
%   run('C:\VD\.worktrees\telemetry-viewer\Full Car Models\entrypoints\run_single_baseline_qss.m')
%
% The result is written to the repository/worktree root. Existing output is preserved;
% exportQSSResults creates a unique filename when needed.

clear classes
entrypointDir = fileparts(mfilename('fullpath'));
run(fullfile(entrypointDir, 'bootstrap.m'));

eventName = "autocross";
outputPath = fullfile(repoRoot, 'baseline.qss.h5');
includeDetail = false;

fprintf('Building one scalar baseline car...\n');
[allCars, eventParams, designTable] = carConfig();
baselineMask = designTable.static_front_ride_height_in == 0 & ...
    designTable.static_rear_ride_height_in == 0;
baselineIndices = find(baselineMask);
if numel(baselineIndices) ~= 1
    error('run_single_baseline_qss:notScalarBaseline', ...
        ['Could not identify exactly one zero-ride-height baseline in the ' ...
         'carConfig() grid.']);
end
baselineIndex = baselineIndices(1);
baselineCell = allCars(baselineIndex, :);

car = baselineCell{1, 1};
accelCar = baselineCell{1, 2};
fprintf('Solving baseline g-g envelope...\n');
car = makeGG(gg2(car, 0), car);
ev = Events2(car, accelCar, eventParams);

solverName = char(eventName);
solverName(1) = upper(solverName(1));
fprintf('Running baseline %s...\n', eventName);
ev.(solverName)();

resultOptions = struct( ...
    'reconstruct_detail', includeDetail, ...
    'car', car, ...
    'case_id', 'baseline', ...
    'case_label', 'Baseline');
result = buildTelemetryResult(ev, resultOptions);
outputPath = exportQSSResults(result, outputPath);

fprintf('Wrote baseline QSS telemetry: %s\n', outputPath);
