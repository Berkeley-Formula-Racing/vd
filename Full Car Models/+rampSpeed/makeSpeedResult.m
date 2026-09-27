function result = makeSpeedResult(task, status, metrics, diagnostics, attempts)
%MAKESPEEDRESULT Package one speed task with an explicit contract status.
if nargin < 2
    error('rampSpeed:invalidSpeedResult', 'A task and status are required.');
end
if nargin < 3 || isempty(metrics)
    metrics = struct();
end
if nargin < 4 || isempty(diagnostics)
    diagnostics = struct();
end
if nargin < 5 || isempty(attempts)
    attempts = struct([]);
end

task = normalizeTask(task);
status = normalizeText(status);
allowedStatuses = ["planned", "running", "converged", "near_feasible", ...
    "infeasible", "solver_failed", "cancelled"];
if ~any(status == allowedStatuses)
    error('rampSpeed:invalidStatus', ...
        'Status must be planned, running, converged, near_feasible, infeasible, solver_failed, or cancelled.');
end
if ~isstruct(metrics) || ~isscalar(metrics)
    error('rampSpeed:invalidMetrics', 'metrics must be a scalar struct.');
end
if ~isstruct(diagnostics) || ~isscalar(diagnostics)
    error('rampSpeed:invalidDiagnostics', 'diagnostics must be a scalar struct.');
end
metrics = normalizeMetricValues(metrics);

result = struct( ...
    'speedIndex', task.speedIndex, ...
    'speed_mps', task.speed_mps, ...
    'origin', task.origin, ...
    'passIndex', task.passIndex, ...
    'status', status, ...
    'metrics', metrics, ...
    'diagnostics', diagnostics, ...
    'attempts', attempts);
end

function task = normalizeTask(task)
if istable(task)
    if height(task) ~= 1 || ~all(ismember( ...
            {'speedIndex', 'speed_mps', 'origin', 'passIndex'}, ...
            task.Properties.VariableNames))
        error('rampSpeed:invalidSpeedTask', ...
            'task must be one row with speedIndex, speed_mps, origin, and passIndex.');
    end
    task = table2struct(task);
elseif ~(isstruct(task) && isscalar(task) && ...
        all(isfield(task, {'speedIndex', 'speed_mps', 'origin', 'passIndex'})))
    error('rampSpeed:invalidSpeedTask', ...
        'task must be a scalar struct or one-row table with speedIndex, speed_mps, origin, and passIndex.');
end

if ~isnumeric(task.speedIndex) || ~isscalar(task.speedIndex) || ...
        ~isfinite(task.speedIndex) || task.speedIndex < 1 || ...
        task.speedIndex ~= fix(task.speedIndex)
    error('rampSpeed:invalidSpeedTask', 'task.speedIndex must be a positive integer.');
end
if ~isnumeric(task.speed_mps) || ~isreal(task.speed_mps) || ...
        ~isscalar(task.speed_mps) || ~isfinite(task.speed_mps) || task.speed_mps <= 0
    error('rampSpeed:invalidSpeedTask', 'task.speed_mps must be finite and positive.');
end
if ~isnumeric(task.passIndex) || ~isscalar(task.passIndex) || ...
        ~isfinite(task.passIndex) || task.passIndex < 1 || ...
        task.passIndex ~= fix(task.passIndex)
    error('rampSpeed:invalidSpeedTask', 'task.passIndex must be a positive integer.');
end

task.speedIndex = double(task.speedIndex);
task.speed_mps = double(task.speed_mps);
task.passIndex = double(task.passIndex);
task.origin = normalizeText(task.origin);
if ~any(task.origin == ["requested", "seed", "refined", "retry"])
    error('rampSpeed:invalidSpeedTask', ...
        'task.origin must be requested, seed, refined, or retry.');
end
end

function metrics = normalizeMetricValues(metrics)
names = fieldnames(metrics);
for index = 1:numel(names)
    name = names{index};
    value = metrics.(name);
    if isnumeric(value) && isempty(value)
        metrics.(name) = NaN;
    elseif isstruct(value) && isscalar(value)
        metrics.(name) = normalizeMetricValues(value);
    end
end
end

function value = normalizeText(value)
if ischar(value) && isrow(value)
    value = string(value);
elseif ~(isstring(value) && isscalar(value) && ~ismissing(value))
    error('rampSpeed:invalidStatus', 'status must be a text scalar.');
end
value = lower(strtrim(value));
if strlength(value) == 0
    error('rampSpeed:invalidStatus', 'status must be a nonempty text scalar.');
end
end
