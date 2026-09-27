function [results,execution] = executeSpeedPlan(request,tasks,pointSolver,control)
%EXECUTESPEEDPLAN Evaluate each stable speed task exactly once.
if nargin < 4 || isempty(control), control = struct(); end
if nargin < 3 || ~isa(pointSolver,'function_handle')
    error("rampSpeed:invalidPointSolver", ...
        "pointSolver must be a function handle.");
end
if ~istable(tasks)
    error("rampSpeed:invalidSpeedTasks","tasks must be a table.");
end
required = {'speedIndex','speed_mps','origin','passIndex'};
if ~all(ismember(required,tasks.Properties.VariableNames))
    error("rampSpeed:invalidSpeedTasks", ...
        "tasks must include speedIndex, speed_mps, origin, and passIndex.");
end
indices = double(tasks.speedIndex(:));
if any(~isfinite(indices)) || any(indices < 1) || ...
        any(indices ~= fix(indices)) || numel(unique(indices)) ~= numel(indices)
    error("rampSpeed:duplicateSpeedIndex", ...
        "Each execution plan must contain one unique positive speedIndex.");
end
speeds = double(tasks.speed_mps(:));
allowInvalidSpeedRows = isfield(request,"allowInvalidSpeedRows") && ...
    logical(request.allowInvalidSpeedRows);
if any(isfinite(speeds) & speeds <= 0) && ~allowInvalidSpeedRows
    error("rampSpeed:invalidSpeedTasks", ...
        "Finite planned speeds must be positive.");
end
if ~isstruct(control) || ~isscalar(control)
    error("rampSpeed:invalidControl","control must be a scalar struct.");
end

n = height(tasks);
results = cell(n,1);
progressEvents = repmat(progressTemplate(),n,1);
execution = struct("status","running","tasks",tasks, ...
    "results",{results},"progressEvents",progressEvents, ...
    "completedTasks",0,"plannedTasks",n,"evaluatedSpeedIndices",zeros(0,1));
previousState = [];
cancelled = false;

for i = 1:n
    task = table2struct(tasks(i,:));
    if shouldCancel(control)
        cancelled = true;
        [results,progressEvents] = fillCancelled(results,progressEvents, ...
            tasks,i,"cancelled before this task was evaluated");
        break
    end

    startEvent = progressTemplate();
    startEvent.speedIndex = task.speedIndex;
    startEvent.speed_mps = task.speed_mps;
    startEvent.completedTasks = execution.completedTasks;
    startEvent.plannedTasks = n;
    startEvent.status = "running";
    notifyTaskStart(control,startEvent);

    invalidSpeed = ~isfinite(task.speed_mps) || ...
        (allowInvalidSpeedRows && task.speed_mps <= 0);
    if invalidSpeed
        pointResult = failedPointResult(task,"solver_failed", ...
            "failed: planned speed is not finite and positive");
        results{i} = pointResult;
        event = progressTemplate();
        event.speedIndex = task.speedIndex;
        event.speed_mps = task.speed_mps;
        event.completedTasks = execution.completedTasks;
        event.plannedTasks = n;
        event.status = string(pointResult.status);
        progressEvents(i) = event;
        notifyProgress(control,event);
        if shouldCancel(control)
            cancelled = true;
            [results,progressEvents] = fillCancelled(results,progressEvents, ...
                tasks,i+1,"cancelled before this task was evaluated");
            break
        end
        continue
    end

    try
        pointResult = pointSolver(task,previousState,request,control);
        pointResult = validatePointResult(pointResult,task);
        if isfield(pointResult,'state') && isnumeric(pointResult.state) && ...
                numel(pointResult.state) == 9 && all(isfinite(pointResult.state(:))) && ...
                isValidResult(pointResult)
            previousState = double(pointResult.state(:).');
        end
    catch ME
        if strcmp(ME.identifier,"rampSpeed:cancelled")
            pointResult = cancelledPointResult(task, ...
                "cancelled during this task evaluation");
            cancelled = true;
        else
            pointResult = failedPointResult(task,"solver_failed", ...
                string(ME.identifier)+": "+string(ME.message));
        end
    end
    results{i} = pointResult;
    execution.completedTasks = i;
    execution.evaluatedSpeedIndices(end+1,1) = task.speedIndex;
    event = progressTemplate();
    event.speedIndex = task.speedIndex;
    event.speed_mps = task.speed_mps;
    event.completedTasks = i;
    event.plannedTasks = n;
    event.status = string(pointResult.status);
    progressEvents(i) = event;
    notifyProgress(control,event);
    if cancelled || shouldCancel(control)
        cancelled = true;
        [results,progressEvents] = fillCancelled(results,progressEvents, ...
            tasks,i+1,"cancelled before this task was evaluated");
        break
    end
end

if cancelled
    execution.status = "cancelled";
else
    execution.status = "completed";
end
execution.results = results;
execution.progressEvents = progressEvents;
execution.cancelled = cancelled;
end

function result = validatePointResult(result,task)
if ~isstruct(result) || ~isscalar(result) || ...
        ~all(isfield(result,{'speedIndex','speed_mps','status'}))
    error("rampSpeed:invalidPointResult", ...
        "The point solver must return a scalar SpeedResult.");
end
if double(result.speedIndex) ~= double(task.speedIndex) || ...
        abs(double(result.speed_mps)-double(task.speed_mps)) > ...
        32*eps(max([1,double(task.speed_mps)]))
    error("rampSpeed:resultTaskMismatch", ...
        "Point solver returned a result for a different stable task.");
end
if ~isfield(result,'metrics'), result.metrics = struct(); end
if ~isfield(result,'diagnostics'), result.diagnostics = struct(); end
if ~isfield(result,'attempts'), result.attempts = struct([]); end
end

function result = failedPointResult(task,status,reason)
diagnostics = struct("reason",string(reason),"errorIdentifier", ...
    "rampSpeed:"+string(status));
diagnostics.reason = "rampSpeed:"+string(status)+": "+diagnostics.reason;
if isfinite(task.speed_mps) && task.speed_mps > 0
    result = rampSpeed.makeSpeedResult(task,status,struct(),diagnostics);
else
    result = struct("speedIndex",double(task.speedIndex), ...
        "speed_mps",double(task.speed_mps),"origin",string(task.origin), ...
        "passIndex",double(task.passIndex),"status",string(status), ...
        "metrics",struct(),"diagnostics",diagnostics,"attempts",struct([]));
end
end

function result = cancelledPointResult(task,reason)
result = rampSpeed.makeSpeedResult(task,"cancelled",struct(), ...
    struct("reason",string(reason),"errorIdentifier","rampSpeed:cancelled"));
end

function [results,events] = fillCancelled(results,events,tasks,startIndex,reason)
for i = startIndex:height(tasks)
    task = table2struct(tasks(i,:));
    results{i} = cancelledPointResult(task,reason);
    event = progressTemplate();
    event.speedIndex = task.speedIndex;
    event.speed_mps = task.speed_mps;
    event.completedTasks = i-1;
    event.plannedTasks = height(tasks);
    event.status = "cancelled";
    events(i) = event;
end
end

function event = progressTemplate()
event = struct("speedIndex",NaN,"speed_mps",NaN, ...
    "completedTasks",0,"plannedTasks",0,"status","");
end

function notifyProgress(control,event)
if ~isfield(control,'onProgress') || isempty(control.onProgress)
    return
end
try
    control.onProgress(event);
catch
    % Progress callbacks are advisory and cannot invalidate a result.
end
end

function notifyTaskStart(control,event)
if ~isfield(control,'onTaskStart') || isempty(control.onTaskStart)
    return
end
try
    control.onTaskStart(event);
catch
    % Progress callbacks are advisory and cannot invalidate a result.
end
end

function value = shouldCancel(control)
value = false;
if isfield(control,'shouldCancel') && ~isempty(control.shouldCancel)
    callback = control.shouldCancel;
    if isa(callback,'function_handle')
        value = callback();
    else
        value = callback;
    end
    if ~isscalar(value)
        error("rampSpeed:invalidCancellation", ...
            "control.shouldCancel must return a scalar logical.");
    end
    value = logical(value);
end
end

function tf = isValidResult(result)
tf = ismember(string(result.status),["converged","near_feasible"]);
end
