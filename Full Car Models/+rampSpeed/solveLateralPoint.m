function result = solveLateralPoint(car,task,settings,profile,control)
%SOLVELATERALPOINT Solve one lateral ramp speed into a SpeedResult.
%
% The legacy rampSweep remains the numerical lateral solver. This boundary
% gives it one speed at a time, retains its detailed point table in diagnostics,
% and converts all failure modes into an explicit SpeedResult status.

if nargin < 3 || isempty(settings), settings = struct(); end
if nargin < 4 || isempty(profile), profile = "accurate"; end
if nargin < 5 || isempty(control), control = struct(); end
if ~isstruct(task) || ~isscalar(task) || ...
        ~all(isfield(task,{'speedIndex','speed_mps','origin','passIndex'}))
    error("rampSpeed:invalidSpeedTask", ...
        "task must contain speedIndex, speed_mps, origin, and passIndex.");
end
if ~isstruct(settings) || ~isscalar(settings)
    error("rampSpeed:invalidSettings","settings must be a scalar struct.");
end
if ~isstruct(control) || ~isscalar(control)
    error("rampSpeed:invalidControl","control must be a scalar struct.");
end
if ~isfinite(task.speed_mps) || task.speed_mps <= 0
    error("rampSpeed:invalidSpeedDomain", ...
        "A lateral point speed must be finite and positive.");
end

profile = rampSpeed.resolveSolverProfile(profile);
if shouldCancel(control)
    result = rampSpeed.makeSpeedResult(task,"cancelled",struct(), ...
        struct("reason","cancelled before lateral point evaluation", ...
        "errorIdentifier","rampSpeed:cancelled"));
    return
end

options = settings;
options.speeds = double(task.speed_mps);
options.verbose = false;
options.progressFcn = [];
options.cancelFcn = @cancelRequested;
options.solverOptions = profile.solverOptions;

try
    raw = rampSweep(car,options);
catch ME
    status = "solver_failed";
    if strcmp(ME.identifier,"rampSweep:cancelled") || ...
            strcmp(ME.identifier,"rampSpeed:cancelled")
        status = "cancelled";
    elseif any(strcmp(ME.identifier,["rampSweep:noSolution", ...
            "rampSpeed:noSolution"]))
        status = "infeasible";
    end
    diagnostics = struct("reason",string(ME.message), ...
        "errorIdentifier",string(ME.identifier),"raw",failureRaw(task,ME), ...
        "pointCount",0);
    result = rampSpeed.makeSpeedResult(task,status,struct(),diagnostics);
    return
end

if ~isfield(raw,"perSpeed") || isempty(raw.perSpeed) || ...
        (istable(raw.perSpeed) && height(raw.perSpeed) == 0)
    result = rampSpeed.makeSpeedResult(task,"infeasible",struct(), ...
        struct("reason","no lateral ramp row was returned", ...
        "errorIdentifier","rampSpeed:noSolution","raw",raw, ...
        "pointCount",0));
    return
end

row = raw.perSpeed(1,:);
metrics = tableRowMetrics(row);
diagnostics = struct("reason","lateral ramp point solved", ...
    "errorIdentifier","","raw",raw,"pointCount",height(raw.points), ...
    "rampStatus",string(raw.status));
attempt = struct("source","rampSweep","status","converged", ...
    "exitflag",selectedExitflag(row),"reason","lateral ramp returned a row");
result = rampSpeed.makeSpeedResult(task,"converged",metrics,diagnostics,attempt);

    function value = cancelRequested()
        value = shouldCancel(control);
    end
end

function metrics = tableRowMetrics(row)
metrics = struct();
names = row.Properties.VariableNames;
for i = 1:numel(names)
    name = names{i};
    value = row.(name);
    if isnumeric(value) && isscalar(value)
        metrics.(name) = double(value);
    elseif islogical(value) && isscalar(value)
        metrics.(name) = logical(value);
    elseif isstring(value) && isscalar(value)
        metrics.(name) = value;
    end
end
end

function value = selectedExitflag(row)
value = NaN;
n = {"n_exit1","n_exit2","exitflag"};
for i = 1:numel(n)
    if ismember(n{i},row.Properties.VariableNames)
        candidate = row.(n{i});
        if isnumeric(candidate) && isscalar(candidate) && isfinite(candidate)
            value = double(candidate);
            return
        end
    end
end
end

function raw = failureRaw(task,ME)
raw = struct("status","failed","settings",struct("speeds",task.speed_mps), ...
    "perSpeed",table(),"points",table(),"error",struct( ...
    "identifier",string(ME.identifier),"message",string(ME.message)));
end

function value = shouldCancel(control)
value = false;
if isfield(control,"shouldCancel") && ~isempty(control.shouldCancel)
    callback = control.shouldCancel;
    if isa(callback,"function_handle")
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
