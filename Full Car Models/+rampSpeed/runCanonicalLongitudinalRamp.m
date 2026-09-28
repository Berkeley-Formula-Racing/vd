function run = runCanonicalLongitudinalRamp(car,settings,caseInfo,callbacks)
%RUNCANONICALLONGITUDINALRAMP Canonical one-pass longitudinal study core.

if nargin < 2 || isempty(settings), settings = struct(); end
if nargin < 3 || isempty(caseInfo), caseInfo = struct(); end
if nargin < 4 || isempty(callbacks), callbacks = struct(); end
if ~isstruct(settings) || ~isscalar(settings)
    error("rampSpeed:invalidSettings","settings must be a scalar struct.");
end
if ~isstruct(caseInfo) || ~isscalar(caseInfo)
    error("rampSpeed:invalidCaseInfo","caseInfo must be a scalar struct.");
end
if ~isstruct(callbacks) || ~isscalar(callbacks)
    error("rampSpeed:invalidCallbacks","callbacks must be a scalar struct.");
end

[profile,settings] = rampSpeed.resolveSolverProfileFromSettings(settings);
car = rampSpeed.applySolverProfileToCar(car,profile);
grid = speedGrid(settings);
requested = requestedSpeeds(settings);
[requested,speedGridMeta] = rampSpeed.resolveFixedSpeedGrid(requested,grid);
rampModel = rampSpeed.buildRampModel(car,setupForModel(caseInfo), ...
    profile,requested(isfinite(requested) & requested > 0));

request = struct("rampType","longitudinal","settings",settings, ...
    "solverProfile",profile.id,"setupKey",caseKey(caseInfo));
request.settings.speeds_mps = requested;
request.settings.speedPolicy = "fixed";
control = makeControl(callbacks,rampModel);

[tasks,~] = makeTasks(requested,ones(numel(requested),1), ...
    speedGridMeta.provenance.source, ...
    speedGridMeta.provenance.reason,1);
[allResults,execution] = rampSpeed.executeSpeedPlan( ...
    request,tasks,@solveTask,control);
allTasks = tasks;
finalSpeeds = requested;
requestedForRun = requested;

settingsForRun = settings;
settingsForRun.speeds = finalSpeeds.';
settingsForRun.solverProfile = profile.id;
settingsForRun.solverOptions = profile.solverOptions;
settingsForRun.residualTolerance = profile.solverOptions.constraintTolerance;
settingsForRun.ceqTol = profile.solverOptions.constraintTolerance;
settingsForRun.constraintTolerance = profile.solverOptions.constraintTolerance;
settingsForRun.inequalityTolerance = profile.solverOptions.constraintTolerance;
request.settings = settingsForRun;
raw = rampSpeed.assembleRun(allTasks,allResults,request,execution);
runMeta = struct("source","rampSpeed.runLongitudinalRamp", ...
    "status",raw.status,"started",datetime('now'), ...
    "completed",datetime('now'),"solver",profile.solverOptions, ...
    "solverProfile",rampSpeed.serializeSolverProfile(profile), ...
    "requestedSpeeds_mps",requestedForRun, ...
    "lateralMetricsApplicable",false, ...
    "speedGrid",speedGridMeta, ...
    "rampModel",rampModel, ...
    "retrySpeeds_mps",speedGridMeta.exactRetrySpeeds_mps, ...
    "speedErrors",raw.raw.speedErrors, ...
    "warnings",strings(0,1),"errors",strings(0,1));
if any(raw.perSpeed.status == "solver_failed")
    runMeta.warnings(end+1,1) = "One or more speeds returned solver_failed.";
end
if any(raw.perSpeed.status == "infeasible")
    runMeta.warnings(end+1,1) = "One or more speeds were physically infeasible.";
end
run = rampSpeed.normalizeRampResult(raw,"longitudinal",settingsForRun, ...
    caseInfo,runMeta);
run.settings = settingsForRun;
run.settings.speeds = finalSpeeds.';
run.runMeta.speedGrid = speedGridMeta;
run.runMeta.retrySpeeds_mps = speedGridMeta.exactRetrySpeeds_mps;
run.runMeta.requestedSpeeds_mps = requestedForRun;
run.runMeta.speedErrors = raw.raw.speedErrors;
run.runMeta.solverProfile = rampSpeed.serializeSolverProfile(profile);
run.runMeta.solver = profile.solverOptions;
run.raw = raw;
run.raw.settings = run.settings;
run.raw.status = raw.status;
run.status = aggregateRunStatus(raw,execution);
run.runMeta.status = run.status;
run.perSpeed = addStableFields(run.perSpeed,raw.perSpeed,finalSpeeds);
invalidRequested = ~isfinite(run.perSpeed.speed_mps);
if any(invalidRequested)
    run.perSpeed.valid(invalidRequested) = false;
    run.perSpeed.status(invalidRequested) = "invalid";
    run.perSpeed.reason(invalidRequested) = run.perSpeed.solver_reason(invalidRequested);
    run.status = aggregateRunStatus(run.raw,execution);
    run.runMeta.status = run.status;
end

function result = solveTask(task,seed,~,controlArg)
result = rampSpeed.solveLongitudinalPoint(car,task.speed_mps,seed, ...
    profile,controlArg);
result.speedIndex = task.speedIndex;
result.speed_mps = task.speed_mps;
result.origin = string(task.origin);
result.passIndex = task.passIndex;
end
end

function [tasks,nextId] = makeTasks(speeds,passes,origins,reasons,nextId)
speeds = double(speeds(:));
n = numel(speeds);
if numel(passes) ~= n || numel(origins) ~= n || numel(reasons) ~= n
    error("rampSpeed:invalidSpeedTasks","Task metadata dimensions do not match.");
end
tasks = table((nextId:nextId+n-1).',speeds,string(origins(:)), ...
    double(passes(:)),string(reasons(:)),repmat("planned",n,1), ...
    'VariableNames',{'speedIndex','speed_mps','origin','passIndex', ...
    'refinementReason','status'});
nextId = nextId+n;
end

function control = makeControl(callbacks,rampModel)
control = struct();
if nargin >= 2 && isstruct(rampModel) && isscalar(rampModel)
    control.rampModel = rampModel;
end
if isfield(callbacks,"onProgress") && ~isempty(callbacks.onProgress)
    control.onTaskStart = @forwardTaskStart;
    control.onProgress = @forwardProgress;
end
if isfield(callbacks,"isCancelled") && ~isempty(callbacks.isCancelled)
    control.shouldCancel = callbacks.isCancelled;
end

    completedSpeeds = 0;

    function forwardTaskStart(event)
        callbacks.onProgress(legacyProgressEvent(event,completedSpeeds));
    end

    function forwardProgress(event)
        % The old app callback contract emitted no completion event for a
        % non-finite planned row. Preserve that behavior at this adapter seam;
        % executeSpeedPlan still records the full event internally.
        if ~isfinite(event.speed_mps) && string(event.status) == "solver_failed"
            return
        end
        if ismember(string(event.status),["converged","near_feasible"])
            completedSpeeds = completedSpeeds + 1;
        end
        callbacks.onProgress(legacyProgressEvent(event,completedSpeeds));
    end

    function legacy = legacyProgressEvent(event,completed)
        legacy = struct("phase","speed", ...
            "speedIndex",event.speedIndex, ...
            "speed_mps",event.speed_mps, ...
            "completedSpeeds",completed, ...
            "requestedSpeeds",event.plannedTasks);
    end
end

function setup = setupForModel(caseInfo)
setup = struct();
if isstruct(caseInfo) && isscalar(caseInfo) && isfield(caseInfo,"setupSpec") && ...
        isstruct(caseInfo.setupSpec) && isscalar(caseInfo.setupSpec)
    setup = caseInfo.setupSpec;
end
end

function run = addStableFields(run,raw,finalSpeeds)
rawSpeeds = double(raw.speed_mps(:));
loc = zeros(numel(finalSpeeds),1);
for i = 1:numel(finalSpeeds)
    if isfinite(finalSpeeds(i))
        finiteRaw = rawSpeeds(isfinite(rawSpeeds));
        scale = max([1;abs(finalSpeeds(i));abs(finiteRaw)]);
        match = find(abs(rawSpeeds-finalSpeeds(i)) <= 32*eps(scale),1);
        if ~isempty(match), loc(i) = match; end
    elseif i <= numel(rawSpeeds)
        loc(i) = i;
    end
    if isempty(loc(i)) || loc(i) == 0, loc(i) = min(i,numel(rawSpeeds)); end
end
loc(loc==0) = 1;
run.speedIndex = raw.speedIndex(loc);
run.origin = raw.origin(loc);
run.passIndex = raw.passIndex(loc);
run.refinementReason = raw.refinementReason(loc);
run.solver_status = raw.status(loc);
run.solver_reason = raw.reason(loc);
end

function speeds = requestedSpeeds(settings)
speeds = zeros(0,1);
if isfield(settings,"speeds") && ~isempty(settings.speeds)
    speeds = double(settings.speeds(:));
end
end

function grid = speedGrid(settings)
grid = struct("mode","fixed");
if isfield(settings,"speedGrid") && ~isempty(settings.speedGrid)
    if isstruct(settings.speedGrid), grid = settings.speedGrid;
    else, grid.mode = settings.speedGrid; end
end
if ~isstruct(grid) || ~isscalar(grid)
    error("rampSpeed:invalidFixedSpeedGrid", ...
        "speedGrid must be a scalar struct or mode.");
end
end

function key = caseKey(caseInfo)
key = "default";
if isstruct(caseInfo) && isfield(caseInfo,"id") && ~isempty(caseInfo.id)
    key = string(caseInfo.id);
end
end

function value = aggregateRunStatus(raw,execution)
if isfield(execution,"status") && string(execution.status) == "cancelled"
    value = "cancelled";
elseif isempty(raw.perSpeed) || ~any(raw.perSpeed.valid)
    value = "failed";
elseif any(~raw.perSpeed.valid)
    value = "partial";
else
    value = "completed";
end
end
