function run = runCanonicalLateralRamp(car,settings,caseInfo,callbacks)
%RUNCANONICALLATERALRAMP Execute lateral speeds through the shared task core.

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
if ~isfield(settings,"mode") || isempty(settings.mode)
    settings.mode = "balanced";
end
requested = requestedSpeeds(settings);
if isempty(requested), requested = (5:2.5:30).'; end
requested = unique(double(requested(:)),"stable");

request = struct("rampType","lateral","settings",settings, ...
    "solverProfile",profile.id,"setupKey",caseKey(caseInfo), ...
    "allowInvalidSpeedRows",true);
request.settings.speeds_mps = requested;
request.settings.speedPolicy = "fixed";
control = makeControl(callbacks);

[tasks,~] = makeTasks(requested,ones(numel(requested),1), ...
    repmat("requested",numel(requested),1), ...
    repmat("",numel(requested),1),1);
[results,execution] = rampSpeed.executeSpeedPlan(request,tasks, ...
    @solveTask,control);

settingsForRun = settings;
settingsForRun.speeds = requested.';
settingsForRun.solverProfile = profile.id;
settingsForRun.solverOptions = profile.solverOptions;
settingsForRun.residualTolerance = profile.solverOptions.constraintTolerance;
settingsForRun.ceqTol = profile.solverOptions.constraintTolerance;
settingsForRun.constraintTolerance = profile.solverOptions.constraintTolerance;
settingsForRun.inequalityTolerance = profile.solverOptions.constraintTolerance;
request.settings = settingsForRun;
raw = rampSpeed.assembleRun(tasks,results,request,execution);
runMeta = struct("source","rampSpeed.runLateralRamp", ...
    "status",raw.status,"started",datetime('now'), ...
    "completed",datetime('now'),"solver",profile.solverOptions, ...
    "solverProfile",rampSpeed.serializeSolverProfile(profile), ...
    "requestedSpeeds_mps",requested,"lateralMetricsApplicable",true, ...
    "speedErrors",raw.raw.speedErrors,"warnings",strings(0,1), ...
    "errors",strings(0,1),"speedGrid",fixedGrid(requested));
if profile.approximate
    runMeta.warnings(end+1,1) = ...
        "Approximate aero preview: ride-height aero iteration is disabled.";
end
if any(raw.perSpeed.status == "solver_failed")
    runMeta.warnings(end+1,1) = "One or more speeds returned solver_failed.";
end
if any(raw.perSpeed.status == "infeasible")
    runMeta.warnings(end+1,1) = "One or more speeds were physically infeasible.";
end

run = rampSpeed.normalizeRampResult(raw,"lateral",settingsForRun, ...
    caseInfo,runMeta);
run.settings = settingsForRun;
run.settings.speeds = requested.';
run.runMeta.requestedSpeeds_mps = requested;
run.runMeta.speedGrid = runMeta.speedGrid;
run.runMeta.speedErrors = raw.raw.speedErrors;
run.runMeta.solverProfile = rampSpeed.serializeSolverProfile(profile);
run.runMeta.solver = profile.solverOptions;
run.raw = raw;
run.raw.settings = run.settings;
run.raw.status = raw.status;
run.status = aggregateRunStatus(raw,execution);
run.runMeta.status = run.status;
run.perSpeed = addStableFields(run.perSpeed,raw.perSpeed,requested);
% Keep the historical lateral CaseRun status vocabulary at the normalized
% boundary; the canonical raw table still distinguishes solver_failed,
% infeasible, and cancelled rows explicitly.
run.perSpeed.status(run.perSpeed.status == "invalid") = "failed";

    function result = solveTask(task,~,~,controlArg)
        result = rampSpeed.solveLateralPoint(car,task,settings,profile,controlArg);
    end
end

function [tasks,nextId] = makeTasks(speeds,passes,origins,reasons,nextId)
n = numel(speeds);
tasks = table((nextId:nextId+n-1).',double(speeds(:)),string(origins(:)), ...
    double(passes(:)),string(reasons(:)),repmat("planned",n,1), ...
    'VariableNames',{'speedIndex','speed_mps','origin','passIndex', ...
    'refinementReason','status'});
nextId = nextId+n;
end

function control = makeControl(callbacks)
control = struct();
if isfield(callbacks,"onProgress") && ~isempty(callbacks.onProgress)
    control.onProgress = @forwardProgress;
end
if isfield(callbacks,"isCancelled") && ~isempty(callbacks.isCancelled)
    control.shouldCancel = callbacks.isCancelled;
end

    function forwardProgress(event)
        callbacks.onProgress(struct("phase","speed", ...
            "speedIndex",event.speedIndex,"speed_mps",event.speed_mps, ...
            "completedSpeeds",sum(event.completedTasks >= 1), ...
            "requestedSpeeds",event.plannedTasks));
    end
end

function grid = fixedGrid(requested)
grid = struct("mode","fixed","requestedSpeeds_mps",requested, ...
    "seedSpeeds_mps",requested,"finalSpeeds_mps",requested,"passes",0, ...
    "provenance",table(requested,repmat("requested",numel(requested),1), ...
    zeros(numel(requested),1),repmat("user_requested",numel(requested),1), ...
    'VariableNames',{'speed_mps','source','pass','reason'}), ...
    "refinementHistory",struct([]),"stopReason","", ...
    "exactRetrySpeeds_mps",zeros(0,1));
end

function speeds = requestedSpeeds(settings)
speeds = zeros(0,1);
if isfield(settings,"speeds") && ~isempty(settings.speeds)
    speeds = double(settings.speeds(:));
end
end

function key = caseKey(caseInfo)
key = "default";
if isfield(caseInfo,"id") && ~isempty(caseInfo.id), key = string(caseInfo.id); end
end

function run = addStableFields(run,raw,requested)
loc = zeros(numel(requested),1);
rawSpeeds = double(raw.speed_mps(:));
for i = 1:numel(requested)
    scale = max([1;abs(requested(i));abs(rawSpeeds(isfinite(rawSpeeds)))]);
    match = find(abs(rawSpeeds-requested(i)) <= 32*eps(scale),1);
    if ~isempty(match), loc(i) = match; end
    if loc(i) == 0, loc(i) = min(i,numel(rawSpeeds)); end
end
loc(loc == 0) = 1;
run.speedIndex = raw.speedIndex(loc);
run.origin = raw.origin(loc);
run.passIndex = raw.passIndex(loc);
run.refinementReason = raw.refinementReason(loc);
run.solver_status = raw.status(loc);
run.solver_reason = raw.reason(loc);
end

function status = aggregateRunStatus(raw,execution)
if isfield(execution,"status") && string(execution.status) == "cancelled"
    status = "cancelled";
elseif isempty(raw.perSpeed) || ~any(raw.perSpeed.valid)
    status = "failed";
elseif any(~raw.perSpeed.valid)
    status = "partial";
else
    status = "completed";
end
end
