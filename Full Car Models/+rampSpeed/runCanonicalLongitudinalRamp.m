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
requested = requestedSpeeds(settings);
grid = speedGrid(settings);
if isfield(grid,"range_mps") && ~isempty(grid.range_mps)
    requested = double(grid.range_mps(:));
end
if isempty(requested)
    requested = (5:2.5:30).';
end
% Fixed studies preserve the caller's requested order. Adaptive planning owns
% its own sorted seed/final grids, while the stable task IDs retain plan order.
requested = unique(double(requested(:)),'stable');
if any(isfinite(requested) & requested <= 0)
    error("rampSpeed:invalidSpeedDomain", ...
        "Finite longitudinal speeds must be positive.");
end
if isAdaptive(grid) && any(~isfinite(requested))
    requested = requested(isfinite(requested) & requested > 0);
    if isempty(requested)
        error("rampSpeed:invalidSpeedDomain", ...
            "Adaptive longitudinal runs require at least one finite positive speed.");
    end
end

request = struct("rampType","longitudinal","settings",settings, ...
    "solverProfile",profile.id,"setupKey",caseKey(caseInfo));
request.settings.speeds_mps = requested;
request.settings.speedPolicy = "fixed";
control = makeControl(callbacks);

if isAdaptive(grid)
    policy = adaptivePolicy(grid);
    speedPlan = rampSpeed.planAdaptiveSpeeds(requested,policy);
    seedReasons = repmat("",numel(speedPlan.seedSpeeds_mps),1);
    [tasks,nextId] = makeTasks(speedPlan.seedSpeeds_mps, ...
        ones(numel(speedPlan.seedSpeeds_mps),1), ...
        seedOrigins(speedPlan),seedReasons,1);
    [results,execution] = rampSpeed.executeSpeedPlan(request,tasks, ...
        @solveTask,control);
    allTasks = tasks;
    allResults = results;
    working = rampSpeed.assembleRun(allTasks,allResults,request,execution);
    while speedPlan.pass < policy.maxPasses && ...
            ~strcmp(execution.status,"cancelled")
        [speedPlan,report] = rampSpeed.refineAdaptiveSpeeds(speedPlan,working);
        if isempty(report.insertedSpeeds_mps)
            break
        end
        nextPass = speedPlan.pass;
        count = numel(report.insertedSpeeds_mps);
        [newTasks,nextId] = makeTasks(report.insertedSpeeds_mps, ...
            repmat(nextPass,count,1),repmat("refined",count,1), ...
            report.insertReasons,nextId);
        [newResults,newExecution] = rampSpeed.executeSpeedPlan( ...
            request,newTasks,@solveTask,control);
        allTasks = [allTasks;newTasks]; %#ok<AGROW>
        allResults = [allResults;newResults]; %#ok<AGROW>
        execution = combineExecution(execution,newExecution);
        working = rampSpeed.assembleRun(allTasks,allResults,request,execution);
    end
    if speedPlan.pass >= policy.maxPasses && strlength(string(speedPlan.stopReason)) == 0
        speedPlan.stopReason = "max_passes";
        speedPlan.status = "limited";
        if ~isempty(speedPlan.refinementHistory)
            speedPlan.refinementHistory(end).stopReason = "max_passes";
        end
    end
    finalSpeeds = speedPlan.speeds_mps(:);
    requestedForRun = requested;
    speedGridMeta = struct( ...
        "policy",speedPlan.policy, ...
        "requestedSpeeds_mps",speedPlan.requestedSpeeds_mps, ...
        "seedSpeeds_mps",speedPlan.seedSpeeds_mps, ...
        "finalSpeeds_mps",finalSpeeds, ...
        "passes",speedPlan.pass, ...
        "provenance",speedPlan.provenance, ...
        "refinementHistory",speedPlan.refinementHistory, ...
        "stopReason",string(speedPlan.stopReason), ...
        "exactRetrySpeeds_mps",requestedInvalid(working, ...
            speedPlan.requestedSpeeds_mps));
else
    [tasks,~] = makeTasks(requested,ones(numel(requested),1), ...
        repmat("requested",numel(requested),1), ...
        repmat("",numel(requested),1),1);
    [allResults,execution] = rampSpeed.executeSpeedPlan( ...
        request,tasks,@solveTask,control);
    allTasks = tasks;
    working = rampSpeed.assembleRun(allTasks,allResults,request,execution);
    finalSpeeds = requested;
    requestedForRun = requested;
    speedGridMeta = struct("mode","fixed", ...
        "requestedSpeeds_mps",requested, ...
        "seedSpeeds_mps",requested, ...
        "finalSpeeds_mps",requested, ...
        "passes",0, ...
        "provenance",table(requested,repmat("requested",numel(requested),1), ...
            zeros(numel(requested),1),repmat("user_requested",numel(requested),1), ...
            'VariableNames',{'speed_mps','source','pass','reason'}), ...
        "refinementHistory",struct([]),"stopReason","", ...
        "exactRetrySpeeds_mps",zeros(0,1));
end

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

function origins = seedOrigins(plan)
origins = repmat("seed",numel(plan.seedSpeeds_mps),1);
for i = 1:numel(origins)
    if i <= height(plan.provenance)
        origins(i) = string(plan.provenance.source(i));
    end
end
end

function execution = combineExecution(first,second)
execution = second;
execution.status = string(second.status);
execution.completedTasks = first.completedTasks + second.completedTasks;
execution.plannedTasks = first.plannedTasks + second.plannedTasks;
execution.evaluatedSpeedIndices = [first.evaluatedSpeedIndices(:); ...
    second.evaluatedSpeedIndices(:)];
execution.cancelled = first.cancelled || second.cancelled;
end

function control = makeControl(callbacks)
control = struct();
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

function values = requestedInvalid(raw,requested)
values = zeros(0,1);
if isempty(raw) || ~istable(raw.perSpeed), return, end
for value = requested(:).'
    index = find(abs(double(raw.perSpeed.speed_mps)-value) <= ...
        32*eps(max([1;abs(value);abs(double(raw.perSpeed.speed_mps(:)))])),1);
    if ~isempty(index) && ~raw.perSpeed.valid(index)
        values(end+1,1) = value; %#ok<AGROW>
    end
end
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
    error("rampSpeed:invalidAdaptivePolicy","speedGrid must be a scalar struct or mode.");
end
end

function tf = isAdaptive(grid)
mode = lower(strtrim(string(grid.mode)));
tf = any(mode == ["preview","accurate","highaccuracy"]);
end

function policy = adaptivePolicy(grid)
mode = "accurate";
if isfield(grid,"mode"), mode = grid.mode; end
overrides = struct();
names = ["baseSpacing_mps","minRefinementSpacing_mps", ...
    "balanceTolerance_fraction","relativeTolerance", ...
    "forceAbsoluteTolerance_N","accelerationAbsoluteTolerance_mps2", ...
    "residualRefinementFraction","maxPasses","maxPoints","stateFields"];
for name = names
    field = char(name);
    if isfield(grid,field), overrides.(field) = grid.(field); end
end
if isfield(grid,"policy") && isstruct(grid.policy)
    names = fieldnames(grid.policy);
    for i = 1:numel(names), overrides.(names{i}) = grid.policy.(names{i}); end
end
policy = rampSpeed.adaptiveSpeedPolicy(mode,overrides);
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
