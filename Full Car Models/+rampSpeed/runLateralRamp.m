function run = runLateralRamp(car,settings,caseInfo,callbacks)
%RUNLATERALRAMP Run rampSweep and return a normalized lateral run.

if nargin < 2 || isempty(settings)
    settings = struct();
end
if nargin < 3 || isempty(caseInfo)
    caseInfo = struct();
end
if nargin < 4 || isempty(callbacks)
    callbacks = struct();
end
if ~isstruct(settings) || ~isscalar(settings)
    error("rampSpeed:invalidSettings", ...
        "settings must be a scalar struct.");
end
if ~isstruct(caseInfo) || ~isscalar(caseInfo)
    error("rampSpeed:invalidCaseInfo", ...
        "caseInfo must be a scalar struct.");
end
if ~isstruct(callbacks) || ~isscalar(callbacks)
    error("rampSpeed:invalidCallbacks", ...
        "callbacks must be a scalar struct.");
end

progressFcn = [];
cancelFcn = [];
if isfield(callbacks,"onProgress") && ~isempty(callbacks.onProgress)
    progressFcn = callbacks.onProgress;
end
if isfield(callbacks,"isCancelled") && ~isempty(callbacks.isCancelled)
    cancelFcn = callbacks.isCancelled;
end

rawOptions = settings;
rawOptions.progressFcn = progressFcn;
rawOptions.cancelFcn = cancelFcn;

started = datetime('now');
runMeta = struct("source","rampSpeed.runLateralRamp", ...
    "started",started,"completed",datetime.empty, ...
    "warnings",strings(0,1),"errors",strings(0,1));
failure = [];

try
    raw = rampSweep(car,rawOptions);
catch ME
    if ~isHandledRampError(ME)
        rethrow(ME)
    end
    failure = ME;
    if strcmp(ME.identifier,"rampSweep:cancelled")
        runMeta.status = "cancelled";
    else
        runMeta.status = "failed";
    end
    runMeta.errors = string(ME.message);
    raw = failurePayload(settings,runMeta.status,ME);
end

if isfield(raw,"status") && ~isempty(raw.status)
    runMeta.status = string(raw.status);
end
runMeta.completed = datetime('now');

requested = requestedSpeeds(settings,raw);
normalizationSettings = settings;
normalizationSettings.speeds = normalizationSpeeds(raw,requested).';
normalizationSettings.mode = effectiveMode(settings,raw);
diagnostics = fillSpeedErrorSpeeds(speedErrorsForRun(raw,failure),requested);
if isstruct(raw)
    raw.speedErrors = diagnostics;
end
runMeta.speedErrors = diagnostics;
run = rampSpeed.normalizeRampResult(raw,"lateral",normalizationSettings, ...
    caseInfo,runMeta);

% The schema's settings describe the request, including the actual speed
% vector used when rampSweep defaults were requested.
run.settings = settings;
run.settings.speeds = requested.';
run.settings.mode = run.mode;
run.runMeta.requestedSpeeds_mps = requested;
run.perSpeed = rebuildPerSpeed(run,raw,requested,diagnostics, ...
    missingReason(failure,run));
if ~isempty(failure)
    run.status = runMeta.status;
end
end

function value = requestedSpeeds(settings,raw)
if isfield(settings,"speeds") && ~isempty(settings.speeds)
    value = double(settings.speeds(:));
    return
end

if isstruct(raw) && isfield(raw,"settings") && ...
        isstruct(raw.settings) && isfield(raw.settings,"speeds") && ...
        ~isempty(raw.settings.speeds)
    value = double(raw.settings.speeds(:));
    return
end

value = rawSpeeds(raw);
if isempty(value)
    value = (5:2.5:30).';
end
end

function value = normalizationSpeeds(raw,requested)
value = rawSpeeds(raw);
if isempty(value)
    value = requested;
end
end

function value = rawSpeeds(raw)
value = zeros(0,1);
if ~isstruct(raw) || ~isscalar(raw) || ~isfield(raw,"perSpeed")
    return
end
data = raw.perSpeed;
if istable(data)
    names = data.Properties.VariableNames;
    for name = ["speed_mps","vCar","speed","long_vel"]
        if any(strcmp(names,name))
            value = double(data.(name)(:));
            return
        end
    end
elseif isstruct(data) && ~isempty(data)
    for name = ["speed_mps","vCar","speed","long_vel"]
        if isfield(data,name)
            value = double([data.(name)].');
            return
        end
    end
end
end

function indices = rawSpeedIndices(raw,n)
indices = (1:n).';
if ~isstruct(raw) || ~isscalar(raw) || ~isfield(raw,"perSpeed")
    indices = zeros(n,1);
    return
end
data = raw.perSpeed;
if isempty(data) || (istable(data) && height(data) == 0)
    indices = zeros(n,1);
    return
end
if istable(data) && any(strcmp(data.Properties.VariableNames,"speed_index"))
    candidate = double(data.speed_index(:));
elseif isstruct(data) && isfield(data,"speed_index")
    candidate = double([data.speed_index].');
else
    return
end
if numel(candidate) == n
    indices = candidate;
end
end

function T = rebuildPerSpeed(run,raw,requested,diagnostics,reason)
templateSettings = run.settings;
templateSettings.speeds = requested.';
template = rampSpeed.makeRun("lateral",run.mode,templateSettings, ...
    struct("id",run.caseId));
T = template.perSpeed;
T.speed_mps = requested;

source = run.perSpeed;
nTarget = numel(requested);
nSource = height(source);
indices = rawSpeedIndices(raw,nSource);
used = false(nTarget,1);
for i = 1:nSource
    target = indices(i);
    if isfinite(target) && target >= 1 && target <= nTarget && ...
            target == fix(target) && ~used(target)
        T(target,:) = source(i,:);
        used(target) = true;
    end
end

missing = ~used;
T.valid(missing) = false;
T.status(missing) = "failed";
missingReasons = repmat(string(reason),nTarget,1);
for i = 1:numel(diagnostics)
    target = diagnostics(i).speed_index;
    if isfinite(target) && target >= 1 && target <= nTarget && ...
            target == fix(target)
        missingReasons(target) = formatSpeedError(diagnostics(i));
    end
end
T.reason(missing) = missingReasons(missing);
T.speed_mps = requested;
end

function reason = missingReason(failure,run)
if ~isempty(failure)
    reason = string(failure.message);
elseif run.status == "cancelled"
    reason = "cancelled before this speed was solved";
else
    reason = "no solved ramp returned for requested speed";
end
end

function raw = failurePayload(settings,status,ME)
raw = struct();
raw.perSpeed = table();
raw.points = table();
raw.settings = settings;
raw.status = status;
raw.speedErrors = exceptionSpeedErrors(ME);
raw.error = struct("identifier",string(ME.identifier), ...
    "message",string(ME.message));
end

function tf = isHandledRampError(ME)
tf = any(strcmp(ME.identifier,["rampSweep:cancelled", ...
    "rampSweep:noSolution","rampSpeed:noSolution"]));
end

function mode = effectiveMode(settings,raw)
mode = "";
if isfield(settings,"mode") && ~isempty(settings.mode)
    mode = string(settings.mode);
elseif isstruct(raw) && isfield(raw,"settings") && ...
        isstruct(raw.settings) && isfield(raw.settings,"mode") && ...
        ~isempty(raw.settings.mode)
    mode = string(raw.settings.mode);
end
if isempty(mode) || strlength(mode(1)) == 0
    mode = "balanced";
else
    mode = lower(mode(1));
end
end

function errors = speedErrorsForRun(raw,failure)
errors = emptySpeedErrors();
if isstruct(raw) && isfield(raw,"speedErrors") && ...
        ~isempty(raw.speedErrors)
    errors = raw.speedErrors;
elseif ~isempty(failure)
    errors = exceptionSpeedErrors(failure);
end
end

function errors = fillSpeedErrorSpeeds(errors,requested)
for i = 1:numel(errors)
    target = errors(i).speed_index;
    if isfinite(target) && target >= 1 && target <= numel(requested) && ...
            target == fix(target) && ~isfinite(errors(i).speed_mps)
        errors(i).speed_mps = requested(target);
    end
end
end

function errors = exceptionSpeedErrors(ME)
errors = emptySpeedErrors();
if isempty(ME) || isempty(ME.cause)
    return
end
for i = 1:numel(ME.cause)
    cause = ME.cause{i};
    identifier = string(cause.identifier);
    if startsWith(identifier,"rampSweep:speed:")
        speedIndex = str2double(extractAfter(identifier,"rampSweep:speed:"));
        original = cause;
        if ~isempty(cause.cause)
            original = cause.cause{1};
        end
        errors(end+1) = struct( ...
            "speed_index",speedIndex, ...
            "speed_mps",NaN, ...
            "identifier",string(original.identifier), ...
            "message",string(original.message), ...
            "stack",original.stack); %#ok<AGROW>
    else
        nested = exceptionSpeedErrors(cause);
        if ~isempty(nested)
            errors(end+1:end+numel(nested)) = nested;
        end
    end
end
end

function reason = formatSpeedError(entry)
identifier = string(entry.identifier);
message = string(entry.message);
if strlength(identifier) > 0
    reason = identifier + ": " + message;
else
    reason = message;
end
end

function errors = emptySpeedErrors()
errors = struct('speed_index',{},'speed_mps',{},'identifier',{}, ...
    'message',{},'stack',{});
end
