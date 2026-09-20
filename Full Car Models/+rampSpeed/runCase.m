function run = runCase(car,caseInfo,request,callbacks)
%RUNCASE Dispatch one setup to the normalized ramp adapter.

if nargin < 3 || isempty(request)
    request = struct();
end
if nargin < 4 || isempty(callbacks)
    callbacks = struct();
end
if ~isstruct(request) || ~isscalar(request)
    error("rampSpeed:invalidRequest", ...
        "request must be a scalar struct.");
end
if ~isstruct(caseInfo) || ~isscalar(caseInfo)
    error("rampSpeed:invalidCaseInfo", ...
        "caseInfo must be a scalar struct.");
end
if ~isstruct(callbacks) || ~isscalar(callbacks)
    error("rampSpeed:invalidCallbacks", ...
        "callbacks must be a scalar struct.");
end

% Direct callers may inject a runner too. runStudy invokes the injected
% handle itself, so the default handle must not recurse here.
if isfield(request,'runCaseFcn') && ~isempty(request.runCaseFcn) && ...
        isa(request.runCaseFcn,'function_handle') && ...
        ~strcmp(func2str(request.runCaseFcn),'rampSpeed.runCase')
    run = request.runCaseFcn(car,caseInfo,request,callbacks);
    return
end

type = normalizeType(getField(request,'rampType',"lateral"));
settings = getField(request,'settings',struct());
if isempty(settings)
    settings = struct();
end
if ~isstruct(settings) || ~isscalar(settings)
    error("rampSpeed:invalidSettings", ...
        "request.settings must be a scalar struct.");
end

callbacks = addProgressQueue(callbacks,request,caseInfo);
switch type
    case "lateral"
        if ~isfield(settings,'mode') || isempty(settings.mode) || ...
                strlength(string(settings.mode)) == 0
            settings.mode = "coast";
        end
        run = rampSpeed.runLateralRamp(car,settings,caseInfo,callbacks);
    case "longitudinal"
        run = rampSpeed.runLongitudinalRamp(car,settings,caseInfo,callbacks);
    otherwise
        error("rampSpeed:unsupportedType", ...
            "Unsupported ramp type: %s",type);
end
end

function callbacks = addProgressQueue(callbacks,request,caseInfo)
if ~isfield(callbacks,'onProgress')
    callbacks.onProgress = [];
end
if ~isfield(request,'progressQueue') || isempty(request.progressQueue)
    return
end

prior = callbacks.onProgress;
queue = request.progressQueue;
callbacks.onProgress = @(event)publishProgress(prior,queue,event, ...
    caseInfo,request);
end

function publishProgress(prior,queue,rawEvent,caseInfo,request)
event = struct( ...
    "phase","case", ...
    "caseId",string(getField(caseInfo,'id',"")), ...
    "speedIndex",numericField(rawEvent,'speedIndex',NaN), ...
    "speed_mps",numericField(rawEvent,'speed_mps',NaN), ...
    "completedCases",numericField(rawEvent,'completedCases',0), ...
    "totalCases",numericField(request,'totalCases',0), ...
    "message",stringField(rawEvent,'message',""));
if ~isempty(prior)
    try
        prior(event);
    catch
        % Progress callbacks are advisory and must not abort a case.
    end
end
sendProgress(queue,event);
end

function sendProgress(queue,event)
try
    if isa(queue,'function_handle')
        queue(event);
    else
        send(queue,event);
    end
catch
    % Progress is advisory; a disconnected UI queue must not abort a solve.
end
end

function value = normalizeType(value)
value = lower(strtrim(string(value)));
if ~isscalar(value)
    error("rampSpeed:unsupportedType", ...
        "rampType must be scalar text.");
end
switch value
    case {"lateral","lateral-limit","lateral limit"}
        value = "lateral";
    case {"longitudinal","pure-longitudinal", ...
            "pure longitudinal","pure_longitudinal"}
        value = "longitudinal";
    otherwise
        error("rampSpeed:unsupportedType", ...
            "Unsupported ramp type: %s",value);
end
end

function value = getField(S,name,default)
value = default;
if isstruct(S) && isfield(S,name) && ~isempty(S.(name))
    value = S.(name);
end
end

function value = numericField(S,name,default)
value = default;
if isstruct(S) && isfield(S,name) && ~isempty(S.(name))
    candidate = double(S.(name));
    if isscalar(candidate)
        value = candidate;
    end
end
end

function value = stringField(S,name,default)
value = string(default);
if isstruct(S) && isfield(S,name) && ~isempty(S.(name))
    candidate = string(S.(name));
    if isscalar(candidate)
        value = candidate;
    end
end
end
