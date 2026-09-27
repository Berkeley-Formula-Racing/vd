function isValid = validateRequest(rawRequest)
%VALIDATEREQUEST Validate a raw Ramp-Speed request at the public boundary.
if ~isstruct(rawRequest) || ~isscalar(rawRequest) || ...
        ~isfield(rawRequest, 'settings') || ~isstruct(rawRequest.settings) || ...
        ~isscalar(rawRequest.settings)
    error('rampSpeed:invalidRequest', ...
        'A request must be a scalar struct with scalar struct settings.');
end

if ~isfield(rawRequest, 'rampType')
    error('rampSpeed:invalidRampType', ...
        'Request rampType must be lateral or longitudinal.');
end
rampType = normalizeText(rawRequest.rampType, 'rampSpeed:invalidRampType', ...
    'Request rampType must be lateral or longitudinal.');
if ~any(rampType == ["lateral", "longitudinal"])
    error('rampSpeed:invalidRampType', ...
        'Request rampType must be lateral or longitudinal.');
end

if ~isfield(rawRequest.settings, 'speeds_mps')
    error('rampSpeed:invalidSpeedDomain', ...
        'settings.speeds_mps must be a strictly increasing vector of finite positive speeds.');
end
speeds = rawRequest.settings.speeds_mps;
if ~isnumeric(speeds) || ~isreal(speeds) || ~isvector(speeds) || isempty(speeds)
    error('rampSpeed:invalidSpeedDomain', ...
        'settings.speeds_mps must be a strictly increasing vector of finite positive speeds.');
end
speeds = double(speeds(:));
if any(~isfinite(speeds)) || any(speeds <= 0) || any(diff(speeds) <= 0)
    error('rampSpeed:invalidSpeedDomain', ...
        'settings.speeds_mps must be a strictly increasing vector of finite positive speeds.');
end

speedPolicy = "fixed";
if isfield(rawRequest.settings, 'speedPolicy')
    speedPolicy = normalizeText(rawRequest.settings.speedPolicy, ...
        'rampSpeed:invalidSpeedPolicy', ...
        'settings.speedPolicy must be fixed or adaptive.');
end
if ~any(speedPolicy == ["fixed", "adaptive"])
    error('rampSpeed:invalidSpeedPolicy', ...
        'settings.speedPolicy must be fixed or adaptive.');
end
if rampType == "lateral" && speedPolicy == "adaptive"
    error('rampSpeed:unsupportedAdaptiveMode', ...
        'Adaptive speed sampling is not supported for lateral ramps.');
end

if isfield(rawRequest, 'execution')
    executionMode = requestExecutionMode(rawRequest.execution);
    if ~any(executionMode == ["serial", "parallel"])
        error('rampSpeed:invalidExecutionMode', ...
            'Execution mode must be serial or parallel.');
    end
end
if isfield(rawRequest, 'solverProfile')
    profile = rawRequest.solverProfile;
    if ~(ischar(profile) && isrow(profile)) && ...
            ~(isstring(profile) && isscalar(profile) && ~ismissing(profile))
        error('rampSpeed:invalidSolverProfile', ...
            'solverProfile must be a text scalar.');
    end
end

isValid = true;
end

function value = requestExecutionMode(execution)
if isstruct(execution) && isscalar(execution) && isfield(execution, 'mode')
    value = normalizeText(execution.mode, 'rampSpeed:invalidExecutionMode', ...
        'Execution mode must be serial or parallel.');
elseif ischar(execution) || (isstring(execution) && isscalar(execution))
    value = normalizeText(execution, 'rampSpeed:invalidExecutionMode', ...
        'Execution mode must be serial or parallel.');
else
    error('rampSpeed:invalidExecutionMode', ...
        'Execution must be a mode string or a scalar struct with a mode field.');
end
end

function value = normalizeText(value, errorId, message)
if ischar(value) && isrow(value)
    value = string(value);
elseif ~(isstring(value) && isscalar(value) && ~ismissing(value))
    error(errorId, '%s', message);
end
value = lower(strtrim(value));
if strlength(value) == 0
    error(errorId, '%s', message);
end
end
