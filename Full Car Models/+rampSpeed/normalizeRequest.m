function request = normalizeRequest(rawRequest)
%NORMALIZEREQUEST Return the canonical Ramp-Speed request representation.
rampSpeed.validateRequest(rawRequest);
request = rawRequest;

request.rampType = normalizeText(rawRequest.rampType, ...
    'rampSpeed:invalidRampType', 'rampType must be lateral or longitudinal.');
request.settings.speeds_mps = double(rawRequest.settings.speeds_mps(:));

if isfield(rawRequest.settings, 'speedPolicy')
    request.settings.speedPolicy = normalizeText(rawRequest.settings.speedPolicy, ...
        'rampSpeed:invalidSpeedPolicy', 'settings.speedPolicy must be fixed or adaptive.');
else
    request.settings.speedPolicy = "fixed";
end

if isfield(rawRequest, 'execution')
    if isstruct(rawRequest.execution)
        request.execution.mode = normalizeText(rawRequest.execution.mode, ...
            'rampSpeed:invalidExecutionMode', 'Execution mode must be serial or parallel.');
    else
        request.execution = struct('mode', normalizeText(rawRequest.execution, ...
            'rampSpeed:invalidExecutionMode', 'Execution mode must be serial or parallel.'));
    end
else
    request.execution = struct('mode', "serial");
end

if isfield(rawRequest, 'caseIds')
    request.caseIds = normalizeCaseIds(rawRequest.caseIds);
else
    request.caseIds = strings(0, 1);
end

if isfield(rawRequest, 'solverProfile')
    request.solverProfile = normalizeText(rawRequest.solverProfile, ...
        'rampSpeed:invalidSolverProfile', 'solverProfile must be a text scalar.');
else
    request.solverProfile = "";
end
end

function ids = normalizeCaseIds(ids)
if ischar(ids) || isstring(ids) || iscellstr(ids)
    ids = string(ids(:));
    if any(ismissing(ids)) || any(strlength(strtrim(ids)) == 0)
        error('rampSpeed:invalidCaseIds', 'caseIds must not contain missing or empty IDs.');
    end
elseif isnumeric(ids) && isreal(ids) && (isvector(ids) || isempty(ids))
    ids = ids(:);
else
    error('rampSpeed:invalidCaseIds', ...
        'caseIds must be a text vector or a real numeric vector.');
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
