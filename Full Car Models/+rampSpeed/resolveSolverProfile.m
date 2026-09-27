function profile = resolveSolverProfile(selection,solverOverrides)
%RESOLVESOLVERPROFILE Resolve a profile name and validated solver overrides.
%
% Passing an empty selection (or omitting it) selects the accurate profile.
% A resolved profile struct may also be passed when restoring saved metadata.

if nargin < 1 || isempty(selection)
    selection = "accurate";
end
if nargin < 2 || isempty(solverOverrides)
    solverOverrides = struct();
end

if isstruct(selection)
    if ~isfield(selection,'powertrainModel')
        selection.powertrainModel = "continuousEnvelope";
    end
    [valid,issues] = rampSpeed.validateSolverProfile(selection);
    if ~valid
        error("rampSpeed:invalidSolverProfile", ...
            "Invalid solver profile: %s",strjoin(issues,"; "));
    end
    profile = selection;
else
    [isText,id] = scalarText(selection);
    if ~isText || strlength(strtrim(id)) == 0
        error("rampSpeed:invalidSolverProfile", ...
            "Solver profile selection must be scalar text or a resolved profile struct.");
    end
    profiles = rampSpeed.solverProfiles();
    profileIds = lower(string({profiles.id}));
    match = find(profileIds == lower(strtrim(id)),1);
    if isempty(match)
        error("rampSpeed:unknownSolverProfile", ...
            "Unknown solver profile '%s'. Choose accurate, fastPreview, or approximateAeroPreview.", ...
            id);
    end
    profile = profiles(match);
end

if ~isstruct(solverOverrides) || ~isscalar(solverOverrides)
    error("rampSpeed:invalidSolverProfile", ...
        "Solver overrides must be a scalar struct.");
end
overrideNames = string(fieldnames(solverOverrides));
allowed = ["maxFunctionEvaluations","constraintTolerance", ...
    "stepTolerance","display"];
for name = overrideNames(:).'
    if ~any(name == allowed)
        error("rampSpeed:invalidSolverProfile", ...
            "Unsupported solver option '%s'.",name);
    end
    profile.solverOptions.(char(name)) = solverOverrides.(char(name));
end

profile.id = string(profile.id);
profile.label = string(profile.label);
profile.aeroMode = string(profile.aeroMode);
profile.powertrainModel = string(profile.powertrainModel);
profile.description = string(profile.description);
profile.solverOptions.display = string(profile.solverOptions.display);
[valid,issues] = rampSpeed.validateSolverProfile(profile);
if ~valid
    error("rampSpeed:invalidSolverProfile", ...
        "Invalid solver profile: %s",strjoin(issues,"; "));
end
end

function [valid,value] = scalarText(candidate)
valid = false;
value = "";
if isstring(candidate) && isscalar(candidate) && ~ismissing(candidate)
    value = string(candidate);
    valid = true;
elseif ischar(candidate) && (isrow(candidate) || isempty(candidate))
    value = string(candidate);
    valid = true;
end
end
