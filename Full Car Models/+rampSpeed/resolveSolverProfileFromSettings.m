function [profile,settings] = resolveSolverProfileFromSettings(settings)
%RESOLVESOLVERPROFILEFROMSETTINGS Resolve profile fields on run settings.
%
% Ramp-speed entry points accept solverProfile in their existing settings
% struct. Legacy explicit solver option fields remain supported as overrides.

if nargin < 1 || isempty(settings)
    settings = struct();
end
if ~isstruct(settings) || ~isscalar(settings)
    error("rampSpeed:invalidSettings", ...
        "settings must be a scalar struct.");
end

selection = "accurate";
if isfield(settings,"solverProfile") && ~isempty(settings.solverProfile)
    selection = settings.solverProfile;
end

overrides = struct();
if isfield(settings,"solverOptions") && ~isempty(settings.solverOptions)
    if ~isstruct(settings.solverOptions) || ~isscalar(settings.solverOptions)
        error("rampSpeed:invalidSolverProfile", ...
            "settings.solverOptions must be a scalar struct.");
    end
    overrides = settings.solverOptions;
end
names = ["maxFunctionEvaluations","constraintTolerance", ...
    "stepTolerance","display"];
for name = names
    fieldName = char(name);
    if isfield(settings,fieldName) && ~isfield(overrides,fieldName)
        overrides.(fieldName) = settings.(fieldName);
    end
end

profile = rampSpeed.resolveSolverProfile(selection,overrides);
settings.solverProfile = profile.id;
settings.solverOptions = profile.solverOptions;
end
