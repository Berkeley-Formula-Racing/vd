function [valid,issues] = validateSolverProfile(profile)
%VALIDATESOLVERPROFILE Check a resolved ramp-speed profile for safe use.

issues = strings(0,1);
if ~isstruct(profile) || ~isscalar(profile)
    issues(end+1,1) = "profile must be a scalar struct.";
    valid = false;
    return
end

required = ["id","label","version","solverOptions", ...
    "aeroMode","powertrainModel","approximate","description"];
for name = required
    if ~isfield(profile,char(name))
        issues(end+1,1) = "profile is missing required field '" + name + "'.";
    end
end
if ~isempty(issues)
    valid = false;
    return
end

[hasId,id] = scalarText(profile.id);
if ~hasId || ~any(id == ["accurate","fastPreview", ...
        "approximateAeroPreview"])
    issues(end+1,1) = "profile.id is not a supported ramp-speed profile.";
end
[hasLabel,label] = scalarText(profile.label);
if ~hasLabel || strlength(strtrim(label)) == 0
    issues(end+1,1) = "profile.label must be nonempty scalar text.";
end
[hasDescription,description] = scalarText(profile.description);
if ~hasDescription || strlength(strtrim(description)) == 0
    issues(end+1,1) = "profile.description must be nonempty scalar text.";
end

if ~isnumeric(profile.version) || ~isscalar(profile.version) || ...
        ~isfinite(profile.version) || profile.version < 1 || ...
        profile.version ~= fix(profile.version)
    issues(end+1,1) = "profile.version must be a positive integer.";
end

[hasAeroMode,aeroMode] = scalarText(profile.aeroMode);
if ~hasAeroMode || ~any(aeroMode == ["coupled","static"])
    issues(end+1,1) = "profile.aeroMode must be 'coupled' or 'static'.";
end

[hasPowertrainModel,powertrainModel] = scalarText(profile.powertrainModel);
if ~hasPowertrainModel || powertrainModel ~= "continuousEnvelope"
    issues(end+1) = "profile.powertrainModel must be 'continuousEnvelope'.";
end
if ~islogical(profile.approximate) || ~isscalar(profile.approximate)
    issues(end+1,1) = "profile.approximate must be a scalar logical.";
elseif hasId && any(id == ["accurate","fastPreview"]) && ...
        (profile.approximate || ~hasAeroMode || aeroMode ~= "coupled")
    issues(end+1,1) = "accurate and fastPreview require coupled, non-approximate aero.";
elseif hasId && id == "approximateAeroPreview" && ...
        (~profile.approximate || ~hasAeroMode || aeroMode ~= "static")
    issues(end+1,1) = "approximateAeroPreview requires static aero and approximate=true.";
end

if ~isstruct(profile.solverOptions) || ~isscalar(profile.solverOptions)
    issues(end+1,1) = "profile.solverOptions must be a scalar struct.";
else
    allowedOptions = ["maxFunctionEvaluations","constraintTolerance", ...
        "stepTolerance","display"];
    optionNames = string(fieldnames(profile.solverOptions));
    for name = optionNames(:).'
        if ~any(name == allowedOptions)
            issues(end+1,1) = "unsupported solver option '" + name + "'.";
            continue
        end
        value = profile.solverOptions.(char(name));
        switch name
            case "maxFunctionEvaluations"
                if ~isnumeric(value) || ~isscalar(value) || ...
                        ~isfinite(value) || value < 1 || value ~= fix(value)
                    issues(end+1,1) = ...
                        "maxFunctionEvaluations must be a positive integer.";
                end
            case {"constraintTolerance","stepTolerance"}
                if ~isnumeric(value) || ~isscalar(value) || ...
                        ~isfinite(value) || value <= 0
                    issues(end+1,1) = name + " must be a positive finite scalar.";
                end
            case "display"
                [hasDisplay,display] = scalarText(value);
                allowedDisplay = ["off","final","notify","notify-detailed", ...
                    "iter","iter-detailed"];
                if ~hasDisplay || ~any(lower(display) == allowedDisplay)
                    issues(end+1,1) = ...
                        "display must be a supported fmincon display mode.";
                end
        end
    end
end

valid = isempty(issues);
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
