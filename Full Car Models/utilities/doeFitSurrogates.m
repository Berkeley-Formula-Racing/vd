function models = doeFitSurrogates(U,metricTable,responseNames)
%DOEFITSURROGATES Fit adaptive DOE response, objective, and feasibility models.

validateattributes(U,{'numeric'},{'2d','real','finite','nonempty'})
if ~istable(metricTable) || height(metricTable) ~= size(U,1)
    error('doeFitSurrogates:badMetrics', ...
        'metricTable must have one row for every normalized design row.')
end
if ~ismember('valid',metricTable.Properties.VariableNames)
    error('doeFitSurrogates:missingValidity', ...
        'metricTable must contain a valid column.')
end
valid = metricTable.valid;
if ~islogical(valid) || ~isvector(valid)
    error('doeFitSurrogates:badValidity', ...
        'metricTable.valid must be a logical column.')
end

responseNames = string(responseNames(:));
if isempty(responseNames) || any(ismissing(responseNames) | strlength(responseNames) == 0) || ...
        numel(unique(responseNames)) ~= numel(responseNames)
    error('doeFitSurrogates:badResponses', ...
        'responseNames must contain unique, nonempty response names.')
end

models = struct();
models.responses = struct();
models.responseSpreads = struct();
models.responseNames = responseNames;
models.skippedResponses = strings(0,1);
for i = 1:numel(responseNames)
    responseName = responseNames(i);
    fieldName = matlab.lang.makeValidName(char(responseName));
    [model,spread,available] = fitResponse(U,metricTable,valid,responseName);
    if available
        models.responses.(fieldName) = model;
        models.responseSpreads.(fieldName) = spread;
    else
        models.skippedResponses(end+1,1) = responseName; %#ok<AGROW>
    end
end
if isempty(fieldnames(models.responses))
    error('doeFitSurrogates:noResponses', ...
        'Every requested adaptive response is unavailable or has zero spread.')
end

[models.objective,~,objectiveAvailable] = ...
    fitResponse(U,metricTable,valid,"objective_score");
if ~objectiveAvailable
    models.objective = [];
end

classes = unique(valid);
if numel(classes) == 2
    models.feasibility = fitcensemble(U,valid);
else
    models.feasibility = double(classes(1));
end
end

function [model,spread,available] = fitResponse(U,metricTable,valid,responseName)
model = [];
spread = NaN;
available = false;
fieldName = char(responseName);
if ~ismember(fieldName,metricTable.Properties.VariableNames)
    return
end
y = metricTable.(fieldName);
if ~isnumeric(y) || ~isreal(y) || ~isvector(y)
    return
end
y = y(:);
ok = valid & isfinite(y);
if nnz(ok) < 2
    return
end
spread = std(y(ok),0);
if ~isfinite(spread) || spread <= eps(max(abs(y(ok))))
    return
end
model = fitrgp(U(ok,:),y(ok), ...
    'KernelFunction','ardmatern32','Standardize',true, ...
    'FitMethod','exact','PredictMethod','exact');
available = true;
end
