function [y,sd] = predictDOEModel(analysis,response,predictorTable)
%PREDICTDOEMODEL Predict a DOE response through its preferred surrogate.

response = char(string(response));
if ~isfield(analysis,'preferredModel') || ~isfield(analysis.preferredModel,response)
    error('predictDOEModel:unknownResponse','No preferred model exists for %s.',response);
end
if ~istable(predictorTable)
    error('predictDOEModel:badPredictors','Predictors must be supplied as a table.');
end
predictors = analysis.settings.predictors;
missing = setdiff(predictors,predictorTable.Properties.VariableNames);
if ~isempty(missing)
    error('predictDOEModel:missingPredictor', ...
        'Predictor table is missing: %s.',strjoin(missing,', '));
end

preferred = string(analysis.preferredModel.(response));
if preferred == "gp"
    if ~isfield(analysis,'gpModels') || ~isfield(analysis.gpModels,response) || ...
            isempty(analysis.gpModels.(response))
        error('predictDOEModel:missingGP','No GP model exists for %s.',response);
    end
    [y,sd] = predict(analysis.gpModels.(response),predictorTable{:,predictors});
else
    if ~isfield(analysis,'models') || ~isfield(analysis.models,response)
        error('predictDOEModel:missingQuadratic', ...
            'No quadratic model exists for %s.',response);
    end
    y = predict(analysis.models.(response),predictorTable(:,predictors));
    sd = nan(size(y));
end
y = y(:);
sd = sd(:);
end
