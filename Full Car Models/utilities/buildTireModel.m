function tire = buildTireModel(tireParams)
%BUILDTIREMODEL Construct the configured vehicle-facing tire model.

if ~isstruct(tireParams) || ~isscalar(tireParams) || ...
        ~isfield(tireParams,'model_type')
    error('buildTireModel:badConfig', ...
        'tireParams.model_type is required.');
end
modelType = lower(string(tireParams.model_type));
switch modelType
    case "legacy"
        tire = Tire2(tireParams.p_i,tireParams.Fx_parameters, ...
            tireParams.Fy_parameters,tireParams.friction_scaling_factor);
    case {"lc0_nd","nondimensional"}
        if ~isfield(tireParams,'model_artifact') || ...
                ~isfile(tireParams.model_artifact)
            error('buildTireModel:missingArtifact', ...
                'Nondimensional tire artifact is missing: %s', ...
                string(tireParams.model_artifact));
        end
        loaded = load(tireParams.model_artifact,'model');
        if ~isfield(loaded,'model')
            error('buildTireModel:badArtifact', ...
                'Nondimensional artifact must contain a variable named model.');
        end
        mode = getField(tireParams,'model_uncertainty',"nominal");
        tire = NondimensionalTire(loaded.model,tireParams.p_i, ...
            tireParams.friction_scaling_factor,mode);
    otherwise
        error('buildTireModel:badType', ...
            'Unsupported tire model type %s.',modelType);
end
end
function value = getField(s,name,defaultValue)
if isfield(s,name) && ~isempty(s.(name))
    value = s.(name);
else
    value = defaultValue;
end
end
