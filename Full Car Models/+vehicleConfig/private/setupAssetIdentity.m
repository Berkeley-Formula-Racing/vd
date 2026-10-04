function identity = setupAssetIdentity(setup)
%SETUPASSETIDENTITY Resolve the source paths that define a QSS asset context.
required = {'profile','fixed','aeroMap'};
if ~isstruct(setup) || ~isscalar(setup) || ~all(isfield(setup,required)) || ...
        ~isstruct(setup.fixed) || ~isscalar(setup.fixed) || ...
        ~isstruct(setup.aeroMap) || ~isscalar(setup.aeroMap)
    error('vehicleConfig:badSchema', ...
        'Asset loading requires one setup with profile, fixed, and aeroMap fields.');
end
fixedRequired = {'model_type','model_artifact'};
mapRequired = {'enabled','path','claScale','cdaScale','copOffset'};
if ~all(isfield(setup.fixed,fixedRequired)) || ...
        ~all(isfield(setup.aeroMap,mapRequired))
    error('vehicleConfig:badSchema', ...
        'Asset loading requires the reviewed fixed and aeroMap fields.');
end
if ~(ischar(setup.profile) || (isstring(setup.profile) && isscalar(setup.profile))) || ...
        ~(ischar(setup.fixed.model_type) || ...
            (isstring(setup.fixed.model_type) && isscalar(setup.fixed.model_type)))
    error('vehicleConfig:badSchema','Profile and tire model type must be scalar text.');
end
if ~(isnumeric(setup.aeroMap.enabled) || islogical(setup.aeroMap.enabled)) || ...
        ~isscalar(setup.aeroMap.enabled) || ~isreal(setup.aeroMap.enabled) || ...
        ~isfinite(double(setup.aeroMap.enabled))
    error('vehicleConfig:badSchema','aeroMap.enabled must be one finite scalar flag.');
end

helperDir = fileparts(mfilename('fullpath'));
packageDir = fileparts(helperDir);
modelRoot = fileparts(packageDir);
componentRoot = fullfile(modelRoot,'carComponents');
modelType = lower(string(setup.fixed.model_type));
modelPath = "";
if any(modelType == ["lc0_nd","nondimensional"])
    modelPath = resolvePath(setup.fixed.model_artifact);
end

identity = struct( ...
    'profile',lower(string(setup.profile)), ...
    'aeroMapEnabled',logical(setup.aeroMap.enabled), ...
    'aeroMapPath',resolvePath(setup.aeroMap.path), ...
    'claScale',setup.aeroMap.claScale, ...
    'cdaScale',setup.aeroMap.cdaScale, ...
    'copOffset',setup.aeroMap.copOffset, ...
    'tireModelType',modelType, ...
    'tireModelPath',modelPath, ...
    'camberModelPath',string(fullfile(componentRoot,'camber_models_fast.mat')), ...
    'camberRatioPath',string(fullfile(componentRoot,'camberratiossmoothed.mat')));
end

function path = resolvePath(value)
if isempty(value)
    path = "";
    return
end
if ~(ischar(value) || (isstring(value) && isscalar(value)))
    error('vehicleConfig:badSchema','Asset paths must be scalar text.');
end
path = char(string(value));
if isempty(path)
    path = "";
    return
end
if ispc
    isAbsolute = ~isempty(regexp(path,'^[A-Za-z]:[\\/]|^\\\\','once'));
else
    isAbsolute = startsWith(path,'/');
end
if ~isAbsolute
    path = fullfile(pwd,path);
end
[found,attributes] = fileattrib(path);
if found
    path = attributes.Name;
end
path = string(path);
end
