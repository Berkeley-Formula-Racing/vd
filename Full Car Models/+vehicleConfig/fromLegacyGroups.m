function setup = fromLegacyGroups(cP,aP,eP,DTp,Bp,tP)
%FROMLEGACYGROUPS Translate the six legacy parameter groups to source data.
names = vehicleConfig.parameterNames();
groups = {cP,aP,eP,DTp,Bp,tP};
parameters = struct();
for k = 1:numel(names)
    name = names{k};
    found = false;
    for groupIndex = 1:numel(groups)
        group = groups{groupIndex};
        if isstruct(group) && isscalar(group) && isfield(group,name)
            parameters.(name) = group.(name);
            found = true;
            break
        end
    end
    if ~found
        error('vehicleConfig:missingLegacyParameter', ...
            'Legacy parameter groups are missing sampled field %s.',name);
    end
end

fixed = struct();
fixed.gears = eP.gears;
fixed.primary_reduction = eP.primary_reduction;
fixed.torque_fn = eP.torque_fn;
fixed.Fx_parameters = tP.Fx_parameters;
fixed.Fy_parameters = tP.Fy_parameters;
fixed.friction_scaling_factor = tP.friction_scaling_factor;
fixed.model_type = fieldOr(tP,'model_type',"legacy");
fixed.model_artifact = fieldOr(tP,'model_artifact',"");
fixed.model_uncertainty = fieldOr(tP,'model_uncertainty',"nominal");

mapPath = fieldOr(aP,'map_path','');
mapPath = absolutePath(mapPath);
aeroMap = struct( ...
    'enabled',aeroMapEnabled(aP), ...
    'path',mapPath, ...
    'claScale',fieldOr(aP,'claScale',1), ...
    'cdaScale',fieldOr(aP,'cdaScale',1), ...
    'copOffset',fieldOr(aP,'copOffset',0));
suspension = struct( ...
    'rollSplitPolicy',"explicit", ...
    'frontArbStiffness_NmPerRad',0, ...
    'rearArbStiffness_NmPerRad',0);

setup = struct( ...
    'schemaVersion',1, ...
    'profile',"qss", ...
    'id',"legacy-qss", ...
    'label',"Legacy QSS", ...
    'parameters',parameters, ...
    'fixed',fixed, ...
    'aeroMap',aeroMap, ...
    'suspension',suspension);
end

function value = fieldOr(source,name,defaultValue)
if isfield(source,name)
    value = source.(name);
else
    value = defaultValue;
end
end

function path = absolutePath(path)
path = char(string(path));
if isempty(path)
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
end
