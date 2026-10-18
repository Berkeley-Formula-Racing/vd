function setup = fromRampSetup(fixed,source,mapInput,assets)
%FROMRAMPSETUP Translate frozen ramp inputs to the shared builder schema.
if nargin < 3 || isempty(mapInput), mapInput = ""; end
if nargin < 4 || isempty(assets), assets = struct(); end
validateattributes(fixed,{'struct'},{'scalar'},mfilename,'fixed');
validateattributes(source,{'struct'},{'scalar'},mfilename,'setup');
if ~isempty(assets)
    validateattributes(assets,{'struct'},{'scalar'},mfilename,'assets');
end
if ~isfield(fixed,'vehicle') || ~isstruct(fixed.vehicle) || ...
        ~isscalar(fixed.vehicle)
    error('vehicleConfig:badSchema', ...
        'fixed.vehicle must be one scalar source struct.');
end

setupNames = {'driverWeight_kg','rearWeightDistribution_percent', ...
    'frontSpringRate_lb_in','rearSpringRate_lb_in', ...
    'frontRideHeight_in','rearRideHeight_in'};
for k = 1:numel(setupNames)
    name = setupNames{k};
    if ~isfield(source,name)
        error('vehicleConfig:badSchema','Ramp setup is missing %s.',name);
    end
end

names = rampParameterNames();
vehicleNames = setdiff(names,setupNames,'stable');
vehicle = fixed.vehicle;
parameters = struct();
for k = 1:numel(vehicleNames)
    name = vehicleNames{k};
    if ~isfield(vehicle,name)
        error('vehicleConfig:badSchema', ...
            'fixed.vehicle is missing ramp source field %s.',name);
    end
    parameters.(name) = vehicle.(name);
end
for k = 1:numel(setupNames)
    name = setupNames{k};
    parameters.(name) = scalarFinite(source.(name),name);
end

fixedNames = {'gears','primary_reduction','torque_fn','Fx_parameters', ...
    'Fy_parameters','friction_scaling_factor'};
fixedConfig = struct();
for k = 1:numel(fixedNames)
    name = fixedNames{k};
    if ~isfield(vehicle,name)
        error('vehicleConfig:badSchema', ...
            'fixed.vehicle is missing ramp component source field %s.',name);
    end
    fixedConfig.(name) = vehicle.(name);
end
fixedConfig.model_type = "legacy";

frontArb = 0;
if isfield(fixed,'frontArbStiffness_NmPerRad')
    frontArb = fixed.frontArbStiffness_NmPerRad;
end

mapPath = "";
mapEnabled = false;
claScale = 1;
cdaScale = 1;
copOffset = 0;
selectedMap = [];
if isa(mapInput,'AeroMap')
    selectedMap = mapInput;
    mapPath = string(mapInput.sourcePath);
    mapEnabled = true;
elseif ~isempty(mapInput) && strlength(string(mapInput)) > 0
    mapPath = string(mapInput);
    mapEnabled = true;
    if isstruct(assets) && isfield(assets,'map') && ...
            isa(assets.map,'AeroMap') && ...
            string(assets.map.sourcePath) == mapPath
        selectedMap = assets.map;
    end
end
if ~isempty(selectedMap)
    claScale = selectedMap.claScale;
    cdaScale = selectedMap.cdaScale;
    copOffset = selectedMap.copOffset;
end

setup = struct( ...
    'schemaVersion',1, ...
    'profile',"ramp", ...
    'id',string(textField(source,'id',"legacy-ramp")), ...
    'label',string(textField(source,'label',"Legacy ramp setup")), ...
    'parameters',parameters, ...
    'fixed',fixedConfig, ...
    'aeroMap',struct('enabled',mapEnabled,'path',mapPath, ...
        'claScale',claScale,'cdaScale',cdaScale,'copOffset',copOffset), ...
    'suspension',struct('rollSplitPolicy',"derived", ...
        'frontArbStiffness_NmPerRad',frontArb, ...
        'rearArbStiffness_NmPerRad', ...
            scalarFinite(source.rearArbStiffness_NmPerRad, ...
                'rearArbStiffness_NmPerRad')));
end

function value = scalarFinite(value,name)
value = double(value);
if ~isscalar(value) || ~isfinite(value)
    error('vehicleConfig:invalidParameter', ...
        'Ramp source %s must be a finite scalar.',name);
end
end

function value = textField(source,name,default)
value = default;
if isfield(source,name) && ~isempty(source.(name))
    candidate = source.(name);
    if isstring(candidate) && isscalar(candidate) && ...
            (ismissing(candidate) || strlength(candidate) == 0)
        return
    end
    value = candidate;
end
end
