function [car,setup,derived] = buildCarFromSetup(config,setup)
%BUILDCARFROMSETUP Build exactly one solver car from a setup specification.

if nargin < 2 || isempty(setup)
    setup = config.defaultSetup;
end
validateattributes(config,{'struct'},{'scalar'},mfilename,'config');
validateattributes(setup,{'struct'},{'scalar'},mfilename,'setup');
[setup,~] = rampSpeed.normalizeSetupSpec(config,setup);
required = {'baselineVersion','id','label','aeroMapId', ...
    'rearArbStiffness_NmPerRad','frontSpringRate_lb_in', ...
    'rearSpringRate_lb_in','frontRideHeight_in','rearRideHeight_in', ...
    'driverWeight_kg','rearWeightDistribution_percent'};
for i = 1:numel(required)
    if ~isfield(setup,required{i})
        error('rampSpeed:invalidSetup', ...
            'Setup is missing required field %s.',required{i});
    end
end
setup = validateDiscrete(setup,config.options,'rearArbStiffness_NmPerRad');
setup = validateDiscrete(setup,config.options,'frontSpringRate_lb_in');
setup = validateDiscrete(setup,config.options,'rearSpringRate_lb_in');
setup.driverWeight_kg = finiteScalar(setup.driverWeight_kg, ...
    'driverWeight_kg');
setup.rearWeightDistribution_percent = finiteScalar( ...
    setup.rearWeightDistribution_percent,'rearWeightDistribution_percent');
if setup.driverWeight_kg < 0
    error('rampSpeed:invalidSetup', ...
        'driverWeight_kg must be nonnegative.');
end
if setup.rearWeightDistribution_percent <= 0 || ...
        setup.rearWeightDistribution_percent >= 100
    error('rampSpeed:invalidSetup', ...
        'rearWeightDistribution_percent must be between 0 and 100.');
end
setup.frontRideHeight_in = finiteScalar(setup.frontRideHeight_in, ...
    'frontRideHeight_in');
setup.rearRideHeight_in = finiteScalar(setup.rearRideHeight_in, ...
    'rearRideHeight_in');
setup.id = string(setup.id);
setup.label = string(setup.label);
if strlength(strtrim(setup.id)) == 0 || strlength(strtrim(setup.label)) == 0
    error('rampSpeed:invalidSetup','Setup ID and label must be non-empty.');
end

catalog = rampSpeed.aeroMapCatalog();
mapIds = string({catalog.id});
mapIndex = find(mapIds == string(setup.aeroMapId),1);
assets = rampSpeed.loadRampAssets(config,setup);
mapPath = assets.mapPath;
map = assets.map;

vehicle = config.fixedParameters.vehicle;
frontOffset = setup.frontRideHeight_in - vehicle.map_reference_front_ride_height_in;
rearOffset = setup.rearRideHeight_in - vehicle.map_reference_rear_ride_height_in;
outside = frontOffset < map.frontRangeIn(1) || ...
    frontOffset > map.frontRangeIn(2) || rearOffset < map.rearRangeIn(1) || ...
    rearOffset > map.rearRangeIn(2);
if outside
    warning('rampSpeed:rideHeightOutsideAeroMap', ...
        ['Setup %s uses ride heights outside the %s aero-map envelope; ' ...
        'bounded nearest extrapolation will be used.'],setup.id,setup.aeroMapId);
end

[car,derived] = buildSingleRampCar(config.fixedParameters,setup,map,assets);
derived.aeroMapId = string(setup.aeroMapId);
derived.aeroMapLabel = string(catalog(mapIndex).label);
derived.aeroMapRelativePath = string(catalog(mapIndex).relativePath);
derived.aeroMapFingerprint = assets.fingerprints.aeroMap;
derived.assetFingerprints = assets.fingerprints;
derived.aeroMapOutsideEnvelope = outside;
setup.source = getStringField(setup,'source',"user");
setup.isBaseline = isBaselineSetup(setup,config);
end

function setup = validateDiscrete(setup,options,name)
value = finiteScalar(setup.(name),name);
allowed = double(options.(name));
if ~any(value == allowed)
    error('rampSpeed:invalidSetupOption', ...
        '%s=%g is not one of the baseline-defined choices [%s].', ...
        name,value,strjoin(string(allowed),', '));
end
setup.(name) = value;
end

function value = finiteScalar(value,name)
value = double(value);
if ~isscalar(value) || ~isfinite(value)
    error('rampSpeed:invalidSetup', ...
        '%s must be a finite scalar.',name);
end
end

function value = getStringField(S,name,default)
value = string(default);
if isfield(S,name) && ~isempty(S.(name))
    value = string(S.(name));
end
end

function tf = isBaselineSetup(setup,config)
tf = string(setup.id) == string(config.id);
end
