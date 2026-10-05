function [car,derived] = buildSingleRampCar(fixed,setup,mapInput,assets)
%BUILDSINGLERAMPCAR Construct one ramp-speed solver car from frozen inputs.
%   FIXED contains the ramp-speed baseline parameters. SETUP contains only
%   the editable setup values. This function deliberately creates one Car.

if nargin < 3 || isempty(mapInput)
    mapInput = "";
end
if nargin < 4 || isempty(assets)
    assets = struct();
end
validateattributes(fixed,{'struct'},{'scalar'},mfilename,'fixed');
validateattributes(setup,{'struct'},{'scalar'},mfilename,'setup');

required = ["rearArbStiffness_NmPerRad","frontSpringRate_lb_in", ...
    "rearSpringRate_lb_in","frontRideHeight_in","rearRideHeight_in", ...
    "driverWeight_kg","rearWeightDistribution_percent"];
for name = required
    if ~isfield(setup,char(name))
        error('buildSingleRampCar:missingSetupField', ...
            'Setup is missing %s.',name);
    end
end

driverWeight = scalarFinite(setup.driverWeight_kg,'driverWeight_kg');
rearWeightDistribution = scalarFinite( ...
    setup.rearWeightDistribution_percent,'rearWeightDistribution_percent');
if driverWeight < 0
    error('buildSingleRampCar:invalidSetupValue', ...
        'driverWeight_kg must be nonnegative.');
end
if rearWeightDistribution <= 0 || rearWeightDistribution >= 100
    error('buildSingleRampCar:invalidSetupValue', ...
        'rearWeightDistribution_percent must be between 0 and 100.');
end
setup.driverWeight_kg = driverWeight;
setup.rearWeightDistribution_percent = rearWeightDistribution;
setup.frontSpringRate_lb_in = scalarFinite(setup.frontSpringRate_lb_in, ...
    'frontSpringRate_lb_in');
setup.rearSpringRate_lb_in = scalarFinite(setup.rearSpringRate_lb_in, ...
    'rearSpringRate_lb_in');
setup.rearArbStiffness_NmPerRad = scalarFinite( ...
    setup.rearArbStiffness_NmPerRad,'rearArbStiffness_NmPerRad');
setup.frontRideHeight_in = scalarFinite(setup.frontRideHeight_in, ...
    'frontRideHeight_in');
setup.rearRideHeight_in = scalarFinite(setup.rearRideHeight_in, ...
    'rearRideHeight_in');
if isfield(fixed,'frontArbStiffness_NmPerRad')
    fixed.frontArbStiffness_NmPerRad = scalarFinite( ...
        fixed.frontArbStiffness_NmPerRad,'frontArbStiffness_NmPerRad');
end

map = [];
if isa(mapInput,'AeroMap')
    map = mapInput;
elseif strlength(string(mapInput)) > 0
    map = AeroMap(char(mapInput));
end
if ~isempty(map)
    if ~isstruct(assets) || ~isscalar(assets), assets = struct(); end
    assets.map = map;
end
spec = vehicleConfig.fromRampSetup(fixed,setup,map,assets);
[car,derived] = vehicleConfig.buildCar(spec,"accel",assets);
end

function value = scalarFinite(value,name)
value = double(value);
if ~isscalar(value) || ~isfinite(value)
    error('buildSingleRampCar:invalidSetupValue', ...
        '%s must be a finite scalar.',name);
end
end