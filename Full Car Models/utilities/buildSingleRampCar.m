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

q = fixed.vehicle;
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
q.weight_dist = rearWeightDistribution/100;
spring = struct();
spring.front_lb_in = scalarFinite(setup.frontSpringRate_lb_in, ...
    'frontSpringRate_lb_in');
spring.rear_lb_in = scalarFinite(setup.rearSpringRate_lb_in, ...
    'rearSpringRate_lb_in');
rearArb = scalarFinite(setup.rearArbStiffness_NmPerRad, ...
    'rearArbStiffness_NmPerRad');
frontRide = scalarFinite(setup.frontRideHeight_in,'frontRideHeight_in');
rearRide = scalarFinite(setup.rearRideHeight_in,'rearRideHeight_in');

if isfield(fixed,'frontArbStiffness_NmPerRad')
    frontArb = scalarFinite(fixed.frontArbStiffness_NmPerRad, ...
        'frontArbStiffness_NmPerRad');
else
    frontArb = 0;
end

derived = rampSpeed.deriveRollStiffness(spring.front_lb_in, ...
    spring.rear_lb_in,q.motion_ratio_front,q.motion_ratio_rear, ...
    q.track_width,frontArb,rearArb);
R_sf = derived.R_sf;

map = [];
if isa(mapInput,'AeroMap')
    map = mapInput;
elseif strlength(string(mapInput)) > 0
    map = AeroMap(char(mapInput));
end
aero = Aero(q.accel_cda,q.accel_cla,q.aero_distribution, ...
    q.accel_cla_p_deg_p,q.accel_D_p_deg_p,map);

powertrain = Powertrain(q.redline,q.shift_point,q.gears, ...
    q.primary_reduction,q.torque_fn,q.shift_time,q.final_drive, ...
    q.wheel_radius,q.drivetrain_efficiency,q.G_d1,q.G_d2_overrun, ...
    q.G_d2_driving,q.brake_distribution,q.max_braking_torque);
camberRatios = [];
if isfield(assets,'camberRatios')
    camberRatios = assets.camberRatios;
end
tire = Tire2(q.p_i,q.Fx_parameters,q.Fy_parameters, ...
    q.friction_scaling_factor,camberRatios);

rideHeightAero = makeRideHeightAero(q,spring,frontRide,rearRide);
camberKinematics = struct( ...
    'roll_gradient_deg_per_g',q.roll_gradient_deg_per_g, ...
    'rear_roll_camber_outer_deg_per_deg',q.rear_roll_camber_outer_deg_per_deg, ...
    'rear_roll_camber_inner_deg_per_deg',q.rear_roll_camber_inner_deg_per_deg, ...
    'ride_camber_front_deg_per_in',q.ride_camber_front_deg_per_in, ...
    'ride_camber_rear_deg_per_in',q.ride_camber_rear_deg_per_in);
if isfield(assets,'camberModelData')
    camberKinematics.camberModelData = assets.camberModelData;
end

car = buildCarFromParameterSet(q,driverWeight,aero,powertrain, ...
    tire,rideHeightAero,camberKinematics,R_sf);

% Keep the selected suspension values visible to other model paths. The
% steady-state ramp equations use R_sf; transient paths can use these fields.
car.k_rf = frontArb;
car.k_rr = rearArb;
car.k_f_b = derived.wheelRateFront_Npm;
car.k_r_b = derived.wheelRateRear_Npm;
car.k_f_r = derived.wheelRateFront_Npm;
car.k_r_r = derived.wheelRateRear_Npm;
car.rs_total = (derived.frontRollStiffness_NmPerRad + derived.rearRollStiffness_NmPerRad)*pi/180;

derived.driverWeight_kg = driverWeight;
derived.rearWeightDistribution_percent = rearWeightDistribution;
derived.frontWeightDistribution_percent = 100 - rearWeightDistribution;
derived.totalMass_kg = car.M;
derived.staticRearLoad_N = car.M*car.g*rearWeightDistribution/100;
derived.staticFrontLoad_N = car.M*car.g*(1 - rearWeightDistribution/100);
if isfield(assets,'fingerprints')
    derived.assetFingerprints = assets.fingerprints;
end

end

function cfg = makeRideHeightAero(q,spring,frontRide,rearRide)
lbInToNpm = 175.126835;
cfg = struct( ...
    'enabled',true, ...
    'spring_rate_front_lb_in',spring.front_lb_in, ...
    'spring_rate_rear_lb_in',spring.rear_lb_in, ...
    'motion_ratio_front',q.motion_ratio_front, ...
    'motion_ratio_rear',q.motion_ratio_rear, ...
    'wheel_rate_front_Npm',spring.front_lb_in*lbInToNpm*q.motion_ratio_front^2, ...
    'wheel_rate_rear_Npm',spring.rear_lb_in*lbInToNpm*q.motion_ratio_rear^2, ...
    'static_front_ride_height_in',frontRide, ...
    'static_rear_ride_height_in',rearRide, ...
    'map_reference_front_ride_height_in',q.map_reference_front_ride_height_in, ...
    'map_reference_rear_ride_height_in',q.map_reference_rear_ride_height_in);
end

function value = scalarFinite(value,name)
value = double(value);
if ~isscalar(value) || ~isfinite(value)
    error('buildSingleRampCar:invalidSetupValue', ...
        '%s must be a finite scalar.',name);
end
end
