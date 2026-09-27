function [car,config] = carConfigBaseline()
%CARCONFIGBASELINE Return the stable single-car ramp-speed baseline.
%   [CAR,CONFIG] = CARCONFIGBASELINE creates exactly one acceleration-oriented
%   Car object. CONFIG is the frozen, serializable source for future setup
%   duplication and editing; it is intentionally independent of carConfig.

root = fileparts(mfilename('fullpath'));
config = baselineConfig(root);
defaultSetup = config.defaultSetup;
assets = rampSpeed.loadRampAssets(config,defaultSetup);
[car,~] = buildSingleRampCar(config.fixedParameters,defaultSetup, ...
    assets.map,assets);
end

function config = baselineConfig(root)
vehicle = struct();
vehicle.mass = 162;
vehicle.accel_driver_weight = 59;
vehicle.wheelbase = 62*0.0254;
vehicle.weight_dist = 0.512;
vehicle.track_width = 47*0.0254;
vehicle.wheel_radius = 0.1956;
vehicle.cg_height = 11.75*0.0254;
vehicle.roll_center_height_front = 3.4*0.0254;
vehicle.roll_center_height_rear = 3.6*0.0254;
vehicle.I_zz = 83.28;
vehicle.ackermann = 1;
vehicle.camber_compliance_f = 0;
vehicle.camber_compliance_r = 0;
vehicle.static_r_toe = 0;
vehicle.roll_gradient_deg_per_g = 0.68;
vehicle.rear_roll_camber_outer_deg_per_deg = 0.58;
vehicle.rear_roll_camber_inner_deg_per_deg = -0.592;
vehicle.ride_camber_front_deg_per_in = 0;
vehicle.ride_camber_rear_deg_per_in = 0;
vehicle.I_wheel = 0.164;
vehicle.I_driveline = 0;
vehicle.Crr = 0.014;
vehicle.motion_ratio_front = 0.847;
vehicle.motion_ratio_rear = 0.984;
vehicle.map_reference_front_ride_height_in = 0;
vehicle.map_reference_rear_ride_height_in = 0;

vehicle.accel_cda = 0.855;
vehicle.accel_cla = 2.37;
vehicle.aero_distribution = 0.418;
vehicle.accel_cla_p_deg_p = 0;
vehicle.accel_D_p_deg_p = 0;

vehicle.redline = 11500;
vehicle.shift_point = 10000;
vehicle.gears = [32/16 30/18 28/20 26/22 24/24];
vehicle.primary_reduction = 76/32;
vehicle.torque_fn = KTM450();
vehicle.shift_time = 0.050;
vehicle.final_drive = 33/11;
vehicle.drivetrain_efficiency = 0.87;
vehicle.G_d1 = 0;
vehicle.G_d2_overrun = 0;
vehicle.G_d2_driving = 0;
vehicle.brake_distribution = 0.75;
vehicle.max_braking_torque = 840;
vehicle.gamma_f = -1;
vehicle.gamma_r = -1;

vehicle.p_i = 12;
componentRoot = fullfile(root,'carComponents');
fx = load(fullfile(componentRoot,'Fx_combined_parameters_run38_30.mat'), ...
    'Xbestcell');
fy = load(fullfile(componentRoot,'Lapsim_Fy_combined_parameters_1965run15.mat'), ...
    'Xbestcell');
vehicle.Fx_parameters = cell2mat(fx.Xbestcell);
vehicle.Fy_parameters = cell2mat(fy.Xbestcell);
vehicle.friction_scaling_factor = 1;
vehicle.grip_scaling_front = 0.6125;
vehicle.grip_scaling_rear = 0.62;

% The current acceleration setup has no front ARB; keep this fixed while the
% rear ARB remains the discrete setup factor.
fixed = struct('vehicle',vehicle,'frontArbStiffness_NmPerRad',0);
config = struct( ...
    'schemaVersion',1, ...
    'baselineVersion',"ramp-baseline-1", ...
    'id',"ramp-baseline", ...
    'label',"Ramp-speed baseline", ...
    'source',"carConfigBaseline", ...
    'fixedParameters',fixed, ...
    'defaultMapPath',fullfile(root,'aeromap_b26.csv'));

% These lists are deliberately owned here, not inherited from carConfig.
config.options = struct( ...
    'rearArbStiffness_NmPerRad',[0 699 986 1493], ...
    'frontSpringRate_lb_in',[250 300 350], ...
    'rearSpringRate_lb_in',[200 250 300]);
config.defaults = struct( ...
    'rearArbStiffness_NmPerRad',986, ...
    'frontSpringRate_lb_in',300, ...
    'rearSpringRate_lb_in',250, ...
    'frontRideHeight_in',0, ...
    'rearRideHeight_in',0, ...
    'aeroMapId',"b26", ...
    'driverWeight_kg',vehicle.accel_driver_weight, ...
    'rearWeightDistribution_percent',vehicle.weight_dist*100);
config.defaultSetup = struct( ...
    'schemaVersion',1, ...
    'baselineVersion',config.baselineVersion, ...
    'id',config.id, ...
    'label',config.label, ...
    'source',config.source, ...
    'aeroMapId',config.defaults.aeroMapId, ...
    'driverWeight_kg',config.defaults.driverWeight_kg, ...
    'rearWeightDistribution_percent',config.defaults.rearWeightDistribution_percent, ...
    'rearArbStiffness_NmPerRad',config.defaults.rearArbStiffness_NmPerRad, ...
    'frontSpringRate_lb_in',config.defaults.frontSpringRate_lb_in, ...
    'rearSpringRate_lb_in',config.defaults.rearSpringRate_lb_in, ...
    'frontRideHeight_in',config.defaults.frontRideHeight_in, ...
    'rearRideHeight_in',config.defaults.rearRideHeight_in, ...
    'isBaseline',true);
end
