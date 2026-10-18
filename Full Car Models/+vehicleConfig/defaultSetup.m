function [setup,eventParams] = defaultSetup(tireModel,uncertaintyMode)
%DEFAULTSETUP Return editable QSS source defaults and event parameters.
if nargin < 1 || isempty(tireModel), tireModel = "legacy"; end
if nargin < 2 || isempty(uncertaintyMode), uncertaintyMode = "nominal"; end

modelDir = fileparts(fileparts(mfilename('fullpath')));

% Preserve the six legacy parameter groups and their original numeric values.
cP = struct();
cP.mass = 162;
cP.driver_weight = 64;
cP.accel_driver_weight = 59;
cP.wheelbase = 62*0.0254;
cP.weight_dist = 0.512;
cP.track_width = 47*0.0254;
cP.wheel_radius = 0.1956;
cP.cg_height = [11.75] * 0.0254;
cP.roll_center_height_front = 3.4 * 0.0254;
cP.roll_center_height_rear = 3.6 * 0.0254;
cP.R_sf = [0.34];
cP.I_zz = [83.28];
cP.ackermann = [1];
cP.camber_compliance_f = 0;
cP.camber_compliance_r = 0;
cP.roll_gradient_deg_per_g = 0.68;
cP.rear_roll_camber_outer_deg_per_deg = 0.58;
cP.rear_roll_camber_inner_deg_per_deg = -0.592;
cP.ride_camber_front_deg_per_in = 0;
cP.ride_camber_rear_deg_per_in = 0;
cP.static_r_toe = [0];
cP.I_wheel = 0.164;
cP.I_driveline = 0;
cP.Crr = 0.014;
cP.qss_steering_limit_deg = 25;
cP.spring_rate_front_lb_in = 300;
cP.spring_rate_rear_lb_in = 250;
cP.motion_ratio_front = 0.847;
cP.motion_ratio_rear = 0.984;
cP.static_front_ride_height_in = 0;
cP.static_rear_ride_height_in = 0;

aP = struct();
aP.cda = [1.48];
aP.cla = 3.969;
aP.cla_p_deg_p = 0;
aP.D_p_deg_p = 0;
aP.accel_cda = [0.855];
aP.accel_cla = [2.37];
aP.acc_cla_p_deg_p = 0;
aP.acc_D_p_deg_p = 0;
aP.distribution = 0.418;
aP.use_aeromap = true;
aP.map_path = fullfile(modelDir,'aeromap_b26.csv');
aP.map_reference_front_ride_height_in = 0;
aP.map_reference_rear_ride_height_in = 0;
aP.claScale = 1;
aP.cdaScale = 1;
aP.copOffset = 0;

eP = struct();
eP.redline = 11500;
eP.shift_point = 10000;
eP.gears = [32/16 30/18 28/20 26/22 24/24];
eP.primary_reduction = 76/32;
eP.torque_fn = KTM450();
eP.shift_time = 0.050;

DTp = struct();
DTp.final_drive = [33/11];
DTp.drivetrain_efficiency = [0.87];
DTp.G_d1 = 0;
DTp.G_d2_overrun = 0;
TBR = 1;
DTp.G_d2_driving = (TBR-1)./(2+2*TBR);

Bp = struct();
Bp.brake_distribution = [0.75];
Bp.max_braking_torque = 840;

tP = struct();
tP.gamma_f = -1;
tP.gamma_r = -1;
tP.p_i = [12];
tP.model_type = lower(string(tireModel));
tP.model_uncertainty = lower(string(uncertaintyMode));
tP.model_artifact = fullfile(fileparts(modelDir), ...
    'Magic Formula','experimental','lc0_nd_tire','results', ...
    'lc0_nd_pacejka_comparison.mat');
fxData = load(fullfile(modelDir,'carComponents', ...
    'Fx_combined_parameters_run38_30.mat'),'Xbestcell');
tP.Fx_parameters = cell2mat(fxData.Xbestcell);
fyData = load(fullfile(modelDir,'carComponents', ...
    'Lapsim_Fy_combined_parameters_1965run15.mat'),'Xbestcell');
tP.Fy_parameters = cell2mat(fyData.Xbestcell);
tP.friction_scaling_factor = 1;
tP.grip_scaling_front = 0.6125;
tP.grip_scaling_rear = 0.62;

setup = vehicleConfig.fromLegacyGroups(cP,aP,eP,DTp,Bp,tP);
setup.id = "qss-default";
setup.label = "QSS default";

eventParams = struct();
eventParams.skidpad_radius = 8.5;
eventParams.accel_length = 75;
eventParams.accel_rollout = 2.19;
eventParams.track_rollout = 6;
eventParams.endurance_laps = 10;
eventParams.driver_factor_autocross = 1.1;
eventParams.driver_factor_endurance = 1.085;
eventParams.winning_time = struct( ...
    'skidpad',4.8, ...
    'accel',4.206, ...
    'autocross',46.911, ...
    'endurance',1389.891);
end
