function cfg = vdVehicleConfig()
%VDVEHICLECONFIG Master defaults for the repository's existing full-car model.
% Values retain the upstream calibration. Units are documented at each field.
% weight_dist is the REAR mass fraction. These are not SN6 specifications.
% Study overrides are applied after loading this configuration.
fullCarRoot = fileparts(fileparts(mfilename('fullpath')));

% car parameters (updated 2/4/21)
carParams = struct();
carParams.mass = 162; % kg, excluding driver
carParams.driver_weight = 64; % kg, lap driver
carParams.accel_driver_weight = 59; % kg, acceleration driver
carParams.wheelbase = 62*0.0254; % 62 in
carParams.weight_dist = 0.512; % percentage of weight in rear
carParams.track_width = 47*0.0254; % (47 in)
carParams.wheel_radius = 0.1956; % loaded
% radius (7.7 in)
carParams.cg_height = [11.75] * 0.0254; % m, converted from 11.75 in
carParams.roll_center_height_front = 3.4 * 0.0254; %
carParams.roll_center_height_rear = 3.6 * 0.0254; %
carParams.R_sf = [0.34]; % proportion of roll stiffness in front (not same as LLTD)
carParams.I_zz = [83.28];%, 82.28]; %kg-m^2
carParams.ackermann = [1]; %expressed as exponent for current ackermann curve
carParams.camber_compliance_f =  0; %lateral deg/G
carParams.camber_compliance_r =  0;
% Camber kinematics. The front roll/steer surface remains the fitted TTC
% model; the rear roll gains below were formerly hidden in Camber_Evaluation.
% Ride-camber gains are signed deg/in of wheel bump compression: negative
% increases negative physical camber under bump. Leave at zero until the
% suspension kinematic model or measured data supplies them.
carParams.roll_gradient_deg_per_g = 0.68;
carParams.rear_roll_camber_outer_deg_per_deg = 0.58;
carParams.rear_roll_camber_inner_deg_per_deg = -0.592;
carParams.ride_camber_front_deg_per_in = 0;
carParams.ride_camber_rear_deg_per_in  = 0;
carParams.static_r_toe = [0]; %toe in deg, toe out - negative

% Rotating inertia (kg-m^2).
%
% I_wheel is PER WHEEL and covers tyre, rim, hub and brake disc. I_driveline
% is referred to the CRANK -- crank, balancer, clutch basket and primary gear
% -- and gets multiplied by the square of the total reduction, so it matters
% far more in first gear than in fifth.
%
carParams.I_wheel = 0.164;

% I_driveline is now ZERO, and that is a deliberate choice, not "no crank
% inertia".
carParams.I_driveline = 0;

% Rolling resistance coefficient. Force = Crr*(M*g + downforce) opposing
% travel, added in equations() alongside aero drag.
%
% 0.014 is FITTED to the 75 m accel (4.80 s) with I_driveline held at 0, not
% measured. It is a plausible warm-slick value (FSAE run 0.012-0.030), but read
% it as the whole "extra longitudinal resistance" the accel needs, absorbing
% rolling resistance, driveline inertia and whatever else.

carParams.Crr = 0.014;
% Shared QSS steering stop.  Keep this at the existing 25 degree bound until
% a measured physical steering stop is supplied.
carParams.qss_steering_limit_deg = 25;

% Ride-height model. Rates are spring rates at the damper, with motion ratio
% defined as shock travel / wheel travel. The map looks up height changes in
% inches relative to its CFD reference, so static and map-reference heights
% both start at zero: the first aero query is the map's (0,0) row.
carParams.spring_rate_front_lb_in = 300;
carParams.spring_rate_rear_lb_in  = 250;
carParams.motion_ratio_front = 0.847;
carParams.motion_ratio_rear  = 0.984;
% Keep the default car as one calibrated baseline.  Ride-height variation is
% generated explicitly by aeroMapStarCases rather than multiplying every
% ordinary carConfig/DOE call into a 3-by-3 full-factorial grid.
carParams.static_front_ride_height_in = 0;
carParams.static_rear_ride_height_in  = 0;

% aero parameters (updated 6/6/22)
aeroParams = struct();
aeroParams.cda = [1.48]; % m^2 (1.88)
aeroParams.cla = 3.969; % m^2 (3.45)

aeroParams.cla_p_deg_p = 0;
aeroParams.D_p_deg_p = 0;

aeroParams.accel_cda = [0.855]; % low drag
aeroParams.accel_cla = [2.37]; %
aeroParams.acc_cla_p_deg_p = 0; % NOT USED - see note above
aeroParams.acc_D_p_deg_p = 0;   % NOT USED - see note above

aeroParams.distribution = 0.418; % proportion of downforce in front
% Main aero-model switch. true uses the ride-height coupled aeromap below;
% false uses the static cla/cda/distribution values immediately above.
aeroParams.use_aeromap = true;
aeroParams.map_path = fullfile(fullCarRoot,'aeromap_b26.csv');
% These are independent of the car's static ride heights. Set them to the
% CAD F/R heights used for CFD once known; leave both at zero to make the
% (0,0) aeromap point the static baseline.
aeroParams.map_reference_front_ride_height_in = 0;
aeroParams.map_reference_rear_ride_height_in  = 0;

% KTM engine parameters (updated 5/1/19)
eParams = struct();
eParams.redline = 11500; % 11500   20000
eParams.shift_point = 10000; % approximate 10000   25000
% these parameters are non-iterable
eParams.gears = [32/16 30/18 28/20 26/22 24/24]; % updated KTM450[32/16 30/18 28/20 26/22 24/24]
eParams.primary_reduction = 76/32; % KTM450 76/32
% Evaluate this in the builders, so geometry-only studies need no engine path.
eParams.torque_fn = @KTM450; %KTM450()
eParams.shift_time = 0.050; % seconds FOR UPSHIFT ONLY; 150ms for downshift

% drivetrain parameters (updated 10/14/23)
DTparams = struct();
DTparams.final_drive = [33/11];% drivetrain sprocket ratio [33/11] 7.2918
DTparams.drivetrain_efficiency = [0.87]; % scales torque value  (0.87)
DTparams.G_d1 = 0; % differential torque transfer offset due to internal friction
DTparams.G_d2_overrun = 0; % differential torque transfer gain in overrun (not used right now)
TBR = 1;%1:0.5:4;
DTparams.G_d2_driving = (TBR-1)./(2+2*TBR); % differential torque transfer gain on power

% brake parameters (updated 8/14/23)
Bparams = struct();
Bparams.brake_distribution = [0.75];% proportion of brake torque applied to front
Bparams.max_braking_torque = 840; % total braking torque (Nm)

%% ---------------- event parameters ----------------

eventParams = struct();

% Course geometry and event format
eventParams.skidpad_radius  = 8.5;   % m, PATH radius the car drives. The inner
                                     % cone circle is 15.25 m across; 8.5 leaves
                                     % roughly a foot of margin to the cones.
eventParams.accel_length    = 75;    % m, timed acceleration run
eventParams.accel_rollout   = 2.19;  % m, run-up before the timing line on accel
eventParams.track_rollout   = 6;     % m, run-up before the line on autocross
                                     % and endurance
eventParams.endurance_laps  = 10;    % laps in the endurance event

eventParams.driver_factor_autocross = 1.1;
eventParams.driver_factor_endurance = 1.085;

% Scoring reference times: the fastest time recorded at the event, which the
% FSAE points formulas score against. Update these per event, not per car.
eventParams.winning_time = struct( ...
    'skidpad',   4.8, ...
    'accel',     4.206, ...     % Michigan 2024 (2023 was 4.174)
    'autocross', 46.911, ...    % Michigan 2024 (2023 was 45.886)
    'endurance', 1389.891);


cfg = struct('car',carParams,'aero',aeroParams,'engine',eParams, ...
    'drivetrain',DTparams,'brakes',Bparams,'events',eventParams);

% Ramp-speed setup choices. Shared vehicle and tire values come from the
% same defaults as carConfig; these are the existing ramp-specific choices.
cfg.ramp.cad_front_x_in = 4.378;
cfg.ramp.cad_rear_x_in = 59.75;
cfg.ramp.cad_front_z_in = 4.0;
cfg.ramp.cad_rear_z_in = 5.7;
cfg.ramp.frontArbStiffness_NmPerRad = 0;
cfg.ramp.rearArbStiffness_NmPerRad = 986;
cfg.ramp.rearArbOptions_NmPerRad = [0 699 986 1493];
cfg.ramp.frontSpringOptions_lb_in = [250 300 350];
cfg.ramp.rearSpringOptions_lb_in = [200 250 300];

% Historical suspension studies used different setups. Keep the choices
% explicit here instead of silently replacing one study's values with another.
s = struct('k',250*4.45*39.37,'Ixx',60,'Iyy',82, ...
    'k_rf',3300*6.474,'k_rr',3200*6.474,'damper_index',9, ...
    'motion_ratio_front','MR18Front','motion_ratio_rear','MR18Front');
cfg.suspension.modal = s;
s.k_rf = 3000*6.474; s.k_rr = 7000*6.474; s.damper_index = 6;
cfg.suspension.decoupled = s;
s.k = 200*4.45*39.37;
cfg.suspension.transient = s;
cfg.suspension.vibration = struct('k',s.k,'c',10000, ...
    'Ixx',s.Ixx,'Iyy',s.Iyy,'k_rf',18000,'k_rr',18000);
cfg.suspension.damping_file = fullfile(fullCarRoot, ...
    'Suspension Transient Model','dampingcurves.mat');

% Standalone historical studies retain their original inputs and native units.
% They are separate named presets, not a claim that they describe one car.
% Bicycle Model/toe_compliance.m
cfg.studies.Bicycle_Model_toe_compliance.weight_dist = 0.51;
cfg.studies.Bicycle_Model_toe_compliance.m = 520*0.453592;
cfg.studies.Bicycle_Model_toe_compliance.L = 62 * 0.0254;
cfg.studies.Bicycle_Model_toe_compliance.I_zz = cfg.car.I_zz;
cfg.studies.Bicycle_Model_toe_compliance.N = 4;
% Bicycle Model/BicycleModel.m
cfg.studies.Bicycle_Model_BicycleModel.weight_dist = 0.54;
cfg.studies.Bicycle_Model_BicycleModel.m = 535*0.453592;
cfg.studies.Bicycle_Model_BicycleModel.L = 1.6;
cfg.studies.Bicycle_Model_BicycleModel.I_zz = cfg.car.I_zz;
cfg.studies.Bicycle_Model_BicycleModel.N = 4;
% Bicycle Model/YawResponseTimePlotter.m
cfg.studies.Bicycle_Model_YawResponseTimePlotter.weight_dist = 0.51;
cfg.studies.Bicycle_Model_YawResponseTimePlotter.m = 535*0.453592;
cfg.studies.Bicycle_Model_YawResponseTimePlotter.L = 1.524;
cfg.studies.Bicycle_Model_YawResponseTimePlotter.I_zz = cfg.car.I_zz;
cfg.studies.Bicycle_Model_YawResponseTimePlotter.N = 4;
% Bicycle Model/LatAccelResponseTimePlotter.m
cfg.studies.Bicycle_Model_LatAccelResponseTimePlotter.weight_dist = 0.51;
cfg.studies.Bicycle_Model_LatAccelResponseTimePlotter.m = 535*0.453592;
cfg.studies.Bicycle_Model_LatAccelResponseTimePlotter.L = 1.524;
cfg.studies.Bicycle_Model_LatAccelResponseTimePlotter.I_zz = cfg.car.I_zz;
cfg.studies.Bicycle_Model_LatAccelResponseTimePlotter.N = 4;
% Bicycle Model/Simulink Models/BicycleModelSimv2.m
cfg.studies.Bicycle_Model_Simulink_Models_BicycleModelSimv2.weight_dist = 0.3;
cfg.studies.Bicycle_Model_Simulink_Models_BicycleModelSimv2.m = 242.6718;
cfg.studies.Bicycle_Model_Simulink_Models_BicycleModelSimv2.L = 1.524;
cfg.studies.Bicycle_Model_Simulink_Models_BicycleModelSimv2.I_zz = cfg.car.I_zz;
% Bicycle Model/Simulink Models/NonLin_BicycleModel.m
cfg.studies.Bicycle_Model_Simulink_Models_NonLin_BicycleModel.weight_dist = 0.55;
cfg.studies.Bicycle_Model_Simulink_Models_NonLin_BicycleModel.m = 242.6718;
cfg.studies.Bicycle_Model_Simulink_Models_NonLin_BicycleModel.L = 1.524;
cfg.studies.Bicycle_Model_Simulink_Models_NonLin_BicycleModel.I_zz = cfg.car.I_zz;
% Bicycle Model/Simulink Models/FullCarModel.m
cfg.studies.Bicycle_Model_Simulink_Models_FullCarModel.weight_dist = 0.6;
cfg.studies.Bicycle_Model_Simulink_Models_FullCarModel.m = 242.6718;
cfg.studies.Bicycle_Model_Simulink_Models_FullCarModel.L = 1.524;
cfg.studies.Bicycle_Model_Simulink_Models_FullCarModel.I_zz = cfg.car.I_zz;
cfg.studies.Bicycle_Model_Simulink_Models_FullCarModel.track_width_f = 1.1938;
cfg.studies.Bicycle_Model_Simulink_Models_FullCarModel.track_width_r = 1.1938;
cfg.studies.Bicycle_Model_Simulink_Models_FullCarModel.cla = 1.771;
cfg.studies.Bicycle_Model_Simulink_Models_FullCarModel.cp = 0.5;
cfg.studies.Bicycle_Model_Simulink_Models_FullCarModel.LLTD = 0.6;
cfg.studies.Bicycle_Model_Simulink_Models_FullCarModel.cg_height = 0.3048;
% Bicycle Model/Simulink Models/BicycleModelSim.m
cfg.studies.Bicycle_Model_Simulink_Models_BicycleModelSim.m = 1500;
cfg.studies.Bicycle_Model_Simulink_Models_BicycleModelSim.L = 2500;
cfg.studies.Bicycle_Model_Simulink_Models_BicycleModelSim.I_zz = 2500;
% Bicycle Model/ODE Formulation/BicycleModelODE.m
cfg.studies.Bicycle_Model_ODE_Formulation_BicycleModelODE.weight_dist = 0.3;
cfg.studies.Bicycle_Model_ODE_Formulation_BicycleModelODE.m = 242.6718;
cfg.studies.Bicycle_Model_ODE_Formulation_BicycleModelODE.L = 1.524;
cfg.studies.Bicycle_Model_ODE_Formulation_BicycleModelODE.I_zz = cfg.car.I_zz;
cfg.studies.Bicycle_Model_ODE_Formulation_BicycleModelODE.N = 16;
% roll.m
cfg.studies.roll.mass = 370 * 0.45359;
cfg.studies.roll.r_g = 0.85;
cfg.studies.roll.k_ride_f = 23177.80762;
cfg.studies.roll.k_ride_r = 26147.44256;
cfg.studies.roll.cLa = 3.969;
cfg.studies.roll.CoP = 0.418;
cfg.studies.roll.rho = 1.2;
% understeer_bode.m
cfg.studies.understeer_bode.L = 62 * 0.0254;
cfg.studies.understeer_bode.m = 535 * 0.453592;
cfg.studies.understeer_bode.Iz = cfg.car.I_zz;
% understeer_budget.m
cfg.studies.understeer_budget.m = 168;
cfg.studies.understeer_budget.L = 62 * 25.4;
cfg.studies.understeer_budget.weight_dist = 0.51;
cfg.studies.understeer_budget.WUF = 17.07;
cfg.studies.understeer_budget.EF = 0.001;
cfg.studies.understeer_budget.GF = 0.789;
cfg.studies.understeer_budget.KF = 2288;
cfg.studies.understeer_budget.EYF = 0.191;
cfg.studies.understeer_budget.GYF = 0.291;
cfg.studies.understeer_budget.ENF = 0.774;
cfg.studies.understeer_budget.GNF = 0.022;
cfg.studies.understeer_budget.ROLLpg = 4.47;
% calcWheelForcesAndDisplacements_pitch.m
cfg.studies.calcWheelForcesAndDisplacements_pitch.front_anti_lift = 0;
cfg.studies.calcWheelForcesAndDisplacements_pitch.front_anti_dive = .545;
cfg.studies.calcWheelForcesAndDisplacements_pitch.rear_anti_lift = 1.02;
cfg.studies.calcWheelForcesAndDisplacements_pitch.rear_anti_squat = .457;
% Full Car Models/CamberCurveGen/generateIdealCamberCurves.m
cfg.studies.camber1_generateIdealCamberCurves.C_front_width = cfg.car.track_width/0.0254;
cfg.studies.camber1_generateIdealCamberCurves.C_rear_width = cfg.car.track_width/0.0254;
cfg.studies.camber1_generateIdealCamberCurves.C_front_spring_roll_stiffness = 3000;
cfg.studies.camber1_generateIdealCamberCurves.C_rear_spring_roll_stiffness = 3000;
cfg.studies.camber1_generateIdealCamberCurves.C_front_ARB_roll_stiffness = 0;
cfg.studies.camber1_generateIdealCamberCurves.C_rear_ARB_roll_stiffness = 0;
cfg.studies.camber1_generateIdealCamberCurves.C_weight_dist = 0.54;
cfg.studies.camber1_generateIdealCamberCurves.C_mass = (395 + 150);
% Full Car Models/CamberCurveGen/Testcamber.m
cfg.studies.camber1_Testcamber.C_front_width = cfg.car.track_width/0.0254;
cfg.studies.camber1_Testcamber.C_rear_width = cfg.car.track_width/0.0254;
cfg.studies.camber1_Testcamber.C_wheelbase = cfg.car.wheelbase/0.0254;
cfg.studies.camber1_Testcamber.C_front_spring_roll_stiffness = 2860;
cfg.studies.camber1_Testcamber.C_rear_spring_roll_stiffness = 2911;
cfg.studies.camber1_Testcamber.C_front_bump_stiffness = 2*191.66;
cfg.studies.camber1_Testcamber.C_rear_bump_stiffness = 2*196.02;
cfg.studies.camber1_Testcamber.C_front_ARB_roll_stiffness = 0;
cfg.studies.camber1_Testcamber.C_rear_ARB_roll_stiffness = 986;
cfg.studies.camber1_Testcamber.C_weight_dist = 0.54;
cfg.studies.camber1_Testcamber.C_mass = (395 + 150)*4.44;
cfg.studies.camber1_Testcamber.front_anti_dive = .545;
cfg.studies.camber1_Testcamber.rear_anti_lift = 1.02;
cfg.studies.camber1_Testcamber.rear_anti_squat = .457;
% Full Car Models/CamberCurveGen/generateIdealCamberCurves_simplified.m
cfg.studies.camber1_generateIdealCamberCurves_simplified.C_front_width = cfg.car.track_width/0.0254;
cfg.studies.camber1_generateIdealCamberCurves_simplified.C_rear_width = cfg.car.track_width/0.0254;
cfg.studies.camber1_generateIdealCamberCurves_simplified.C_front_spring_roll_stiffness = 3000;
cfg.studies.camber1_generateIdealCamberCurves_simplified.C_rear_spring_roll_stiffness = 0;
cfg.studies.camber1_generateIdealCamberCurves_simplified.C_front_ARB_roll_stiffness = 0;
cfg.studies.camber1_generateIdealCamberCurves_simplified.C_rear_ARB_roll_stiffness = 0;
cfg.studies.camber1_generateIdealCamberCurves_simplified.C_weight_dist = 0.54;
cfg.studies.camber1_generateIdealCamberCurves_simplified.C_mass = (395 + 150);
% Full Car Models/CamberCurveGen2/generateIdealCamberCurves.m
cfg.studies.camber2_generateIdealCamberCurves.C_front_width = cfg.car.track_width/0.0254;
cfg.studies.camber2_generateIdealCamberCurves.C_rear_width = cfg.car.track_width/0.0254;
cfg.studies.camber2_generateIdealCamberCurves.C_front_spring_roll_stiffness = 3000;
cfg.studies.camber2_generateIdealCamberCurves.C_rear_spring_roll_stiffness = 3000;
cfg.studies.camber2_generateIdealCamberCurves.C_front_ARB_roll_stiffness = 0;
cfg.studies.camber2_generateIdealCamberCurves.C_rear_ARB_roll_stiffness = 0;
cfg.studies.camber2_generateIdealCamberCurves.C_weight_dist = 0.54;
cfg.studies.camber2_generateIdealCamberCurves.C_mass = (395 + 150);
% Full Car Models/CamberCurveGen2/Testcamber.m
cfg.studies.camber2_Testcamber.C_front_width = cfg.car.track_width/0.0254;
cfg.studies.camber2_Testcamber.C_rear_width = cfg.car.track_width/0.0254;
cfg.studies.camber2_Testcamber.C_front_spring_roll_stiffness = 2860;
cfg.studies.camber2_Testcamber.C_rear_spring_roll_stiffness = 2911;
cfg.studies.camber2_Testcamber.C_front_ARB_roll_stiffness = 0;
cfg.studies.camber2_Testcamber.C_rear_ARB_roll_stiffness = 986;
cfg.studies.camber2_Testcamber.C_weight_dist = 0.54;
cfg.studies.camber2_Testcamber.C_mass = 395 + 150;
% Full Car Models/CamberCurveGen2/generateIdealCamberCurves_simplified.m
cfg.studies.camber2_generateIdealCamberCurves_simplified.C_front_width = cfg.car.track_width/0.0254;
cfg.studies.camber2_generateIdealCamberCurves_simplified.C_rear_width = cfg.car.track_width/0.0254;
cfg.studies.camber2_generateIdealCamberCurves_simplified.C_front_spring_roll_stiffness = 3000;
cfg.studies.camber2_generateIdealCamberCurves_simplified.C_rear_spring_roll_stiffness = 0;
cfg.studies.camber2_generateIdealCamberCurves_simplified.C_front_ARB_roll_stiffness = 0;
cfg.studies.camber2_generateIdealCamberCurves_simplified.C_rear_ARB_roll_stiffness = 0;
cfg.studies.camber2_generateIdealCamberCurves_simplified.C_weight_dist = 0.54;
cfg.studies.camber2_generateIdealCamberCurves_simplified.C_mass = (395 + 150);

% Remaining paired and embedded inputs from the same historical studies.
cfg.studies.understeer_budget.FYRperg = 8430.9;
cfg.studies.understeer_budget.NFperg = -403.16;
cfg.studies.understeer_budget.NRperg = -317.93;
cfg.studies.understeer_budget.WR = 885.000;
cfg.studies.understeer_budget.WUR = 99.000;
cfg.studies.understeer_budget.ER = 0.055;
cfg.studies.understeer_budget.GR = -0.720;
cfg.studies.understeer_budget.KR = 367;
cfg.studies.understeer_budget.EYR = 0.045;
cfg.studies.understeer_budget.GYR = -0.328;
cfg.studies.understeer_budget.ENR = -0.113;
cfg.studies.understeer_budget.GNR = 0.023;
cfg.studies.Bicycle_Model_Simulink_Models_BicycleModelSim.a = 1.1;
cfg.studies.Bicycle_Model_Simulink_Models_BicycleModelSim.b = 1.6;
cfg.studies.understeer_bode.weight_dist = 0.512;

end
