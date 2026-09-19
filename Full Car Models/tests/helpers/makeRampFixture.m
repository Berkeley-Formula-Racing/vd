function fixture = makeRampFixture()
%MAKERAMPFIXTURE Deterministic legacy-style data for schema tests.

fixture.frontRideHeightIn = 1.25;
fixture.rearRideHeightIn = 0.95;

fixture.cars = { ...
    struct('id',"baseline",'label',"baseline"), ...
    struct('id',"accel",'label',"acceleration")};
fixture.cases(1) = struct('id',"baseline",'label',"baseline", ...
    'source',"fixture",'designRow',1,'carRole',"lap");
fixture.cases(2) = struct('id',"accel",'label',"acceleration", ...
    'source',"fixture",'designRow',2,'carRole',"acceleration");

speeds = [5;10];
fixture.legacyRampResult = struct();
fixture.legacyRampResult.settings = struct('speeds',speeds.', ...
    'nRamp',4,'mode',"coast",'nBisect',0);

S = table(speeds,'VariableNames',{'vCar'});
S.gLat_max = [1.60;1.45];
S.gLat_top = [1.52;1.36];
S.ramp_complete = [0.95;0.94];
S.K_linear = [2.0;2.2];
S.K_r2 = [0.99;0.98];
S.K_at_limit = [2.7;2.9];
S.CUOsteerFromYaw_linear_deg = [0.010;0.012];
S.CUOsteerFromYaw_limit_deg = [0.020;0.021];
S.mech_balance = [0.52;0.53];
S.aero_balance = [0.48;0.47];
S.grip_balance_mid = [0.70;0.71];
S.grip_balance_limit = [0.30;0.29];
S.alpha_balance_mid = [0.50;0.51];
S.alpha_balance_limit = [0.48;0.47];
S.LLT_norm_balance_mid = [0.05;0.06];
S.LLT_norm_balance_limit = [0.02;0.03];
S.front_Fz_frac_mid = [0.50;0.51];
S.front_Fz_frac_limit = [0.52;0.53];
S.aero_downforce_front_N = [120;130];
S.aero_downforce_rear_N = [80;90];
S.downforce = [200;220];
S.drag = [20;25];
S.ClA = [0.02;0.025];
S.CdA = [0.018;0.020];
S.LoD = [10;8.8];
S.front_ride_height_in = [fixture.frontRideHeightIn;fixture.frontRideHeightIn];
S.rear_ride_height_in = [fixture.rearRideHeightIn;fixture.rearRideHeightIn];
S.pitch_angle_deg = [0.15;0.16];
S.front_shock_travel_in = [0.08;0.09];
S.rear_shock_travel_in = [0.02;0.03];
S.front_camber_deg = [0.10;0.12];
S.rear_camber_deg = [0.11;0.13];
S.min_Fz_limit = [-2;-1];
S.max_ceq = [0.005;0.006];
S.n_exit1 = [4;4];
S.n_exit2 = [0;0];
S.power_limited = [false;false];
S.throttle_top = [0;0];
S.lateral_metrics_applicable = [true;true];
S.traction_limited = [false;false];
fixture.legacyRampResult.perSpeed = S;

P = table([5;10],'VariableNames',{'vCar'});
P.gLat = [1;2];
P.gLong = [0.01;0.015];
P.steer_avg = [1.02;1.03];
P.lat_vel = [0.20;0.25];
P.yaw_rate = [0.10;0.12];
P.engine_rpm = [1000;1200];
P.current_gear = [2;2];
P.throttle = [0.70;0.72];
P.downforce = [100;120];
P.drag = [20;25];
P.ClA = [0.50;0.51];
P.CdA = [0.02;0.03];
P.CoP = [0.50;0.51];
P.aero_residual_in = [0.02;0.03];
P.aero_outside_map = [false;false];
P.Fz_front_axle = [500;520];
P.Fz_rear_axle = [400;410];
P.min_Fz = [200;205];
P.LLTD = [0.52;0.52];
P.LLT_front = [30;35];
P.LLT_rear = [20;22];
P.long_load_transfer = [5;6];
P.exitflag = [1;1];
P.max_ceq = [0.01;0.01];
P.lat_accel_residual = [0.01;0.01];

P.Fz_1 = [450;455]; P.Fz_2 = [440;445];
P.Fz_3 = [430;435]; P.Fz_4 = [425;430];
P.Fx_1 = [100;110]; P.Fx_2 = [100;105];
P.Fx_3 = [90;95]; P.Fx_4 = [88;92];
P.Fy_1 = [220;225]; P.Fy_2 = [215;220];
P.Fy_3 = [210;215]; P.Fy_4 = [205;210];
P.alpha_1 = [1;1]; P.alpha_2 = [1;1];
P.alpha_3 = [1;1]; P.alpha_4 = [1;1];
P.gamma_1 = [0.01;0.01]; P.gamma_2 = [0.01;0.01];
P.gamma_3 = [0.012;0.012]; P.gamma_4 = [0.012;0.012];
P.kappa_1 = [0.02;0.021]; P.kappa_2 = [0.02;0.021];
P.kappa_3 = [0.03;0.031]; P.kappa_4 = [0.03;0.031];
P.T_1 = [80;90]; P.T_2 = [80;90];
P.T_3 = [100;110]; P.T_4 = [100;110];
P.omega_1 = [120;130]; P.omega_2 = [120;130];
P.omega_3 = [130;140]; P.omega_4 = [130;140];
fixture.legacyRampResult.points = P;

longitudinalRaw = struct();
longitudinalRaw.settings = struct('speeds',speeds.');
longitudinalRaw.perSpeed = table(speeds,[2.5;2.2],[0;0],[0;0], ...
    [0;0],[0;0],[0.08;0.09],[0.08;0.09],[1;1],[0.001;0.001], ...
    [0;0],'VariableNames',{ ...
    'long_vel','long_accel','lat_accel','steer_angle','lat_vel', ...
    'yaw_rate','kappa_3','kappa_4','exitflag','max_ceq', ...
    'max_inequality_violation'});
longitudinalRaw.points = longitudinalRaw.perSpeed;
fixture.legacyLongitudinalResult = longitudinalRaw;

fixedTimestamp = datetime(2026,1,1,0,0,0);
fixture.lateralRun = rampSpeed.normalizeRampResult( ...
    fixture.legacyRampResult,"lateral",fixture.legacyRampResult.settings, ...
    fixture.cases(1),struct('source',"makeRampFixture"));
fixture.longitudinalRun = rampSpeed.normalizeRampResult( ...
    fixture.legacyLongitudinalResult,"longitudinal", ...
    fixture.legacyLongitudinalResult.settings, ...
    fixture.cases(2),struct('source',"makeRampFixture"));
fixture.lateralRun.runMeta.created = fixedTimestamp;
fixture.lateralRun.runMeta.completed = fixedTimestamp;
fixture.longitudinalRun.runMeta.created = fixedTimestamp;
fixture.longitudinalRun.runMeta.completed = fixedTimestamp;
end
