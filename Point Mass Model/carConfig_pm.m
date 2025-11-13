function carCell = carConfig_pm()
%% car parameters
carParams = struct();
carParams.mass           = [155 160 168 175 180];
carParams.driver_weight  = [64];
carParams.weight_dist    = [0.512];
carParams.wheel_radius   = [0.1956];
carParams.c_roll         = [0.025];
carParams.track_width = [47] * 0.0254;

%% aero parameters
aeroParams = struct();
aeroParams.cda          = [1.48];
aeroParams.cla          = [3.97];
aeroParams.distribution = [0.418];

%% engine / powertrain
eParams = struct();
eParams.redline      = 11500;
eParams.rpm_idle     = 2000;
eParams.gears        = [32/16 30/18 28/20 26/22 24/24];  
eParams.engineSource = "KTM450";  

%% drivetrain
DTparams = struct();
DTparams.final_drive           = [33/11];
DTparams.drivetrain_efficiency = [0.97];

%% brakes
Bparams = struct();
Bparams.max_braking_torque = [850];
Bparams.tire_limited_brake = [false];

%% simple μ
muParams = struct();
muParams.track_friction_scale = [0.55];
muParams.mu_model_coeff       = [2.727];
muParams.mu_model_exp         = [-0.096];

% pass to the factorial builder
carCell = parameters_loop_pm(carParams, aeroParams, eParams, ...
                             DTparams, Bparams, muParams);
end
