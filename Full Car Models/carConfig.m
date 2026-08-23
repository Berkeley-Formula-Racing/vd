function [carCell,eventParams,designTable,baselineTable] = carConfig(samplingType,numSamples)
% outputs carCell:     one car per grid combination or DOE sample,
%                      column 1 the lap car, column 2 the accel car
%        eventParams:  everything the dynamic events need that is NOT a
%                      property of the car -- course geometry, lap counts,
%                      driver calibration and the scoring reference times.
%                      Events2 requires it; ask for it whenever you build one:
%                          [carCell,eventParams] = carConfig();
%                          ev = Events2(carCell{1,1},carCell{1,2},eventParams);
%                      Asking for only carCell still works, so nothing that
%                      does not build an Events2 needs changing.
%        designTable:  exact sampled parameter values, one row per car.
%        baselineTable: calibrated scalar parameter values, one row.
%
% Default is the calibrated baseline/full-factorial configuration. For DOE:
%   [cars,eventParams,X] = carConfig("Explicit",overrideTable);
if nargin < 1 || isempty(samplingType), samplingType = "FullFactorial"; end
if nargin < 2, numSamples = []; end

% car parameters (updated 2/4/21)
carParams = struct();
carParams.mass = 162; % not including driver (366 lb)
carParams.driver_weight = 64; %
carParams.accel_driver_weight = 59; % (130 lb)
carParams.wheelbase = 62*0.0254; % 62 in
carParams.weight_dist = 0.512; % percentage of weight in rear
carParams.track_width = 47*0.0254; % (47 in)
carParams.wheel_radius = 0.1956; % loaded
% radius (7.7 in)
carParams.cg_height = [11.75] * 0.0254; % (12 in) % 0.2965
carParams.roll_center_height_front = 3.4 * 0.0254; %
carParams.roll_center_height_rear = 3.6 * 0.0254; %
carParams.R_sf = [0.34]; % proportion of roll stiffness in front (not same as LLTD)
carParams.I_zz = [83.28];%, 82.28]; %kg-m^2
carParams.ackermann = [1]; %expressed as exponent for current ackermann curve
carParams.camber_compliance_f =  0; %lateral deg/G
carParams.camber_compliance_r =  0;
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

% KTM engine parameters (updated 5/1/19)
eParams = struct();
eParams.redline = 11500; % 11500   20000
eParams.shift_point = 10000; % approximate 10000   25000
% these parameters are non-iterable
eParams.gears = [32/16 30/18 28/20 26/22 24/24]; % updated KTM450[32/16 30/18 28/20 26/22 24/24]
eParams.primary_reduction = 76/32; % KTM450 76/32
eParams.torque_fn = KTM450(); %KTM450()
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

% tire parameters (updated 5/1/19)
tireParams = struct();
tireParams.gamma_f = -1; %linspace(0, -1.5, 8); % camber angle
tireParams.gamma_r = -1; %linspace(0, -1.5, 8); % camber angle

tireParams.p_i = [12]; % pressure
% these parameters are non-iterable
load('Fx_combined_parameters_run38_30.mat'); % F_x combined magic formula parameters
tireParams.Fx_parameters = cell2mat(Xbestcell);
load('Lapsim_Fy_combined_parameters_1965run15.mat'); % F_y combined magic formula parameters
tireParams.Fy_parameters = cell2mat(Xbestcell);
tireParams.friction_scaling_factor = 1; % scales tire forces to account for test/road surface difference
tireParams.grip_scaling_front = 0.6125;
tireParams.grip_scaling_rear  = 0.62;

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

%% calibrated scalar baseline and requested cases
[baselineCars,baselineTable] = parameters_loop(carParams,aeroParams,eParams, ...
    DTparams,Bparams,tireParams,"FullFactorial",[]);

request = lower(string(samplingType));
switch request
    case {"fullfactorial","grid"}
        carCell = baselineCars;
        designTable = baselineTable;

    case "explicit"
        [carCell,designTable] = parameters_loop(carParams,aeroParams,eParams, ...
            DTparams,Bparams,tireParams,"Explicit",numSamples);

    case {"lhs","random"}
        warning('carConfig:legacyDOE', ...
            ['carConfig legacy DOE sampling is deprecated; use Explicit ' ...
             'with a resolved DOE design table.']);
        study = DOEStudyConfig();
        resolved = doeResolveStudy(study,baselineTable);
        [~,overrides] = doeInitialDesign(resolved,numSamples,study.randomSeed,request);
        [carCell,designTable] = parameters_loop(carParams,aeroParams,eParams, ...
            DTparams,Bparams,tireParams,"Explicit",overrides);

    otherwise
        error('carConfig:badRequest', ...
            'samplingType must be FullFactorial, Explicit, LHS, or Random.')
end
