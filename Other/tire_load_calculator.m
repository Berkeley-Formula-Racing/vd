%% Calculate front tire steering angles, and center of cornering for Aero Team CFD Simulation Radius Sweep
format compact;
clc
clear all
setup_paths
carCell = carConfig(); % generate all cars to sim over
numCars = size(carCell,1);

car = carCell{1,1};
%%--------%% Input Radius %%--------%%

% initial guesses
steer_angle_guess = 0; % degrees
throttle_guess = 1;
long_vel_guess = 5; 
lat_vel_guess = 0;
yaw_rate_guess = 0.0;

kappa_1_guess = 0;
kappa_2_guess = 0;
kappa_3_guess = 0.001;
kappa_4_guess = 0.001;

x0 = [steer_angle_guess,throttle_guess,long_vel_guess,lat_vel_guess,yaw_rate_guess,kappa_1_guess,...
    kappa_2_guess,kappa_3_guess,kappa_4_guess];

% bounds
steer_angle_bounds = [0,25];
throttle_bounds = [0,1]; 
long_vel_bounds = [0,55];
lat_vel_bounds = [-3,3];
yaw_rate_bounds = [0,2];
kappa_1_bounds = [0,0];
kappa_2_bounds = [0,0];
kappa_3_bounds = [0,0.5];
kappa_4_bounds = [0,0.5];

long_vel_bounds = [long_vel_guess,long_vel_guess];

A = [];
b = [];
Aeq = [0 0 0 0 0 1 0 0 0
       0 0 0 0 0 0 1 0 0
       0 0 0 0 0 0 0 1 -1]; % symmetric rear slip: kappa_4 = kappa_3
beq = [0 0 0];
lb = [steer_angle_bounds(1),throttle_bounds(1),long_vel_bounds(1),lat_vel_bounds(1),...
    yaw_rate_bounds(1),kappa_1_bounds(1),kappa_2_bounds(1),kappa_3_bounds(1),kappa_4_bounds(1)];
ub = [steer_angle_bounds(2),throttle_bounds(2),long_vel_bounds(2),lat_vel_bounds(2),...
    yaw_rate_bounds(2),kappa_1_bounds(2),kappa_2_bounds(2),kappa_3_bounds(2),kappa_4_bounds(2)];

% scaling
scaling_factor = ones(1,9);
x0 = x0./scaling_factor;
lb = lb./scaling_factor;
ub = ub./scaling_factor;

% maximize long_decel
f = @(P) -get_long_accel(car, P);                                

% longitudinal acceleration constrained to equal zero
% velocity divided by yaw rate constrained to equal inputted radius
% used for solving skidpad (optimizing velocity for zero longitudinal acceleration
constraint = @(P) car.constraint1(P); 

% default algorithm is interior-point
options = optimoptions('fmincon','MaxFunctionEvaluations',5000,'ConstraintTolerance',1e-2,...
    'StepTolerance',1e-10,'Display','notify-detailed');

[x,fval,exitflag] = fmincon(f,x0,A,b,Aeq,beq,lb,ub,constraint,options);

[engine_rpm,beta,lat_accel,long_accel,yaw_accel,wheel_accel,omega,current_gear,...
Fzvirtual,Fz,alpha,T,Fy] = car.equations(x);

            
obj = car;
P = x;

steer_angle = P(1);
throttle = P(2); % -1 for full braking, 1 for full throttle
long_vel = P(3); % m/s
lat_vel = P(4); % m/s
yaw_rate = P(5); % equal to long_vel/radius (v/r)            
kappa = P(6:9);
% note: 1 = front left tire, 2 = front right tire
%       3 = rear left tire, 4 = rear right tire

% beta = atan(lat_vel/long_vel)*180/pi
long_vel
steer_angle
lat_accel
long_accel
throttle

[Fx,Fy,Fz] = get_all_tire_forces(obj,P);


Fx
Fy
Fz

 %% funcitons
 
 function A_y = get_long_accel(obj,P)
    [~,~,~,A_y,~,~,~,~,~,~,~,~,~] = obj.equations(P);
 end
 
 function A_x = get_lat_accel(obj,P)
    [~,~,A_x,~,~,~,~,~,~,~,~,~,~] = obj.equations(P);
 end
 
 function [Fx,Fy,Fz] = get_all_tire_forces(obj,P)           

    % inputs: vehicle parameters
    % outputs: vehicle accelerations and other properties

    % state and control matrix
    steer_angle = P(1);
    throttle = P(2); % -1 for full braking, 1 for full throttle
    long_vel = P(3); % m/s
    lat_vel = P(4); % m/s
    yaw_rate = P(5); % equal to long_vel/radius (v/r)            
    kappa = P(6:9);
    % note: 1 = front left tire, 2 = front right tire
    %       3 = rear left tire, 4 = rear right tire

    % Powertrain
    omega = zeros(1,4);
    omega(1) = (kappa(1)+1)/obj.R*(long_vel+yaw_rate*obj.t_f/2);
    omega(2) = (kappa(2)+1)/obj.R*(long_vel-yaw_rate*obj.t_f/2);
    omega(3) = (kappa(3)+1)/obj.R*(long_vel+yaw_rate*obj.t_f/2);
    omega(4) = (kappa(4)+1)/obj.R*(long_vel-yaw_rate*obj.t_f/2);

    [engine_rpm,current_gear] = obj.powertrain.engine_rpm(omega(3),omega(4),long_vel);
    [T_1,T_2,T_3,T_4] = obj.powertrain.wheel_torques(engine_rpm, omega(3), omega(4), throttle, current_gear, long_vel);
    T = [T_1,T_2,T_3,T_4];

    [Fz, Fzvirtual] = ssForces(obj,long_vel,yaw_rate,T,steer_angle*pi/180);

    % Tire Slips
    beta = atan(lat_vel/long_vel)*180/pi; % vehicle slip angle in deg
    steer_angle_1 = steer_angle; % could be modified for ackermann steering 
    steer_angle_2 = steer_angle;

    % slip angles (small angle assumption)
    alpha(1) = -steer_angle_1+(lat_vel+obj.l_f*yaw_rate)/(long_vel+yaw_rate*obj.t_f/2)*180/pi; %deg
    alpha(2) = -steer_angle_2+(lat_vel+obj.l_f*yaw_rate)/(long_vel-yaw_rate*obj.t_f/2)*180/pi; %deg
    alpha(3) = (lat_vel-obj.l_r*yaw_rate)/(long_vel+yaw_rate*obj.t_r/2)*180/pi;
    alpha(4) = (lat_vel-obj.l_r*yaw_rate)/(long_vel-yaw_rate*obj.t_r/2)*180/pi;

    % Tire Forces
    steer_angle = steer_angle_1*pi/180;
    [Fx,Fy,~] = obj.tireForce(steer_angle,alpha,kappa,Fz);
    Fx = Fx';
    Fy = Fy';
 end

 


