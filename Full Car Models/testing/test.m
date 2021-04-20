load('fmincon inputs.mat');
long_vel_guess = C.long_vel_guess;
lat_accel_value = C.lat_accel_value;
car = C.car;


% initial guesses
steer_angle_guess = -2;%1
throttle_guess = 1;%0.1
lat_vel_guess = -0.1;
yaw_rate_guess = lat_accel_value/long_vel_guess;

kappa_1_guess = 0;
kappa_2_guess = 0;
kappa_3_guess = 0.01;
kappa_4_guess = 0.01;

x0 = [steer_angle_guess;
    throttle_guess;
    long_vel_guess;
    lat_vel_guess;
    yaw_rate_guess;
    kappa_1_guess;
    kappa_2_guess;
    kappa_3_guess;
    kappa_4_guess]';


x0(3) = long_vel_guess;

% bounds
steer_angle_bounds = [0,25];
throttle_bounds = [0,1];
long_vel_bounds = [long_vel_guess,long_vel_guess];
lat_vel_bounds = [-3,3];
yaw_rate_bounds = [lat_accel_value/long_vel_guess,lat_accel_value/long_vel_guess];
kappa_1_bounds = [0,0];
kappa_2_bounds = [0,0];
kappa_3_bounds = [0,0.2];
kappa_4_bounds = [0,0.2];

A = [];
b = [];
Aeq = [0 0 1 0 0 0 0 0 0
       0 0 0 0 0 1 0 0 0
       0 0 0 0 0 0 1 0 0];
beq = [long_vel_guess 0 0];
lb = [steer_angle_bounds(1),throttle_bounds(1),long_vel_bounds(1),lat_vel_bounds(1),...
    yaw_rate_bounds(1),kappa_1_bounds(1),kappa_2_bounds(1),kappa_3_bounds(1),kappa_4_bounds(1)];
ub = [steer_angle_bounds(2),throttle_bounds(2),long_vel_bounds(2),lat_vel_bounds(2),...
    yaw_rate_bounds(2),kappa_1_bounds(2),kappa_2_bounds(2),kappa_3_bounds(2),kappa_4_bounds(2)];

% scaling
%scaling_factor = [20 1 10 1 1 0.01 0.01 0.01 0.01];

% objective function: longitudinal acceleration (forwards)
f = @(P) -car.long_accel(P);

% constrained to lateral acceleration value
constraint = @(P) car.constraint4(P,lat_accel_value);

% default algorithm is interior-point

%options = setOptimoptions(1000);
options = optimoptions('fmincon','MaxFunctionEvaluations',5000,'ConstraintTolerance',1e-2,...
    'StepTolerance',1e-10,'Display','notify-detailed');
% fval: objective function value (v^2/r) 
[x,fval,exitflag] = fmincon(f,x0,A,b,Aeq,beq,lb,ub,constraint,options);

[engine_rpm,beta,lat_accel,long_accel,yaw_accel,wheel_accel,omega,current_gear,...
Fzvirtual,Fz,alpha,T] = car.equations(x);

long_accel_guess = x;

% generate vector of control variable values
x_accel = [exitflag long_accel lat_accel x omega(1:4) engine_rpm current_gear beta...
    Fz(1:4) alpha(1:4) T(1:4)]

long_accel = long_accel;
