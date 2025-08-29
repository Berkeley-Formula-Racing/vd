%% Inputs
clear;
close all;
setup_paths

weight_dist = 0.51; % rearwards
 m = 520*0.453592; % mass, kg
L = 62 * 0.0254; % wheelbase, m

load('Fy_pure_parameters_run24_new2.mat')
C_f = 0.55*2*4.448*cornering_stiffness(Xbestcell,0,535*(1-weight_dist),12,0);
C_r = 0.55*2*4.448*cornering_stiffness(Xbestcell,0,535*(weight_dist),12,0);
 
u = 10; % forward velocity, m/s

I_zz = 83.28; % yaw inertia, kg-m^2
%C_f = 136*4.448; % front cornering stiffness, N/deg @ 150 lb
%C_r = 136*4.448; % rear cornering stiffness, N/deg @ 150 lb
g = 9.806; % gravitational constant, m/sec^2

N = 4; % steering ratio


%% Outputs
a = L*weight_dist; % front axle to cg, m
b = L*(1-weight_dist); % rear axle to cg, m

D_f = 1.5*(g*m*b)/(L*C_f); % front cornering compliance, deg/g 
D_r = 1.5*(g*m*a)/(L*C_r); % rear cornering compliance, deg/g

Fyf_per_g = m * g * b / L; % N/g, front axle lateral load
Fyr_per_g = m * g * a / L; % N/g, rear axle lateral load

front_toe_vals = linspace(0.8, 0, 25);
rear_toe_vals = linspace(0.8, 0, 25);
output = [];

% Toe compliance (deg/G of lateral force)
toe_comp_front = 0;%  % deg/G
toe_comp_rear  = 0;%  % deg/G

D_f_total = D_f + toe_comp_front ;  % deg/g at individual wheel -> assume symmetrical -> deg/g at wheel = deg/g for axle
D_r_total = D_r + toe_comp_rear  ;  % deg/g

% r = yaw velocity, rad/sec
% beta = sideslip angle, rad
% ay = lateral acceleration, g
% delta = steer angle, rad

%% Bundorf Cornering Compliance Effects

%E_af = % front deflection steer coefficient due to aligning torque, deg/ft-lb, positive
%E_ar = % rear deflection steer coefficient due to aligning torque, deg/ft-lb, negative


%% Transfer Functions

% yaw velocity by steer
r_delta = tf([(57.3*g*a*b*u)/(D_f_total*L),((57.3*g)^2*a*b)/(D_f_total*D_r_total*L)],...
    [(I_zz*u/m),(57.3*g*(a^2*b*m*D_r_total+a*D_f_total*(b^2*m+I_zz)+b*D_r_total*I_zz))/(D_f_total*D_r_total*m*L),...
    (57.3*g*a*b*(57.3*g*L+u^2*(D_f_total-D_r_total)))/(D_f_total*D_r_total*u*L)]);
    
% sideslip by steer
beta_delta = tf([(57.3*g*b*I_zz)/(D_f_total*m*L),(57.3*g*a*b*(57.3*g*b-D_r_total*u^2))/(D_f_total*D_r_total*u*L)],...
    [(I_zz*u/m),(57.3*g*(a^2*b*m*D_r_total+a*D_f_total*(b^2*m+I_zz)+b*D_r_total*I_zz))/(D_f_total*D_r_total*m*L),...
    (57.3*g*a*b*(57.3*g*L+u^2*(D_f_total-D_r_total)))/(D_f_total*D_r_total*u*L)]);

% lateral acceleration by steer

ay_delta = tf([(57.3*g*b*u*I_zz)/(D_f_total*m*L),((57.3*g)^2*a*b^2)/(D_f_total*D_r_total*L),...
    ((57.3*g)^2*a*b*u)/(D_f_total*D_r_total*L)],...
    [(I_zz*u/m),(57.3*g*(a^2*b*m*D_r_total+a*D_f_total*(b^2*m+I_zz)+b*D_r_total*I_zz))/(D_f_total*D_r_total*m*L),...
    (57.3*g*a*b*(57.3*g*L+u^2*(D_f_total-D_r_total)))/(D_f_total*D_r_total*u*L)]);

% sideslip by lateral acceleration

beta_ay = tf([I_zz*u*D_r_total,a*m*(57.3*g*b-D_r_total*u^2)],...
        [I_zz*u^2*D_r_total,57.3*g*a*b*u,57.3*g*a*m*u^2]);

%% Normalization

% divide by SSG to normalize to unit gain
SSG_r_delta = dcgain(r_delta);
SSG_ay_delta = dcgain(ay_delta);
SSG_beta_delta = dcgain(beta_delta);

r_delta = r_delta/SSG_r_delta;
ay_delta = ay_delta/SSG_ay_delta;
beta_delta = beta_delta/SSG_beta_delta;
  
% control sensitivity (deg/100 deg)
steering_sensitivity = (100/N)*SSG_ay_delta/g*(pi/180); 

