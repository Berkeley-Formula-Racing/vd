%% Inputs
clear
setup_paths

weight_dist = 0.54; % rearwards
m = 535*0.453592; % mass, kg
L = 1.6; % wheelbase, m

a = L*weight_dist; % front axle to cg, m
b = L*(1-weight_dist); % rear axle to cg, m

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

% r = yaw velocity, rad/sec
% beta = sideslip angle, rad
% ay = lateral acceleration, g
% delta = steer angle, rad

%% Bundorf Cornering Compliance Effects

%E_af = % front deflection steer coefficient due to aligning torque, deg/ft-lb, positive
%E_ar = % rear deflection steer coefficient due to aligning torque, deg/ft-lb, negative


%% Transfer Functions




% yaw velocity by steer
r_delta = tf([(57.3*g*a*b*u)/(D_f*L),((57.3*g)^2*a*b)/(D_f*D_r*L)],...
    [(I_zz*u/m),(57.3*g*(a^2*b*m*D_r+a*D_f*(b^2*m+I_zz)+b*D_r*I_zz))/(D_f*D_r*m*L),...
    (57.3*g*a*b*(57.3*g*L+u^2*(D_f-D_r)))/(D_f*D_r*u*L)]);
    
% sideslip by steer
beta_delta = tf([(57.3*g*b*I_zz)/(D_f*m*L),(57.3*g*a*b*(57.3*g*b-D_r*u^2))/(D_f*D_r*u*L)],...
    [(I_zz*u/m),(57.3*g*(a^2*b*m*D_r+a*D_f*(b^2*m+I_zz)+b*D_r*I_zz))/(D_f*D_r*m*L),...
    (57.3*g*a*b*(57.3*g*L+u^2*(D_f-D_r)))/(D_f*D_r*u*L)]);

% lateral acceleration by steer

ay_delta = tf([(57.3*g*b*u*I_zz)/(D_f*m*L),((57.3*g)^2*a*b^2)/(D_f*D_r*L),...
    ((57.3*g)^2*a*b*u)/(D_f*D_r*L)],...
    [(I_zz*u/m),(57.3*g*(a^2*b*m*D_r+a*D_f*(b^2*m+I_zz)+b*D_r*I_zz))/(D_f*D_r*m*L),...
    (57.3*g*a*b*(57.3*g*L+u^2*(D_f-D_r)))/(D_f*D_r*u*L)]);

% sideslip by lateral acceleration

beta_ay = tf([I_zz*u*D_r,a*m*(57.3*g*b-D_r*u^2)],...
        [I_zz*u^2*D_r,57.3*g*a*b*u,57.3*g*a*m*u^2]);

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

%% Plots
bodeopt = bodeoptions;
bodeopt.MagScale = 'linear';
bodeopt.MagUnits = 'abs';
bodeopt.FreqScale = 'linear';
bodeopt.FreqUnits = 'Hz';
bodeopt.Title.String = 'Normalized Frequency Response';
bodeopt.Title.FontSize = 18;
bodeopt.Xlabel.FontSize = 14;
bodeopt.Ylabel.FontSize = 14;

bode(ay_delta,r_delta,beta_delta,bodeopt)
xlim([0 4]);
legend('Lateral Acceleration','Yaw Velocity','Sideslip Angle')

figure
step(r_delta)
xlabel('Time (sec)','FontSize',15)
ylabel('Yaw Velocity (rad/sec)','FontSize',15)
title('Yaw Velocity Step Response','FontSize',18)

figure
step(ay_delta)
xlabel('Time (sec)','FontSize',15)
ylabel('Lateral Acceleration (g)','FontSize',15)
title('Lateral Acceleration Step Response')

figure
step(beta_delta)
title('Sideslip Angle Step Response')

%% Time Response Computations

% response time is time required to first reach 90% of steady-state
% referenced to 50% of ramp-steer input
% it is closely approximated by rise time

r_stepinfo = stepinfo(r_delta);
ay_stepinfo = stepinfo(ay_delta);
beta_stepinfo = stepinfo(beta_delta);

r_rT = r_stepinfo.RiseTime;
ay_rT = ay_stepinfo.RiseTime;
r_tT = r_stepinfo.TransientTime;
ay_tT = ay_stepinfo.TransientTime;
r_sT = r_stepinfo.SettlingTime;
ay_sT = ay_stepinfo.SettlingTime;
r_os = r_stepinfo.Overshoot;
ay_os = ay_stepinfo.Overshoot;

% yaw velocity peak to steady state ratio
r_peak_ss_ratio = r_stepinfo.Peak/dcgain(r_delta);

%% Phase Difference (Rear Cornering Compliance Calculation)

% phase difference
freq = logspace(-2,3,1000)*(2*pi);
[~,phase] = bode(r_delta,freq);
r_phase = reshape(phase,size(freq));
semilogx(freq/(2*pi),r_phase)
hold on

[~,phase] = bode(ay_delta,freq);
ay_phase = reshape(phase,size(freq));
semilogx(freq/(2*pi),ay_phase)

% phase difference (R-Ay)
r_ay = r_phase-ay_phase;
semilogx(freq/(2*pi),r_ay)
xlabel('Frequency (Hz)','FontSize',15)
ylabel('Phase (deg)','FontSize',15)
legend('Yaw Velocity','Lateral Acceleration','Phase Difference')

% slope of R-Ay function 
d_r_ay = diff(r_ay)./diff(freq/(2*pi));

% rear cornering compliance (check)
D_r_approx = g/(2*pi*u)*d_r_ay(1)+57.3*g*b/u^2;

%% Yaw Velocity Lag

D_f_phi = D_r/(1-57.3*g*L/(D_r*u^2));

% If D_f>D_f_phi, then yaw velocity leads steering angle for low frequency
% steer inputs (negative phase). If D_f<D_f_phi, yaw velocity lags steering
% angle (increasing steering hysteresis)

%% Natural Frequency and Damping

[omega_n,zeta] = damp(ay_delta);

% understeer influences SSGs
% understeer increases natural frequency

% vehicle is unstable in open-loop if K < ackermann gain
ackermann_gain = -(57.3*g*L)/u^2;

% understeer is not independent of front/rear cornering compliances
% understeer is defined by free response and can be derived from any SSG

response_table = table(ay_rT,  ay_tT, ay_sT, ay_os, r_rT, r_tT, r_sT, r_os)