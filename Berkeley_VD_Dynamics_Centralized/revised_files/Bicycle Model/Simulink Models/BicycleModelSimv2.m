%% Inputs
addpath(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'Full Car Models','config'));
sharedStudyVehicle = vdVehicleConfig();
sharedStudyTire = vdTireConfig([],[],false);
weight_dist = sharedStudyVehicle.studies.Bicycle_Model_Simulink_Models_BicycleModelSimv2.weight_dist; % rearwards
m = sharedStudyVehicle.studies.Bicycle_Model_Simulink_Models_BicycleModelSimv2.m; % mass, kg
L = sharedStudyVehicle.studies.Bicycle_Model_Simulink_Models_BicycleModelSimv2.L; % wheelbase

a = L*weight_dist; % front axle to cg, m
b = L*(1-weight_dist); % rear axle to cg, m

u = 20; % forward velocity, m/s

I_zz = sharedStudyVehicle.studies.Bicycle_Model_Simulink_Models_BicycleModelSimv2.I_zz; % yaw inertia, kg-m^2
C_f = sharedStudyTire.studies.Bicycle_Model_Simulink_Models_BicycleModelSimv2.C_f; % front cornering stiffness, N/deg @ 150 lb
C_r = sharedStudyTire.studies.Bicycle_Model_Simulink_Models_BicycleModelSimv2.C_r; % rear cornering stiffness, N/deg @ 150 lb

step_steer = 1; % rad

%% Plotting

sim BicycleModelSimulink
close_system

figure
plot(beta.time,beta.signals.values)
title('Sideslip Angle Step Response')
figure
plot(r.time,r.signals.values)
title('Yaw Velocity Step Response')
figure
plot(ay.time,ay.signals.values)
title('Lateral Acceleration Step Response')

%% Frequency Response Estimation

% Specify portion of model to estimate:
beta_io(1)=linio('BicycleModelSimulink/Steer Angle',1,'input');
beta_io(2)=linio('BicycleModelSimulink/Car',2,'output');

ay_io(1)=linio('BicycleModelSimulink/Steer Angle',1,'input');
ay_io(2)=linio('BicycleModelSimulink/Car',3,'output');

r_io(1)=linio('BicycleModelSimulink/Steer Angle',1,'input');
r_io(2)=linio('BicycleModelSimulink/Car',4,'output');

% Specify operating point for linearization and estimation:
car_spec = operspec('BicycleModelSimulink');
op = findop('BicycleModelSimulink',car_spec);

% Linearize the model:
% sys = linearize('BicycleModelSimulink',op,io);

% Estimate the frequency response of the car model
input = frest.Sinestream('Frequency',logspace(-2,3,50));
[beta_sysest] = frestimate('BicycleModelSimulink',op,beta_io,input);
[ay_sysest] = frestimate('BicycleModelSimulink',op,ay_io,input);
[r_sysest] = frestimate('BicycleModelSimulink',op,r_io,input);

% Bode plot
bodeopt = bodeoptions;
bodeopt.MagScale = 'linear';
bodeopt.MagUnits = 'abs';
bodeopt.FreqUnits = 'Hz';

figure
bode(ay_sysest,r_sysest,beta_sysest,bodeopt)
legend('Lateral Acceleration','Yaw Velocity','Sideslip Angle')

%% Transfer Function Estimation

data = iddata(beta.signals.values,delta.signals.values,0.01);
beta_sys = tfest(data,2,1);

data = iddata(ay.signals.values,delta.signals.values,0.01);
ay_sys = tfest(data,2,2);

data = iddata(r.signals.values,delta.signals.values,0.01);
r_sys = tfest(data,2,1);

% plot comparison
opt = stepDataOptions;
opt.StepAmplitude = step_steer;
step(r_sys,[], opt);
hold on
plot(r.time,r.signals.values,'o')

%% Linearization

% Specify portion of model to estimate:
beta_io(1)=linio('BicycleModelSimulink/Steer Angle',1,'input');
beta_io(2)=linio('BicycleModelSimulink/Car',2,'output');

ay_io(1)=linio('BicycleModelSimulink/Steer Angle',1,'input');
ay_io(2)=linio('BicycleModelSimulink/Car',3,'output');

r_io(1)=linio('BicycleModelSimulink/Steer Angle',1,'input');
r_io(2)=linio('BicycleModelSimulink/Car',4,'output');

% Specify operating point for linearization and estimation:
car_spec = operspec('BicycleModelSimulink');
op = findop('BicycleModelSimulink',car_spec);

% Linearize the model:
r_linsys = linearize('BicycleModelSimulink',op,r_io);
ay_linsys = linearize('BicycleModelSimulink',op,ay_io);
beta_linsys = linearize('BicycleModelSimulink',op,beta_io);

plot(r.time,r.signals.values,'o');
hold on
opt = stepDataOptions;
opt.StepAmplitude = step_steer;
step(r_linsys,[], opt);
