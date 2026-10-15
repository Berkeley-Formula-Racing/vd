%% Inputs
addpath(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'Full Car Models','config'));
sharedStudyVehicle = vdVehicleConfig();
sharedStudyTire = vdTireConfig([],[],false);
m = sharedStudyVehicle.studies.Bicycle_Model_Simulink_Models_BicycleModelSim.m; % mass, kg
L = sharedStudyVehicle.studies.Bicycle_Model_Simulink_Models_BicycleModelSim.L; % wheelbase
a = sharedStudyVehicle.studies.Bicycle_Model_Simulink_Models_BicycleModelSim.a; % front axle to cg, m
b = sharedStudyVehicle.studies.Bicycle_Model_Simulink_Models_BicycleModelSim.b; % rear axle to cg, m

u = 140/3.6; % forward velocity, m/s

I_zz = sharedStudyVehicle.studies.Bicycle_Model_Simulink_Models_BicycleModelSim.I_zz; % yaw inertia, kg-m^2
C_f = sharedStudyTire.studies.Bicycle_Model_Simulink_Models_BicycleModelSim.C_f; % front cornering stiffness, N/rad
C_r = sharedStudyTire.studies.Bicycle_Model_Simulink_Models_BicycleModelSim.C_r; % rear cornering stiffness, N/rad

%% Plotting

sim BicycleModelSimulink
close_system

figure
plot(beta.time,beta.signals.values)
figure
plot(r.time,r.signals.values)
figure
plot(ay.time,ay.signals.values)

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
sys = linearize('BicycleModelSimulink',op,io);

% Estimate the frequency response of the car model
input = frest.Sinestream('Frequency',logspace(-1,2,50));
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

data = iddata(beta.signals.values,delta.signals.values)
sys = tfest(data,2)

