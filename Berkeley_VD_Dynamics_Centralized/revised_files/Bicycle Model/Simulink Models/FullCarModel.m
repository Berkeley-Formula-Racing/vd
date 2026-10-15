%% Inputs
addpath(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'Full Car Models','config'));
sharedStudyVehicle = vdVehicleConfig();
sharedStudyTire = vdTireConfig([],[],false);
weight_dist = sharedStudyVehicle.studies.Bicycle_Model_Simulink_Models_FullCarModel.weight_dist; % rearwards
m = sharedStudyVehicle.studies.Bicycle_Model_Simulink_Models_FullCarModel.m; % mass, kg
L = sharedStudyVehicle.studies.Bicycle_Model_Simulink_Models_FullCarModel.L; % wheelbase

a = L*weight_dist; % front axle to cg, m
b = L*(1-weight_dist); % rear axle to cg, m

u = 20; % forward velocity, m/s

I_zz = sharedStudyVehicle.studies.Bicycle_Model_Simulink_Models_FullCarModel.I_zz; % yaw inertia, kg-m^2
C_f = sharedStudyTire.studies.Bicycle_Model_Simulink_Models_FullCarModel.C_f; % front cornering stiffness, N/deg @ 150 lb
C_r = sharedStudyTire.studies.Bicycle_Model_Simulink_Models_FullCarModel.C_r; % rear cornering stiffness, N/deg @ 150 lb

sharedFitData = load(sharedStudyTire.studies.Bicycle_Model_Simulink_Models_FullCarModel.fit_1,'Xbestcell');
Xbestcell = sharedFitData.Xbestcell;
X = cell2mat(Xbestcell);

track_width_f = sharedStudyVehicle.studies.Bicycle_Model_Simulink_Models_FullCarModel.track_width_f;
track_width_r = sharedStudyVehicle.studies.Bicycle_Model_Simulink_Models_FullCarModel.track_width_r;

cla = sharedStudyVehicle.studies.Bicycle_Model_Simulink_Models_FullCarModel.cla;
cp = sharedStudyVehicle.studies.Bicycle_Model_Simulink_Models_FullCarModel.cp; % proportion of downforce in front

LLTD = sharedStudyVehicle.studies.Bicycle_Model_Simulink_Models_FullCarModel.LLTD;
cg_height = sharedStudyVehicle.studies.Bicycle_Model_Simulink_Models_FullCarModel.cg_height;

step_steer = 3*pi/180; % rad

%% Plotting

sim FullCarModelSim
close_system

figure
plot(beta.time,beta.signals.values)
title('Sideslip Angle Step Response')
figure
plot(r.time,r.signals.values)
title('Yaw Velocity Step Response')
figure
plot(ay.time,ay.signals.values)
hold on
title('Lateral Acceleration Step Response')

%% Frequency Response Estimation

% Specify portion of model to estimate:
beta_io(1)=linio('FullCarModelSim/Steer Angle',1,'input');
beta_io(2)=linio('FullCarModelSim/iAngle',1,'input');
ay_io(2)=linio('FullCarModelSim/Car',3,'output');

r_io(1)=linio('FullCarModelSim/Steer Angle',1,'input');
r_io(2)=linio('FullCarModelSim/Car',4,'output');

% Specify operating point for linearization and estimation:
car_spec = operspec('FullCarModelSim');
op = findop('FullCarModelSim',car_spec);

% Linearize the model:
% sys = linearize('BicycleModelSimulink',op,io);

% Estimate the frequency response of the car model
input = frest.Sinestream('Frequency',logspace(-2,3,50));
[beta_sysest] = frestimate('FullCarModelSim',op,beta_io,input);
[ay_sysest] = frestimate('FullCarModelSim',op,ay_io,input);
[r_sysest] = frestimate('FullCarModelSim',op,r_io,input);

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

% Bode plot
bodeopt = bodeoptions;
bodeopt.MagScale = 'linear';
bodeopt.MagUnits = 'abs';
bodeopt.FreqUnits = 'Hz';

figure
bode(ay_sys,r_sys,beta_sys,bodeopt)
legend('Lateral Acceleration','Yaw Velocity','Sideslip Angle')

%% Linearization

ops_tsnapshot = linspace(0,10,50);
T = linearize('FullCarModelSim',r_io,ops_tsnapshot);
stepplot(T)

%%
% Specify portion of model to estimate:
beta_io(1)=linio('FullCarModelSim/Steer Angle',1,'input');
beta_io(2)=linio('FullCarModelSim/Car',2,'output');

ay_io(1)=linio('FullCarModelSim/Steer Angle',1,'input');
ay_io(2)=linio('FullCarModelSim/Car',3,'output');

r_io(1)=linio('FullCarModelSim/Steer Angle',1,'input');
r_io(2)=linio('FullCarModelSim/Car',4,'output');

% Specify operating point for linearization and estimation:
car_spec = operspec('FullCarModelSim');
op = findop('FullCarModelSim',car_spec);

%opspec.States(1).SteadyState = 1;
%opspec.States(1).x = beta.signals.values(end);
%opspec.States(2).SteadyState = 1;
%opspec.States(2).x = r.signals.values(end);

% Linearize the model:
r_linsys = linearize('FullCarModelSim',op,r_io);

plot(r.time,r.signals.values,'o');
hold on
opt = stepDataOptions;
opt.StepAmplitude = step_steer;
step(r_linsys,[], opt);

