%% Inputs
addpath(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'Full Car Models','config'));
sharedStudyVehicle = vdVehicleConfig();
sharedStudyTire = vdTireConfig([],[],false);
weight_dist = sharedStudyVehicle.studies.Bicycle_Model_ODE_Formulation_BicycleModelODE.weight_dist; % rearwards
m = sharedStudyVehicle.studies.Bicycle_Model_ODE_Formulation_BicycleModelODE.m; % mass, kg
L = sharedStudyVehicle.studies.Bicycle_Model_ODE_Formulation_BicycleModelODE.L; % wheelbase

u = 20; % forward velocity, m/s

I_zz = sharedStudyVehicle.studies.Bicycle_Model_ODE_Formulation_BicycleModelODE.I_zz; % yaw inertia, kg-m^2
C_f = sharedStudyTire.studies.Bicycle_Model_ODE_Formulation_BicycleModelODE.C_f; % front cornering stiffness, N/deg @ 150 lb
C_r = sharedStudyTire.studies.Bicycle_Model_ODE_Formulation_BicycleModelODE.C_r; % rear cornering stiffness, N/deg @ 150 lb
g = 9.806; % gravitational constant, m/sec^2

N = sharedStudyVehicle.studies.Bicycle_Model_ODE_Formulation_BicycleModelODE.N; % steering ratio, made up for now lol

a = L*weight_dist; % front axle to cg, m
b = L*(1-weight_dist); % rear axle to cg, m

D_f = (g*m*b)/(L*C_f); % front cornering compliance, deg/g
D_r = (g*m*a)/(L*C_r); % rear cornering compliance, deg/g
%%

sharedFitData = load(sharedStudyTire.studies.Bicycle_Model_ODE_Formulation_BicycleModelODE.fit_1,'Xbestcell');
Xbestcell = sharedFitData.Xbestcell;
Xcell = Xbestcell;

BicycleODE2 = @(t,x) BicycleODE(t,x,a,b,u,delta,C_f,C_r,m,I_zz,Xcell);

[t,x] = ode45(BicycleODE2,linspace(0,0.4,1000),[0;0]);

r = x(:,1);
beta = x(:,2);

figure
plot(r)

figure
plot(beta)
