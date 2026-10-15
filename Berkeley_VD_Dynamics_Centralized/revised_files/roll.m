%%% INPUTS %%%
addpath(fullfile(fileparts(mfilename('fullpath')),'Full Car Models','config'));
sharedStudyVehicle = vdVehicleConfig();
mass = sharedStudyVehicle.studies.roll.mass; %370 lbs
g = 9.8;
r_g = sharedStudyVehicle.studies.roll.r_g; %roll gradient
k_ride_f = sharedStudyVehicle.studies.roll.k_ride_f; %n/m
k_ride_r = sharedStudyVehicle.studies.roll.k_ride_r;

cLa = sharedStudyVehicle.studies.roll.cLa;
CoP = sharedStudyVehicle.studies.roll.CoP; %percent forward
rho = sharedStudyVehicle.studies.roll.rho;

skidpad_r = 8;
skidpad_time = 5; %s

%%% OUTPUTS %%%
skidpad_vel = 2 * pi * skidpad_r / skidpad_time; %m/s
ay = ((skidpad_vel^2) / skidpad_r)/g;
roll_deg = ay*r_g;
df = rho/2*(skidpad_vel^2)*cLa;
df_f = df * CoP;
df_r = df * (1 - CoP);

comp_f = df_f / k_ride_f;
comp_r = df_r / k_ride_r;