%%% INPUTS %%%
mass = 370 * 0.45359; %370 lbs
g = 9.8;
r_g = 0.85; %roll gradient
k_ride_f = 23177.80762; %n/m
k_ride_r = 26147.44256;

cLa = 3.969;
CoP = 0.418; %percent forward
rho = 1.2;

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