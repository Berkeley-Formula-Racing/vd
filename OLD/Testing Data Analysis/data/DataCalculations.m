%% Load Data

clear all;clc

% import mat file from motec
load('FrequencyResponse_processed.mat');

shockpot_FL = -shockpot_FL;
shockpot_RL = -shockpot_RL;

%% Inputs

shockpot_length = 2.5; % in
track_width = 50; % in
wheelbase = 60; % in
weight_dist = 0.56; % percentage rear

riderate_front = 150; %lb/in
riderate_rear = 150; %lb/in

% change to equation
MR_shock_front = 0.88;
MR_shock_rear = 1;

% change to equation
MR_pot_front = 0.65;
MR_pot_rear = 0.65;

%% Calculations

velocity = (wheelspeed_FL+wheelspeed_FR)/2*0.44704; % m/s

% divide voltage by supply voltage
shockpot_position_FL = shockpot_FL/8*shockpot_length;
shockpot_position_FR = shockpot_FR/5*shockpot_length;
shockpot_position_RL = shockpot_RL/5*shockpot_length;
shockpot_position_RR = shockpot_RR/5*shockpot_length;

% at ride height
shockpot_rideheight_FL = shockpot_position_FL(100);
shockpot_rideheight_FR = shockpot_position_FR(100);
shockpot_rideheight_RL = shockpot_position_RL(100);
shockpot_rideheight_RR = shockpot_position_RR(100);

% displacement (sign flipped since decrease in shockpot length = positive 
%   upwards bump)
shockpot_bump_FL = shockpot_rideheight_FL - shockpot_position_FL;
shockpot_bump_FR = shockpot_rideheight_FR - shockpot_position_FR;
shockpot_bump_RL = shockpot_rideheight_RL - shockpot_position_RL;
shockpot_bump_RR = shockpot_rideheight_RR - shockpot_position_RR;

% positive is bump (upwards wheel movement)
wheelbump_FL = shockpot_bump_FL/MR_pot_front;
wheelbump_FR = shockpot_bump_FR/MR_pot_front;
wheelbump_RL = shockpot_bump_RL/MR_pot_rear;
wheelbump_RR = shockpot_bump_RR/MR_pot_rear;

% shock bump
shock_bump_FL = wheelbump_FL*MR_shock_front;
shock_bump_FR = wheelbump_FR*MR_shock_front;
shock_bump_RL = wheelbump_RL*MR_shock_rear;
shock_bump_RR = wheelbump_RR*MR_shock_rear;

% vertical chassis ride
ride_displacement_front = (wheelbump_FL+wheelbump_FR)/2;
ride_displacement_rear = (wheelbump_RL+wheelbump_RR)/2;
ride_displacement = ride_displacement_front*(1-weight_dist)...
    +ride_displacement_rear*weight_dist;

% roll angle (SAE convention: positive in left hand turn - right in bump,
%   left in droop)
rollangle_front = atand((wheelbump_FR-wheelbump_FL)/track_width);
rollangle_rear = atand((wheelbump_RR-wheelbump_RL)/track_width);
rollangle = rollangle_front*(1-weight_dist)...
    +rollangle_rear*weight_dist;

% pitch angle (SAE convention: positive is nose up)
pitchangle = atand((ride_displacement_rear-ride_displacement_front)/wheelbase);

% wheel normal loads
normalload_FL = wheelbump_FL*riderate_front;
normalload_FR = wheelbump_FR*riderate_front;
normalload_RL = wheelbump_RL*riderate_rear;
normalload_RR = wheelbump_RR*riderate_rear;


%%
close all

figure
ax1 = subplot(4,1,1);
plot(time,shockpot_FL)
ylim([-inf inf])
xlabel('Time(s)');
ylabel('FL');

ax2 = subplot(4,1,2);
plot(time,shockpot_FR)
xlabel('Time(s)');
ylabel('FR');
ylim([-inf inf])

ax3 = subplot(4,1,3);
plot(time,shockpot_RL)
xlabel('Time(s)');
ylabel('RL');
ylim([-inf inf])

ax4 = subplot(4,1,4);
plot(time,shockpot_RR)
xlabel('Time(s)');
ylabel('RR');
ylim([-inf inf])

linkaxes([ax1,ax2,ax3,ax4],'x')

%%

figure
plot(time,wheelbump_FL)
hold on
plot(time,wheelbump_FR)
plot(time,wheelbump_RL)
plot(time,wheelbump_RR)
ylim([-inf inf])
xlabel('Time(s)');
ylabel('Wheel Bump');
legend('FL','FR','RL','RR');

%%

figure
scatter(lat_accel,rollangle,'.')
p = polyfit(lat_accel,rollangle,1);
y = polyval(p,linspace(-1.5,1.5,1000))
xlabel('Lateral Acceleration (g)');
ylabel('Roll Angle (deg)');
hold on
plot(linspace(-1.5,1.5,1000),y)
legend(['Roll Gradient = ' num2str(p(1))])

figure
plot(long_accel,pitchangle)
xlabel('Lateral Acceleration (g)');
ylabel('Pitch Angle (deg)');

%%
Fs = 1/(time(2)-time(1));

[pxx,f] = pwelch(shockpot_FL,[],[],[],Fs);

plot(f,log10(pxx))
xlim([0 10])
xlabel('Frequency (Hz)')
ylabel('PSD')

%%
shock_bump_FL = wheelbump_FL*MR_shock_front;
shock_bump_FR = wheelbump_FR*MR_shock_front;
shock_bump_RL = wheelbump_RL*MR_shock_rear;
shock_bump_RR = wheelbump_RR*MR_shock_rear;

shock_bump_FL_vel = diff(shock_bump_FL)/(time(2)-time(1));
shock_bump_FR_vel = diff(shock_bump_FR)/(time(2)-time(1));
shock_bump_RL_vel = diff(shock_bump_RL)/(time(2)-time(1));
shock_bump_RR_vel = diff(shock_bump_RR)/(time(2)-time(1));

histogram(shock_bump_RR_vel)
title('Damper Velocity Histogram')
xlabel('Damper Velocity (in/s)')



