%% Load Data

clear all;clc

% import mat file from motec
load('ConstantRadius2_processed.mat');

% input variable names of channels 
data = {engine_rpm,lat_accel,long_accel,shockpot_FR,shockpot_RR,shockpot_RL,...
    shockpot_FL,wheelspeed_FL,wheelspeed_FR,yaw_rate,steer_angle,time};

%% Constant Radius Calculations

indices = (time > 49 & time < 52)...
    | (time > 62 & time < 64)...
    | (time > 66.5 & time < 68.5)...
    | (time > 73.5 & time < 74.5)...
    | (time > 77 & time < 79)...
    | (time > 84 & time < 86)...
    | (time > 150 & time < 152)...
    | (time > 195 & time < 197)...
    | (time > 202.5 & time < 203.5);

for i = 1:numel(data)
    variable = data{i};
    data{i} = variable(indices);
end
% 
% remove_indices = abs(wheelspeed_FL)<5;
% 
% for i = 1:numel(data)
%     data{i}(remove_indices) = [];
% end

% rename variables
engine_rpm = data{1};
lat_accel = data{2};
long_accel = data{3};
shockpot_FR = data{4};
shockpot_RR = data{5};
shockpot_RL = data{6};
shockpot_FL = data{7};
wheelspeed_FL = data{8};
wheelspeed_FR = data{9};
yaw_rate = data{10};
steer_angle = data{11};
time = data{12};

% rezero
yaw_rate = yaw_rate+40;
%%

avg_index = {};
avg_index{1} = (time > 49 & time < 52);
avg_index{2} = (time > 62 & time < 64);
avg_index{3} = (time > 66.5 & time < 68.5);
avg_index{4} = (time > 73.5 & time < 74.5);
avg_index{5} = (time > 84 & time < 86);
avg_index{6} = (time > 150 & time < 152);

avg_steer_angle = [];
avg_lat_accel = [];
for i = 1:numel(avg_index)
    avg_steer_angle(i) = mean(steer_angle(avg_index{i}));
    avg_lat_accel(i) = mean(lat_accel(avg_index{i}));
end

scatter(avg_lat_accel,avg_steer_angle)

%% Plotting
close all

figure
ax1 = subplot(4,1,1);
plot(time,steer_angle)
ylim([-inf inf])
xlabel('Time(s)');
ylabel('Steer Angle (deg)');

ax2 = subplot(4,1,2);
plot(time,lat_accel)
xlabel('Time(s)');
ylabel('Lateral Acceleration (g)');
ylim([-inf inf])

ax3 = subplot(4,1,3);
plot(time,yaw_rate)
xlabel('Time(s)');
ylabel('Yaw Rate (deg/s)');
ylim([-inf inf])

ax4 = subplot(4,1,4);
plot(time,wheelspeed_FL)
hold on
plot(time,wheelspeed_FR)
xlabel('Time(s)');
ylabel('Wheelspeeds (deg/s)');
ylim([-inf inf])

linkaxes([ax1,ax2,ax3,ax4],'x')

figure
ax1 = subplot(2,1,1);
plot(time,steer_angle)
ylim([-inf inf])
xlabel('Time(s)');
ylabel('Steer Angle (deg)');

ax2 = subplot(2,1,2);
velocity = (wheelspeed_FL+wheelspeed_FR)/2*0.44704; % m/s
pseudo_lat_accel = -velocity.^2/30/9.81;
plot(time,pseudo_lat_accel);
hold on
plot(time,lat_accel);

linkaxes([ax1,ax2],'x')

figure
scatter(lat_accel,steer_angle)
    

