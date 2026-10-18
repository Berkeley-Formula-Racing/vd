%% Load Data

clear all;clc

% import mat file from motec
load('10_13_2018_Autocross1.mat');

% input variable names of channels 
data = {Engine_RPM,G_Force_Lat,G_Force_Long,Ground_Speed_Left,Ground_Speed_Right...
    Gyro_Yaw_Velocity,Steering_Angle,Drive_Speed_Left,GP_Volts_4};

load('10_13_2018_Autocross3.mat');
data2 = {Engine_RPM,G_Force_Lat,G_Force_Long,Ground_Speed_Left,Ground_Speed_Right...
    Gyro_Yaw_Velocity,Steering_Angle,Drive_Speed_Left,GP_Volts_4};
for i = 1:numel(data)
    data{i} = [data{i}; data2{i}];
end

load('10_13_2018_AutocrossJake2_processed.mat');
data2 = {Engine_RPM,G_Force_Lat,G_Force_Long,Ground_Speed_Left,Ground_Speed_Right...
    Gyro_Yaw_Velocity,Steering_Angle,Drive_Speed_Left,GP_Volts_4};
for i = 1:numel(data)
    data{i} = [data{i}; data2{i}];
end

% remove_indices = abs(wheelspeed_FL)<5;
% 
% for i = 1:numel(data)
%     data{i}(remove_indices) = [];
% end

% rename variables
engine_rpm = data{1};
lat_accel = data{2};
long_accel = data{3};
wheelspeed_FL = data{4};
wheelspeed_FR = data{5};
yaw_rate = data{6};
steer_angle = data{7};
wheelspeed_RL = data{8};
steering_voltage = data{9};

%% Plotting
figure
scatter(lat_accel,steer_angle)
xlim([-inf 0]);
ylim([-inf 0]);
xlabel('Lateral Acceleration (g)')
ylabel('Steer Angle (deg)')
title('From accelerometer')

velocity = (wheelspeed_FL+wheelspeed_FR)/2*0.44704; % m/s
pseudo_lat_accel = yaw_rate.*velocity*0.44704*pi/180/9.81; %-velocity.^2/30/9.81;
figure
scatter(pseudo_lat_accel,steer_angle)
xlim([-1.5 0]);
ylim([-inf 0]);
xlabel('Lateral Acceleration (g)')
ylabel('Steer Angle (deg)')
title('From yaw rate * velocity')

pseudo_lat_accel = -velocity.^2/30/9.81;
figure
scatter(pseudo_lat_accel,steer_angle)
xlim([-inf 0]);
ylim([-inf 0]);
xlabel('Lateral Acceleration (g)')
ylabel('Steer Angle (deg)')
title('From v^2/r')
