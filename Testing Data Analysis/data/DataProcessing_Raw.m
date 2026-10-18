%% Load Data

clear all;clc

% import mat file from motec
load('18in_lowrideheight.mat');

% input variable names of channels 
data = {Engine_RPM,G_Force_Lat,G_Force_Long,Ground_Speed_Left,Ground_Speed_Right...
    Gyro_Yaw_Velocity,Steering_Angle,Drive_Speed_Left,GP_Volts_4};

%%
for i = 1:numel(data)
    data_filtered{i} = data{i}.Value;
    time_cell{i} = data{i}.Time;
end

% find minimum sample size
[minsize, minindex] = min(cellfun('size', data_filtered, 2));
time = time_cell{minindex};

for i = 1:numel(data)
    % construct timeseries
    data_ts{i} = timeseries(data_filtered{i}',time_cell{i}');

    % resample data
    data_resampled{i} = resample(data_ts{i},time);
end

% rename variables
engine_rpm = data_resampled{1}.Data;
lat_accel = data_resampled{2}.Data;
long_accel = data_resampled{3}.Data;
wheelspeed_FL = data_resampled{4}.Data;
wheelspeed_FR = data_resampled{5}.Data;
yaw_rate = data_resampled{6}.Data;
steer_angle = data_resampled{7}.Data;
wheelspeed_RL = data_resampled{8}.Data;
steering_voltage = data_resampled{9}.Data;

%% Save data
save('18in_lowrideheight_raw.mat','engine_rpm','lat_accel','long_accel','wheelspeed_FL',...
    'wheelspeed_FR','wheelspeed_RL','yaw_rate','steer_angle','time','steering_voltage');
