%% Load Data

clear all;clc

% import mat file from motec
load('AlexEndurance2.mat');

% input variable names of channels 
data = {Engine_RPM,G_Force_Lat,G_Force_Long,Ground_Speed_Left,Ground_Speed_Right...
    Gyro_Yaw_Velocity,Steering_Angle,Drive_Speed_Left,Susp_Pos_FL,Susp_Pos_FR,Susp_Pos_RL,Susp_Pos_RR};

%% Pre-process Data

for i = 1:numel(data)
    % sampling frequency
    Fs = 1/(data{i}.Time(2)-data{i}.Time(1))';

    % filter data
    d = fdesign.lowpass('Fp,Fst,Ap,Ast',4,5,0.001,50,Fs);
    Hd = design(d,'equiripple'); 
    data_filtered{i} = filter(Hd,data{i}.Value);
    
    % shift data by delay caused by filtering
    grp = grpdelay(Hd);
    delay = ceil(mean(grp));
    data_filtered{i}(1:delay) = [];
    
    % remove shifted elements from end of time vector
    time_cell{i} = data{i}.Time(1:end-delay);
end

% find minimum sample size
[minsize, minindex] = min(cellfun('size', data_filtered, 2));
minindex = 2;
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
save('AlexEndurance2_processed.mat','engine_rpm','lat_accel','long_accel','wheelspeed_FL',...
    'wheelspeed_FR','wheelspeed_RL','yaw_rate','steer_angle','time','steering_voltage');

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
plot(time,wheelspeed_RL)
%hold on
%plot(time,wheelspeed_RL)
xlabel('Time(s)');
ylabel('Wheelspeeds (mph)');
legend('FR','RL')
ylim([-inf inf])

linkaxes([ax1,ax2,ax3,ax4],'x')
%figure
%scatter(-lat_accel,-steer_angle)
% hold on
% 
% ackermann = 1.524*(yaw_rate+40)./(wheelspeed_FL*0.447);
% scatter(-lat_accel,-steer_angle+ackermann)
% ylim([-20 20]);

%velocity = (wheelspeed_FL+wheelspeed_FR)/2*pi/180*0.2286;


