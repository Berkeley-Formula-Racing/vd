%% Load Data

clear all;clc

% import mat file from motec
load('AlexEndurance_processed.mat');

% input variable names of channels 
data = {engine_rpm,lat_accel,long_accel,wheelspeed_FL,wheelspeed_FR,yaw_rate,steer_angle,time};
data{end} = data{end}';
%%

% sinusoidal3
% indices = (time > 27 & time < 38);
% indices = (time > 44 & time < 51);
% indices = (time > 58 & time < 64);
% % 
% indices = (time > 27 & time < 38)...
%     |(time > 44 & time < 51)...
%     |(time > 58 & time < 64);

% % sinusoidal4
% indices = (time > 180 & time < 207);
% indices = (time > 245 & time < 255);
% %indices = (time > 303 & time < 312);
% 
% indices = (time > 137 & time < 147)...
%     |(time > 180 & time < 207)...
% |(time > 245 & time < 255)...
% |(time > 303 & time < 312);

% indices = (time > 137);
% 
% for i = 1:numel(data)
%     variable = data{i};
%     data{i} = variable(indices);
% end


% rename variables
engine_rpm = data{1};
lat_accel = data{2};
long_accel = data{3};
wheelspeed_FL = data{4};
wheelspeed_FR = data{5};
yaw_rate = data{6};
steer_angle = data{7};
time = data{8};

% %% Save data
% save('18in_lowrideheight_slalom_only.mat','engine_rpm','lat_accel','long_accel','wheelspeed_FL',...
%     'wheelspeed_FR','wheelspeed_RL','yaw_rate','steer_angle','time');

%% Plotting
%close all

ax1 = subplot(4,1,1);
plot(time,steer_angle)
hold on
ylim([-inf inf])
xlabel('Time(s)');
ylabel('Steer Angle (deg)');

ax2 = subplot(4,1,2);
plot(time,lat_accel)
hold on
xlabel('Time(s)');
ylabel('Lateral Acceleration (g)');
ylim([-inf inf])

ax3 = subplot(4,1,3);
plot(time,yaw_rate)
hold on
xlabel('Time(s)');
ylabel('Yaw Rate (deg/s)');
ylim([-inf inf])

% ax4 = subplot(4,1,4);
% plot(wheelspeed_FL)
% hold on
% plot(wheelspeed_FR)
% xlabel('Time(s)');
% ylabel('Wheelspeeds (mph)');
% ylim([-inf inf])
% legend('FL','FR')

ax4 = subplot(4,1,4);
plot(time,(wheelspeed_FL+wheelspeed_FR)/2)
hold on
xlabel('Time(s)');
ylabel('Velocity (mph)');
ylim([-inf inf])


linkaxes([ax1,ax2,ax3,ax4],'x')


%% Transfer Function Estimation

data = iddata(lat_accel,steer_angle,time(2)-time(1));
ay_sys = tfest(data,2,2);

data = iddata(yaw_rate,steer_angle,time(2)-time(1));
r_sys = tfest(data,2,1);

% Bode plot
bodeopt = bodeoptions;
bodeopt.MagScale = 'linear';
bodeopt.MagUnits = 'abs';
bodeopt.FreqScale = 'linear';
bodeopt.FreqUnits = 'Hz';
bodeopt.Xlim = [0 4];

figure
bode(r_sys,ay_sys,bodeopt)
%% PSD
Fs = 1/(time(2)-time(1));

[pxx,f] = pwelch(steer_angle,[],[],[],Fs);

plot(f,10*log10(pxx))
xlim([0 5])
xlabel('Frequency (Hz)')
ylabel('PSD (dB/Hz)')

%% Coherence
Fs = 1/(time(2)-time(1));

windowsize = round(length(steer_angle)/40);
[cxy,fc] = mscohere(steer_angle,lat_accel,windowsize,[],[],Fs);
[cxy2,fc2] = mscohere(steer_angle,yaw_rate,windowsize,[],[],Fs);

plot(fc,cxy)
hold on
plot(fc2,cxy2)
xlim([0 5])
legend('Lat Accel','Yaw Rate')

%% Frequency Response
Fs = 1/(time(2)-time(1));

windowsize = round(length(steer_angle)/40);
[sstxy,f] = tfestimate(steer_angle,lat_accel,windowsize,[],[],Fs);
[sstxy2,f] = tfestimate(steer_angle,-yaw_rate,windowsize,[],[],Fs);
figure
subplot(2,1,1)
plot(f,abs(sstxy/sstxy(1)))
hold on
plot(f,abs(sstxy2/sstxy2(1)))
xlim([0 2]);
subplot(2,1,2)
plot(f,180/pi*angle(sstxy/sstxy(1)))
hold on
plot(f,180/pi*angle(sstxy2/sstxy2(1)))
xlim([0 2]);
legend('Lateral Acceleration','Yaw Rate')
