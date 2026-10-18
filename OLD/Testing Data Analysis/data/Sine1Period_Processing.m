%% Load Data

clear all;clc

% import mat file from motec
load('FrequencyResponse_processed.mat');

% input variable names of channels 
data = {engine_rpm,lat_accel,long_accel,shockpot_FR,shockpot_RR,shockpot_RL,...
    shockpot_FL,wheelspeed_FL,wheelspeed_FR,yaw_rate,steer_angle,time};

%%

%sine1period
% indices = (time > 83 & time < 114) | (time > 120 & time < 148); %(time > 142 & time < 146);
%indices = (time > 89 & time < 92);
% indices = (time > 97 & time < 99);
% indices = (time > 125 & time < 128);
%indices = (time > 128 & time < 131);

% % frequencyresponse
indices = (time > 62 & time < 86) | (time > 22 & time < 53) | (time > 142 & time < 146);

% % amplituderesponse
% indices = (time > 28 & time < 65) | (time > 75 & time < 117);
% indices = (time > 32 & time < 34);
% indices = (time > 35 & time < 38);
% indices = (time > 39 & time < 41.5); %good
% indices = (time > 80 & time < 82);
% indices = (time > 83 & time < 85.5); %good
% indices = (time > 91 & time < 94);
% indices = (time > 95 & time < 98); %good

for i = 1:numel(data)
    variable = data{i};
    data{i} = variable(indices);
end

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
ylabel('Wheelspeeds (mph)');
ylim([-inf inf])
legend('FL','FR')

linkaxes([ax1,ax2,ax3,ax4],'x')

figure
scatter(lat_accel,steer_angle)
hold on

ackermann = 1.524*(yaw_rate+40)./(wheelspeed_FL*0.447);
scatter(-lat_accel,-steer_angle+ackermann)
ylim([-20 20]);

%% Calculations

lat_accel_peak = (max(lat_accel)+abs(min(lat_accel)))/2
steer_angle_peak = (max(steer_angle)+abs(min(steer_angle)))/2
yaw_rate_peak = (max(yaw_rate)+abs(min(yaw_rate)))/2

lat_accel_gain = lat_accel_peak*9.81/steer_angle_peak
yaw_rate_gain = yaw_rate_peak*9.81/steer_angle_peak

D = finddelay(lat_accel,steer_angle)
D = finddelay(yaw_rate,steer_angle)


%% Transfer Function Estimation

data = iddata(lat_accel,steer_angle,time(2)-time(1));
ay_sys = tfest(data,2,2);

data = iddata(yaw_rate,steer_angle,time(2)-time(1));
r_sys = tfest(data,2,1);

%%
nfft = pow2(round(log2(length(steer_angle))));
[sstxy,f] = tfestimate(steer_angle,lat_accel);
plot(f,abs(sstxy))

%%

[sstxy,f] = tfestimate(steer_angle,yaw_rate);
plot(f,abs(sstxy))

%%
% w = linspace(0,2*pi*5);
% data = iddata(lat_accel,steer_angle,time(2)-time(1));
% g = spa(data,[],w); 
% 
% % Bode plot
% bodeopt = bodeoptions;
% bodeopt.MagScale = 'linear';
% bodeopt.MagUnits = 'abs';
% bodeopt.FreqScale = 'linear';
% bodeopt.FreqUnits = 'Hz';
% bodeopt.Xlim = [0 5];
% 
% figure
% bode(g,bodeopt)

%%

% Bode plot
bodeopt = bodeoptions;
bodeopt.MagScale = 'linear';
bodeopt.MagUnits = 'abs';
bodeopt.FreqScale = 'linear';
bodeopt.FreqUnits = 'Hz';
bodeopt.Xlim = [0 4];

figure
bode(ay_sys,bodeopt)

%%
Fs = 1/(time(2)-time(1));
L = numel(lat_accel);
f = Fs*(1:(L/2))/L;
fourier_transform = fft(steer_angle);

plot(f,abs(fourier_transform(1:L/2)));

%%
t = linspace(0, 5); % Time Vector
u = sin(10*t);        % Forcing Function
y = lsim(r_sys, u, t);    % Calculate System Response
figure
plot(t,u)
hold on
plot(t,y/5)

%%
u = steer_angle;        % Forcing Function
y = lsim(ay_sys, u, time);    % Calculate System Response
figure
plot(time,lat_accel)
hold on
plot(time,y)

%%
figure
plot(time,steer_angle)
hold on
plot(time,lat_accel*5)

