%% Load Data

clear all;clc

% import mat file from motec
load('Step4_processed.mat');

% input variable names of channels 
data = {engine_rpm,lat_accel,long_accel,shockpot_FR,shockpot_RR,shockpot_RL,...
    shockpot_FL,wheelspeed_FL,wheelspeed_FR,yaw_rate,steer_angle,time};

%%
% %step6 - not great, doesn't start from 0 yaw rate, short time
% indices = (time > 41.8 & time < 42.7);

% %step5
% indices = (time > 84 & time < 87);

% | (time > 107 & time < 111) ...
%     | (time > 128 & time < 132) | (time > 192 & time < 194) | ...
%     (time > 213  & time < 216);

%step4
indices = (time > 25 & time < 27);

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
xlabel('Time');
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
ylim([-200 200])

ax4 = subplot(4,1,4);
plot(time,wheelspeed_FL)
hold on
plot(time,wheelspeed_FR)
xlabel('Time(s)');
ylabel('Wheelspeeds (deg/s)');
ylim([-inf inf])

linkaxes([ax1,ax2,ax3,ax4],'x')


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

% %%
% Fs = 1/(time(2)-time(1));
% L = numel(lat_accel);
% f = Fs*(1:(L/2))/L;
% fourier_transform = fft(lat_accel);
% 
% plot(f,abs(fourier_transform(1:L/2)));

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

