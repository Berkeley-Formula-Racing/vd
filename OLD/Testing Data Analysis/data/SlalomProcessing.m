%% Load Data

clear all;clc

% import mat file from motec
load('18in_lowrideheight_steering_processed.mat');

% input variable names of channels 
data = {engine_rpm,lat_accel,long_accel,wheelspeed_FL,wheelspeed_FR,yaw_rate,steer_angle,steering_voltage,time};

%%

mean_steering_voltage = mean(steering_voltage);

%16 in slalom
 %indices = (time > 91 & time < 98); %| %...
%  (time > 161 & time < 167)|...
% (time > 176 & time < 182)|...
% (time > 194 & time < 200)|...
% (time > 208 & time < 216)|...
% (time > 226 & time < 232);

% %18 in slalom (low ride height)
% indices = (time > 149 & time < 157); %|...
%  (time > 166 & time < 174)|...
% (time > 184 & time < 191)|...
% (time > 270 & time < 278)|...
% (time > 288 & time < 294)|...
% (time > 320 & time < 326)|...
% (time > 372 & time < 380)|...
% (time > 390 & time < 397);

% %18 in slalom (low ride height)
%indices = (time > 166 & time < 174);

% %18 in skidpad (low ride height)
%indices = (time > 195 & time < 215); %...
     %(time > 235 & time < 265);

%18 in slalom (low ride height)
%indices = (time > 120);
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
steering_voltage = data{8};
time = data{9};


% %% Save data
% save('18in_lowrideheight_slalom_only.mat','engine_rpm','lat_accel','long_accel','wheelspeed_FL',...
%     'wheelspeed_FR','wheelspeed_RL','yaw_rate','steer_angle','time');

%% Plotting
%close all

ax1 = subplot(4,1,1);
plot(steer_angle)
hold on
ylim([-inf inf])
xlabel('Time(s)');
ylabel('Steer Angle (deg)');

ax2 = subplot(4,1,2);
plot(lat_accel)
hold on
xlabel('Time(s)');
ylabel('Lateral Acceleration (g)');
ylim([-inf inf])

ax3 = subplot(4,1,3);
plot(yaw_rate)
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

% ax4 = subplot(4,1,4);
% plot((wheelspeed_FL+wheelspeed_FR)/2)
% hold on
% xlabel('Time(s)');
% ylabel('Velocity (mph)');
% ylim([-inf inf])

ax4 = subplot(4,1,4);
plot(steering_voltage)
hold on
xlabel('Time(s)');
ylabel('Steering Torque Voltage');
ylim([-inf inf])

% ft-lb
steer_torque = 410/12*(steering_voltage-mean_steering_voltage);
steer_torque = steer_torque - 3;
figure
scatter(lat_accel,steer_torque)
xlabel('Lateral Acceleration (g)')
ylabel('Steering Effort (ft-lb)')

linkaxes([ax1,ax2,ax3,ax4],'x')

%%

[xData, yData] = prepareCurveData( lat_accel, steer_torque );

% Set up fittype and options.
ft = fittype( 'rat12' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.890903252535798 0.959291425205444 0.547215529963803 0.138624442828679];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

indices = false(size(xData));
for i = 1:numel(xData)
    random = normrnd(1,2);
    tony = abs(fitresult(xData)-yData);
    if tony(i) < random
        indices(i) = true;
    end
end
xData = xData(indices);
yData = yData(indices);

%xData = xData(1:8:end);
%yData = yData(1:8:end);

% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( fitresult, xData, yData );
% Label axes
xlabel('Lateral Acceleration (g)')
ylabel('Steering Effort (ft-lb)')
grid on
%%
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
fourier_transform = fft(lat_accel);

%[~,x] = max(abs(fourier_transform(1:L/2)));
%frequency = f(x)
figure(1)
findpeaks(abs(fourier_transform(1:L/2)))
%[pks,locs] = findpeaks(abs(fourier_transform(1:L/2)))

figure(2)
plot(f,abs(fourier_transform(1:L/2)));
hold on
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Lateral Acceleration')
legend('16 in tires','18 in tires');

%%
Fs = 1/(time(2)-time(1));

[pxx,f] = pwelch(yaw_rate);

plot(f,10*log10(pxx))

xlabel('Frequency (Hz)')
ylabel('PSD (dB/Hz)')

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

