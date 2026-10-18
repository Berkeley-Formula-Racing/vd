%% Load Data

clear all;clc

% import mat file from motec
%load('18in_lowrideheight_steering_1Hz.mat');
load('steering_calibration_processed.mat');
load('SteeringEffortCalibration3_processed.mat');

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
% indices = (time > 149 & time < 157)|...
%  (time > 166 & time < 174)|...
% (time > 184 & time < 191)|...
% (time > 270 & time < 278)|...
% (time > 288 & time < 294)|...
% (time > 320 & time < 326)|...
% (time > 372 & time < 380)|...
% (time > 390 & time < 397)|...
% (time > 195 & time < 215)|...
% (time > 235 & time < 265);

% %18 in slalom (low ride height)
%indices = (time > 166 & time < 174);

% %18 in skidpad (low ride height)
% indices = (time > 195 & time < 215)|...
%      (time > 235 & time < 265);

%18 in slalom (low ride height)
% indices = (time > 120);
% 
% 

% steeringcalibration
indices = (time>216 & time<222);
%indices = (time>445 & time<455); %|(time>720 & time<740);
%indices = (time>720 & time<740);

indices = (time>659 & time<663);
indices = (time>695 & time<697);
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
close all

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
steer_torque = 410/12*(steering_voltage-mean(steering_voltage));
%steer_torque = steer_torque - 2;
figure
%plot(lat_accel,steer_torque)
%plot(lat_accel(4:end),steer_torque(1:end-3))
scatter(lat_accel(4:end),steer_torque(1:end-3))
xlabel('Lateral Acceleration (g)')
ylabel('Steering Effort (ft-lb)')

linkaxes([ax1,ax2,ax3,ax4],'x')

%%
plot(steer_angle+5,1.2*(steer_torque+0.5))
xlim([-20 20])
xlabel('Steer Angle (deg)')
ylabel('Steering Effort (ft-lb)')
title('Static Steering Effort')

%%
Fs = 1/(time(2)-time(1));

windowsize = round(length(steer_angle)/5);
[sstxy,f] = tfestimate(steer_angle,lat_accel,windowsize,[],[],Fs);
figure
subplot(2,1,1)
plot(f,abs(sstxy/sstxy(1)))
xlim([0 2]);
subplot(2,1,2)
plot(f,180/pi*angle(sstxy/sstxy(1)))
hold on
plot(f,180/pi*angle(sstxy2/sstxy2(1)))
xlim([0 2]);
legend('Lateral Acceleration','Yaw Rate')

%%

[xData, yData] = prepareCurveData(lat_accel(4:end),1.2*(steer_torque(1:end-3)-2));

% Set up fittype and options.
ft = fittype( 'rat12' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.890903252535798 0.959291425205444 0.547215529963803 0.138624442828679];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

indices = false(size(xData));
for i = 1:numel(xData)
    random = normrnd(1,4);
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
%h = plot(fitresult, xData, yData );


scatter(xData,yData,'.','b');
% Label axes
xlabel('Lateral Acceleration (g)','FontSize',10)
ylabel('Steering Effort (ft-lb)','FontSize',10)
legend('hide')

%%
figure
scatter(lat_accel,steer_angle)
hold on

ackermann = 1.524*(yaw_rate+40)./(wheelspeed_FL*0.447);
scatter(-lat_accel,-steer_angle+ackermann)
ylim([-20 20]);
