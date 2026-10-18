%% Load Data

%clear all;clc

% import mat file from motec
load('ConstantRadiusRun1_processed.mat');

% input variable names of channels 
data = {engine_rpm,lat_accel,long_accel,shockpot_FR,shockpot_RR,shockpot_RL,...
    shockpot_FL,wheelspeed_FL,wheelspeed_FR,yaw_rate,steer_angle,time};

%% Constant Radius Calculations

%Run 1
indices = (time > 125 & time < 147)...
    | (time > 170 & time < 210)...
    | (time > 241 & time < 307);

% % % % % %Run 2
% indices = (time > 40 & time < 165)...
%         | (time > 210 & time < 255);
% % % %| (time > 295 & time < 320)... bad
% 
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
ax1 = subplot(5,1,1);
plot(time,steer_angle)
ylim([-inf inf])
xlabel('Time(s)');
ylabel('Steer Angle (deg)');

ax2 = subplot(5,1,2);
plot(time,lat_accel)
xlabel('Time(s)');
ylabel('Lateral Acceleration (g)');
ylim([-inf inf])

ax3 = subplot(5,1,3);
plot(time,yaw_rate)
xlabel('Time(s)');
ylabel('Yaw Rate (deg/s)');
ylim([-inf inf])

ax4 = subplot(5,1,4);
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
pseudo_lat_accel = -velocity.^2/.625/9.81;
plot(time,pseudo_lat_accel);
hold on
plot(time,lat_accel);

linkaxes([ax1,ax2],'x')

figure
%scatter(lat_accel,steer_angle)
%hold on
%scatter(pseudo_lat_accel,steer_angle)
xlim([0 inf]);
ylim([5 inf]);
xlabel('Lateral Acceleration (g)')
ylabel('Steer Angle (deg)')

hold on
steer_angle2 = -movmean(steer_angle,400);
pseudo_lat_accel2 = -movmean(pseudo_lat_accel,400);

steer_angle2 = steer_angle2(1:40:end);
pseudo_lat_accel2 = pseudo_lat_accel2(1:40:end);

scatter(pseudo_lat_accel2/9.81,steer_angle2)

[fitresult, gof] = createFit(pseudo_lat_accel2/9.81, steer_angle2);

%%
% close all
% 
% ipt = findchangepts(lat_accel,'MaxNumChanges',50,'MinDistance',20);
% 
% keep_indices = logical(zeros(size(lat_accel)));
% for i = 1:numel(ipt)-1
%     
%     if std(lat_accel(ipt(i):ipt(i+1))) > 0.1
%         continue
%     end
%     
%     if std(steer_angle(ipt(i):ipt(i+1))) > 2
%         continue
%     end
%     
%     if std(pseudo_lat_accel(ipt(i):ipt(i+1))) > 0.2
%         continue
%     end
%     
%     if abs(mean(lat_accel(ipt(i):ipt(i+1))) - mean(pseudo_lat_accel(ipt(i):ipt(i+1)))) < 0.2
%         keep_indices(ipt(i):ipt(i+1)) = true;
%     end
%     
%     avg_lat_accel(i) = mean(lat_accel(ipt(i):ipt(i+1)));
%     avg_steer_angle(i) = mean(steer_angle(ipt(i):ipt(i+1)));
% end
% 
% lat_accel_new = lat_accel(keep_indices);
% pseudo_lat_accel_new = pseudo_lat_accel(keep_indices);
% steer_angle_new = steer_angle(keep_indices);
% time_new = time(keep_indices);
% 
% figure
% findchangepts(lat_accel,'MaxNumChanges',50,'MinDistance',20)
% hold on
% plot(pseudo_lat_accel)
% 
% figure
% plot(time,lat_accel)
% hold on
% scatter(time_new,lat_accel_new)
% scatter(time_new,pseudo_lat_accel_new)
% 
% figure
% scatter(lat_accel_new,steer_angle_new)
% % hold on
% % scatter(pseudo_lat_accel_new,steer_angle_new)
% xlim([-inf 0]);
% ylim([-inf 0]);
% xlabel('Lateral Acceleration (g)','FontSize',15)
% ylabel('Steer Angle (deg)','FontSize',15)
% title('Filtered','FontSize',18)
% 
% figure
% scatter(avg_lat_accel,avg_steer_angle)
% xlim([-inf 0]);
% ylim([-inf 0]);
% xlabel('Lateral Acceleration (g)','FontSize',15)
% ylabel('Steer Angle (deg)','FontSize',15)
% title('Filtered','FontSize',18)

%%
% close all
% 
% figure
% ax1 = subplot(3,1,1);
% plot(time,steer_angle)
% ylim([-inf inf])
% xlabel('Time(s)');
% ylabel('Steer Angle (deg)');
% 
% ax2 = subplot(3,1,2);
% plot(time,lat_accel)
% hold on
% scatter(time_new,lat_accel_new);
% xlabel('Time(s)');
% ylabel('Lateral Acceleration (g)');
% ylim([-inf inf])
% 
% ax3 = subplot(3,1,3);
% plot(time,wheelspeed_FL)
% hold on
% plot(time,wheelspeed_FR)
% xlabel('Time(s)');
% ylabel('Wheelspeeds (deg/s)');
% ylim([-inf inf])
% 
% linkaxes([ax1,ax2,ax3],'x')
% 
% figure
% scatter(lat_accel_new,steer_angle_new)

