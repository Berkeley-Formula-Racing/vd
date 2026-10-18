clear all;close all;clc

%% Front Camber

% columns: force (lb), dial indicator, camber gauge (deg)
front_camber = {[50	0.069	0.3
100	0.111	0.4
150	0.142	0.6
200	0.19	0.7
250	0.21	0.9
300	0.24	1],...
[50	0.085	0.3
100	0.11	0.4
150	0.149	0.5
200	0.181	0.7
250	0.214	0.8
300	0.242	0.9],...
[50	0.08	0.4
100	0.124	0.5
150	0.14	0.6
200	0.171	0.8
250	0.2	0.9
300	0.22	1]};

moment_arm = 9; % in
tire_radius = 8; % in

figure
for i = 1:3
    force = front_camber{i}(:,1); % lb
    moment = force*moment_arm; % in-lb
    tire_force = moment/tire_radius; % lb
    
    camber = front_camber{i}(:,3); % deg
    scatter(camber,tire_force)    
    x = linspace(0.2,1.1,1000);
    
    p = polyfit(camber,tire_force,1);
    y = polyval(p,x);
    
    stiffness(i) = p(1);
    
    hold on
    plot(x,y)
    xlabel('Camber Deflection (deg)')
    ylabel('Force (lb)')
    title('Front Camber Stiffness')
    %legend(['Stiffness = ' num2str(round(stiffness(i),0)) ' (lb/deg)'],'Location','Northwest')
    xlim([min(x) max(x)])
end

figure
scatter(camber,tire_force)
hold on
plot(x,y)
xlabel('Camber Deflection (deg)')
ylabel('Force (lb)')
title('Front Camber Stiffness')
%legend(['Stiffness = ' num2str(round(stiffness(i),0)) ' (lb/deg)'],'Location','Northwest')
xlim([min(x) max(x)])

average_front_camber_stiffness = mean(stiffness) %lb/deg
front_lateral_force = 471; % front, 2G, 0.5 LLTD
front_camber_compliance = front_lateral_force/(2*mean(stiffness)) % deg/g

%% Rear Camber

% columns: force (lb), dial indicator, camber gauge (deg)
rear_camber = {[50	0.3
100	0.5
150	0.6
200	0.7
250	0.8
300	1],...
[50	0.4
100	0.5
150	0.6
200	0.8
250	0.9
300	1],...
[50	0.3
100	0.5
150	0.6
200	0.7
250	0.8
300	1]};

moment_arm = 9; % in
tire_radius = 8; % in

figure
for i = 1:3
    force = rear_camber{i}(:,1); % lb
    moment = force*moment_arm; % in-lb
    tire_force = moment/tire_radius; % lb
    
    camber = rear_camber{i}(:,2); % deg
    scatter(camber,tire_force)    
    x = linspace(0.2,1.1,1000);
    
    p = polyfit(camber,tire_force,1);
    y = polyval(p,x);
    
    stiffness(i) = p(1);
    
    hold on
    plot(x,y)
    xlabel('Camber Deflection (deg)')
    ylabel('Force (lb)')
    title('Rear Camber Stiffness')
    %legend(['Stiffness = ' num2str(round(stiffness(i),0)) ' (lb/deg)'],'Location','Northwest')
    xlim([min(x) max(x)])
end

figure
scatter(camber,tire_force)
hold on
plot(x,y)
xlabel('Camber Deflection (deg)')
ylabel('Force (lb)')
title('Rear Camber Stiffness')
%legend(['Stiffness = ' num2str(round(stiffness(i),0)) ' (lb/deg)'],'Location','Northwest')
xlim([min(x) max(x)])

average_rear_camber_stiffness = mean(stiffness) %lb/deg

rear_lateral_force = 516; % rear, 2G, 0.5 LLTD
rear_camber_compliance = rear_lateral_force/(2*mean(stiffness)) % deg/g

%% Front Toe Compliance

% columns: force (lb), dial indicator, calculated toe angle
front_toe = {[13	0.04	0.4365434961
20	0.078	0.8512829012
30	0.12	1.309731879
40	0.15	1.637245078
50	0.17	1.855616281
60	0.19	2.074014453],...
[10	0.025	0.2728380764
20	0.045	0.4911126953
30	0.08	0.8731123366
40	0.14	1.528068595
50	0.174	1.899293656
60	0.205	2.237832694],...
[10	0.015	0.1637024499
20	0.055	0.6002524797
30	0.085	0.9276864863
40	0.115	1.255150796
50	0.165	1.801021077
60	0.19	2.074014453]};

moment_arm = 9; % in

figure
for i = 1:3
    force = front_toe{i}(:,1); % lb
    moment = force*moment_arm/12; % ft-lb
    
    toe = front_toe{i}(:,3); % deg
    scatter(toe,moment)    
    x = linspace(0.1,2.2,1000);
    
    p = polyfit(toe,moment,1);
    y = polyval(p,x);
    
    stiffness(i) = p(1);
    
    hold on
    plot(x,y)
    xlabel('Toe Deflection (deg)')
    ylabel('Moment (ft-lb)')
    title('Front Toe Stiffness')
    %legend(['Stiffness = ' num2str(round(stiffness(i),0)) ' (lb/deg)'],'Location','Northwest')
    xlim([min(x) max(x)])
end

figure
scatter(toe,moment)
hold on
plot(x,y)
xlabel('Toe Deflection (deg)')
ylabel('Moment (ft-lb)')
title('Front Toe Stiffness')
%legend(['Stiffness = ' num2str(round(stiffness(i),0)) ' (lb/deg)'],'Location','Northwest')
xlim([min(x) max(x)])

average_front_toe_stiffness = mean(stiffness) %ft-lb/deg

front_aligning_torque = 471/12; % (ft-lb), front, 2G, 0.5 LLTD
front_toe_compliance = front_aligning_torque/(2*mean(stiffness)) % deg/g

%% Rear Toe Compliance

% columns: force (lb), dial indicator, calculated toe angle
rear_toe = {[10	0.011	0.1200483878
20	0.021	0.2291837292
30	0.033	0.3601472714
40	0.041	0.4474573027
50	0.054	0.5893384098
60	0.062	0.6766516012],...
[10	0.009	0.09822138442
20	0.018	0.1964430575
30	0.027	0.2946653079
40	0.038	0.4147159301
50	0.045	0.4911126953
60	0.056	0.6111665714],...
[10	0.008	0.08730788828
20	0.017	0.185529515
30	0.028	0.3055789394
40	0.037	0.4038021698
50	0.049	0.534768373
60	0.056	0.6111665714]};

moment_arm = 9; % in

figure
for i = 1:3
    force = rear_toe{i}(:,1); % lb
    moment = force*moment_arm/12; % ft-lb
    
    toe = rear_toe{i}(:,3); % deg
    scatter(toe,moment)    
    x = linspace(0.1,0.7,1000);
    
    p = polyfit(toe,moment,1);
    y = polyval(p,x);
    
    stiffness(i) = p(1);
    
    hold on
    plot(x,y)
    xlabel('Toe Deflection (deg)')
    ylabel('Moment (ft-lb)')
    title('Rear Toe Stiffness')
    %legend(['Stiffness = ' num2str(round(stiffness(i),0)) ' (lb/deg)'],'Location','Northwest')
    xlim([min(x) max(x)])
end

figure
scatter(toe,moment)
hold on
plot(x,y)
xlabel('Toe Deflection (deg)')
ylabel('Moment (ft-lb)')
title('Rear Toe Stiffness')
%legend(['Stiffness = ' num2str(round(stiffness(i),0)) ' (lb/deg)'],'Location','Northwest')
xlim([min(x) max(x)])

average_rear_toe_stiffness = mean(stiffness) %ft-lb/deg

rear_aligning_torque = 516/12; % (ft-lb), front, 2G, 0.5 LLTD
rear_toe_compliance = rear_aligning_torque/(2*mean(stiffness)) % deg/g





