%% Inputs
clear all
clc

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'Full Car Models','config'));
sharedStudyVehicle = vdVehicleConfig();
weight_dist = sharedStudyVehicle.studies.Bicycle_Model_YawResponseTimePlotter.weight_dist; % rearwards
m = sharedStudyVehicle.studies.Bicycle_Model_YawResponseTimePlotter.m; % mass, kg
L = sharedStudyVehicle.studies.Bicycle_Model_YawResponseTimePlotter.L; % wheelbase, m

u = 20; % forward velocity, m/s

I_zz = sharedStudyVehicle.studies.Bicycle_Model_YawResponseTimePlotter.I_zz; % yaw inertia, kg-m^2

%C_f = 136*4.448; % front cornering stiffness, N/deg @ 150 lb
%C_r = 136*4.448; % rear cornering stiffness, N/deg @ 150 lb
g = 9.806; % gravitational constant, m/sec^2

N = sharedStudyVehicle.studies.Bicycle_Model_YawResponseTimePlotter.N; % steering ratio

%% Outputs
a = L*weight_dist; % front axle to cg, m
b = L*(1-weight_dist); % rear axle to cg, m

D_r_range = 1.5:0.5:3;
D_f_range = 1.5:0.5:6;

D_rsweep = D_r_range;
for j = 1:numel(D_rsweep)
    D_r = D_rsweep(j);
    D_fsweep = D_r:0.1:D_f_range(end);
    for i = 1:numel(D_fsweep)
        D_f = D_fsweep(i);
        r_delta = tf([(57.3*g*a*b*u)/(D_f*L),((57.3*g)^2*a*b)/(D_f*D_r*L)],...
            [(I_zz*u/m),(57.3*g*(a^2*b*m*D_r+a*D_f*(b^2*m+I_zz)+b*D_r*I_zz))/(D_f*D_r*m*L),...
            (57.3*g*a*b*(57.3*g*L+u^2*(D_f-D_r)))/(D_f*D_r*u*L)]);
        
        % divide by SSG to normalize to unit gain
        SSG_r_delta = dcgain(r_delta);
        r_delta = r_delta/SSG_r_delta;
        
        r_stepinfo = stepinfo(r_delta);
        r_responsetime(i) = r_stepinfo.RiseTime;
    end
    D_total = D_fsweep+D_r;
    plot(D_total,r_responsetime,'k','LineWidth',3);
    hold on
    
    clear r_responsetime
end

%%

D_fsweep = D_f_range;
for j = 1:numel(D_fsweep)
    D_f = D_fsweep(j);
    D_rsweep = D_r_range(1):0.1:min(D_f,D_r_range(end));
    for i = 1:numel(D_rsweep)
        D_r = D_rsweep(i);
        r_delta = tf([(57.3*g*a*b*u)/(D_f*L),((57.3*g)^2*a*b)/(D_f*D_r*L)],...
            [(I_zz*u/m),(57.3*g*(a^2*b*m*D_r+a*D_f*(b^2*m+I_zz)+b*D_r*I_zz))/(D_f*D_r*m*L),...
            (57.3*g*a*b*(57.3*g*L+u^2*(D_f-D_r)))/(D_f*D_r*u*L)]);
        
        % divide by SSG to normalize to unit gain
        SSG_r_delta = dcgain(r_delta);
        r_delta = r_delta/SSG_r_delta;
        
        r_stepinfo = stepinfo(r_delta);
        r_responsetime(i) = r_stepinfo.RiseTime;
        %end
    end
    D_total = D_rsweep+D_f;
    plot(D_total,r_responsetime,'k','LineWidth',3);
    hold on
    
    startpoint(j) = r_responsetime(1);
    endpoint(j) = r_responsetime(end);
    clear r_responsetime
end

for i = 1:(numel(D_f_range)-numel(D_r_range)+1)
    plot([D_f_range(i)+D_r_range(1) D_r_range(end)+D_f_range(numel(D_r_range)+i-1)]...
        ,[startpoint(i) endpoint(numel(D_r_range)+i-1)],'--','Color','b','LineWidth',3);
    xlim([D_f_range(1)+D_r_range(1)-0.75 D_f_range(end)+D_r_range(end)+0.85]);
    ylim([startpoint(end)-0.04 max(endpoint)+0.02]);
end

%%
for i = 1:numel(D_f_range)
    x = [D_r_range(1)+D_f_range(i) D_r_range(1)+D_f_range(i)];
    y = [startpoint(i)-0.01 startpoint(i)];
    
    xnorm = x_to_norm_v2(x(1),x(2)); %[0.01 0.01];%(x-x_min)/(x_max-x_min)
    ynorm = y_to_norm_v2(y(1),y(2)); %[0.01 0.01]; %(y-y_min)/(y_max-y_min)
    annotation('textarrow',xnorm,ynorm,'String',num2str(D_f_range(i)));
  
    hold on
end

for i = 1:numel(D_r_range)
    x = [D_f_range(i)+D_r_range(i)-0.25 D_f_range(i)+D_r_range(i)];
    y = [endpoint(i) endpoint(i)];
    
    xnorm = x_to_norm_v2(x(1),x(2)); %[0.01 0.01];%(x-x_min)/(x_max-x_min)
    ynorm = y_to_norm_v2(y(1),y(2)); %[0.01 0.01]; %(y-y_min)/(y_max-y_min)
    annotation('textarrow',xnorm,ynorm,'String',num2str(D_r_range(i)));
  
    hold on
end

counter = 0;
for i = numel(D_r_range):numel(endpoint)
    x = [D_f_range(i)+D_r_range(end)+0.25 D_f_range(i)+D_r_range(end)];
    y = [endpoint(i) endpoint(i)];
    
    xnorm = x_to_norm_v2(x(1),x(2)); %[0.01 0.01];%(x-x_min)/(x_max-x_min)
    ynorm = y_to_norm_v2(y(1),y(2)); %[0.01 0.01]; %(y-y_min)/(y_max-y_min)
    annotation('textarrow',xnorm,ynorm,'String',['K = ' num2str(counter)]);
    counter = counter+1;
    hold on
end

annotation('textbox',[.15 .7 .1 .1],'String','D_r: Rear Cornering Compliance')
annotation('textbox',[.25 .1 .1 .1],'String','D_f: Front Cornering Compliance')
annotation('textbox',[.65 .7 .1 .1],'String','Understeer Gradient')

xlabel('Total Cornering Compliance (deg/g)','FontSize',15);
ylabel('Yaw Rate Response Time (s)','FontSize',15);
title(['Yaw Rate Response Time, ' num2str(u) ' m/s'],'FontSize',18);

%%
% figure
% step(ay_delta)
% xlabel('Time (sec)','FontSize',15)
% ylabel('Lateral Acceleration (g)','FontSize',15)
% title('Lateral Acceleration Step Response')
