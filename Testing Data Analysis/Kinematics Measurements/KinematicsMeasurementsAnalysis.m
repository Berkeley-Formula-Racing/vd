%% Ride Rates

clear all;close all;clc

% columns: FL FR RL RR
weights = [68	69.2	87.7	88.3
80.5	81.2	94.5	98.2
86.7	87.7	103.3	109.3
91.9	89.7	115.2	117.6
94.4	92.2	123.8	126.8
97.7	94.9	130.2	131.6		
117.6	117.2	116.2	121.7
108.9	105.7	126.7	129.3];

for i = 1:4
    weights(:,i) = weights(:,i)-weights(1,i);
end

displacement = [0	0	0	0
0.048	0.041	0.017	0.019
0.074	0.082	0.049	0.052
0.103	0.103	0.102	0.1
0.121	0.117	0.151	0.148
0.139	0.135	0.184	0.18			
0.264	0.262	0.164	0.165
0.216	0.214	0.2	0.196];

x = linspace(0,0.2,1000);
names = {'Front Left','Front Right','Rear Left','Rear Right'};
for i = 1:4
    figure
    scatter(displacement(1:end-2,i),weights(1:end-2,i))
    stiffness = displacement(1:end-2,i)\weights(1:end-2,i);
    hold on
    plot(x,stiffness*x)
    xlabel('Corner Deflection (in)')
    ylabel('Change in Corner Weight (lb)')
    title(names{i})
    legend(['Stiffness = ' num2str(round(stiffness,0)) ' (lb/in)'],'Location','Northwest')
end

%% Bump Measurements
close all;clc;

supply_voltage = 4.97;
shockpot_travel = 2.165;

% columns: vertical displacement, front mid, rear mid, bottom, top,
% shockpot V
RR = [0	0	0	0	0	0.279
0.125	0.006	0.008	0.022	-0.009	0.568
0.25	0.012	0.019	0.045	-0.017	0.898
0.375	0.017	0.028	0.067	-0.027	1.229
0.5	0.021	0.037	0.089	-0.038	1.545
0.625	0.024	0.047	0.11	-0.049	1.888
0.75	0.026	0.055	0.131	-0.06	2.236
0.875	0.03	0.063	0.151	-0.073	2.594
1	0.03	0.07	0.171	-0.086	2.984
1.125	0.031	0.078	0.189	-0.099	3.351
1.25	0.03	0.085	0.209	-0.114	3.711
1.375	0.029	0.09	0.225	-0.128	4
1.5	0.027	0.095	0.24	-0.158	4.311
1.625	0.023	0.101	0.258	-0.159	4.675];

RL = [0	0	0	0	0	0.042
0.125	-0.001	-0.002	0.012	-0.015	0.198
0.25	-0.002	0.001	0.026	-0.029	0.529
0.375	-0.007	0.001	0.038	-0.048	0.853
0.5	-0.012	0.001	0.049	-0.068	1.182
0.625	-0.017	0.001	0.059	-0.086	1.507
0.75	-0.023	-0.002	0.069	-0.106	1.885
0.875	-0.03	-0.005	0.077	-0.128	2.263
1	-0.038	-0.005	0.087	-0.149	2.697
1.125	-0.049	-0.007	0.094	-0.172	3.074
1.25	-0.056	-0.012	0.102	-0.196	3.403
1.375	-0.066	-0.015	0.11	-0.22	3.741
1.5	-0.078	-0.019	0.116	-0.247	4.07
1.625	-0.09	-0.024	0.122	-0.272	4.377];

FL = [0	0	0	0	0	0.381
0.125	-0.007	-0.007	0.007	-0.022	0.618
0.25	-0.011	-0.011	0.014	-0.042	0.847
0.375	-0.015	-0.014	0.023	-0.061	1.1
0.5	-0.021	-0.018	0.031	-0.082	1.336
0.625	-0.028	-0.021	0.038	-0.103	1.572
0.75	-0.034	-0.023	0.046	-0.122	1.799
0.875	-0.045	-0.027	0.05	-0.146	2.036
1	-0.053	-0.029	0.059	-0.168	2.304
1.125	-0.067	-0.034	0.062	-0.194	2.542
1.25	-0.081	-0.038	0.067	-0.22	2.796
1.375	-0.093	-0.04	0.073	-0.245	3.031
1.5	-0.106	-0.041	0.078	-0.269	3.249
1.625	-0.119	-0.042	0.084	-0.296	3.495
1.75	-0.142	-0.05	0.084	-0.329	3.763
1.875	-0.159	-0.052	0.087	-0.357	3.983
2	-0.174	-0.053	0.092	-0.384	4.221
2.125	-0.207	-0.063	0.121	-0.42	4.555
2.25	-0.215	-0.069	0.152	-0.445	4.804];

FR = [0	0	0.025	0	0	0.415
0.125	-0.003	0.019	0.007	-0.02	0.569
0.25	-0.004	0.014	0.017	-0.038	0.808
0.375	-0.007	0.011	0.026	-0.058	1.065
0.5	-0.012	0.007	0.035	-0.079	1.286
0.625	-0.016	0.003	0.044	-0.099	1.538
0.75	-0.025	-0.001	0.053	-0.123	1.766
0.875	-0.03	-0.005	0.06	-0.145	2.05
1	-0.04	-0.008	0.067	-0.171	2.278
1.125	-0.049	-0.013	0.075	-0.194	2.55
1.25	-0.06	-0.015	0.081	-0.219	2.807
1.375	-0.07	-0.019	0.086	-0.241	3.02
1.5	-0.083	-0.024	0.093	-0.27	3.26
1.625	-0.097	-0.027	0.1	-0.299	3.528
1.75	-0.112	-0.031	0.105	-0.328	3.767
1.875	-0.126	-0.034	0.112	-0.357	4.018
2	-0.143	-0.037	0.119	-0.388	4.252
2.125	-0.163	-0.041	0.123	-0.421	4.505
2.25	-0.18	-0.046	0.129	-0.452	4.748];

info = {FL,FR,RL,RR};

%% MR

for i = 1:4
    wheel = info{i};
    shockpot_displacement = wheel(:,6)/supply_voltage*shockpot_travel;
    shockpot_displacement = shockpot_displacement-shockpot_displacement(1);
    wheel_displacement = wheel(:,1);
    
    figure
    x = linspace(0,max(wheel_displacement));
    scatter(wheel_displacement,shockpot_displacement)
    hold on
    MR = wheel_displacement\shockpot_displacement;
    
    p = polyfit(wheel_displacement,shockpot_displacement,3);
    y = polyval(p,x);
    MR_equation = polyder(p);
    y2 = polyval(MR_equation,x);
    
    pp=csaps(wheel_displacement,shockpot_displacement,0.999);
    p_der=fnder(pp,1);
    y_prime=ppval(p_der,x);
    
    plot(wheel_displacement(2:end)-0.125/2,diff(shockpot_displacement)./diff(wheel_displacement))
    
    plot(x,MR*x)
    plot(x,y)
    plot(x,y2)
    
    plot(x,ppval(pp,x))
    plot(x,y_prime)
    xlabel('Wheel Vertical Displacement (in)')
    %ylabel('Shock Displacement (lb)')
    ylabel('MR')
    title(names{i})
    legend(['MR = ' num2str(round(MR,3)) ' (lb/in)'],'Location','Northwest')
end

%% Bump Steer
toe_dist = (5.5+5.5);
camber_dist = (6.375+4.875);
camber_dist_FR = (6.375+4.375);

for i = 1:4
    wheel = info{i};
    steer = asind((wheel(:,2)-wheel(:,3))/toe_dist);
    wheel_displacement = wheel(:,1);
    
    figure
    x = linspace(0,max(wheel_displacement));
    scatter(wheel_displacement,steer)
    %hold on
    %MR = wheel_di-0splacement\shockpot_displacement;
    xlabel('Wheel Vertical Displacement (in)')
    ylabel('Steer Angle (deg)')
    title(names{i})
end

%% Bump Camber
for i = 1:4
    wheel = info{i};
    camber = asind((wheel(:,5)-wheel(:,4))/camber_dist);
    wheel_displacement = wheel(:,1);
    
    figure
    x = linspace(0,max(wheel_displacement));
    scatter(wheel_displacement,camber)
    %hold on
    %MR = wheel_di-0splacement\shockpot_displacement;
    xlabel('Wheel Vertical Displacement (in)')
    ylabel('Camber Angle (deg)')
    title(names{i})
end








