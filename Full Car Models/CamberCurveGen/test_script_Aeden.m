%taken from carConfig
load('Fx_combined_parameters_run38_30.mat'); % F_x combined magic formula parameters
Fx_parameters = cell2mat(Xbestcell);
load('Lapsim_Fy_combined_parameters_1965run15.mat'); % F_y combined magic formula parameters
Fy_parameters = cell2mat(Xbestcell);
friction_scaling_factor = 1.05*0.55;
p_i = 12; % pressure

tire = Tire2(0,p_i,Fx_parameters,Fy_parameters,friction_scaling_factor);

%% run
clc
F_z_L = 100;
F_z_R = 100;

camber = -1;
gamma_L = -camber;
gamma_R = camber;

clc
[F_y_tot, F_y_L, F_y_R, M_x_L, M_x_R, alpha_val] = ...
    singleAxleCamberEvaluation(F_z_L, F_z_R, gamma_L, gamma_R, tire);
F_y_tot
F_y_L
F_y_R

%% test a
F_z_R_vector = 0:300;
F_z_L = 100;
F_y_tot_vector = zeros(size(F_z_R_vector));
F_y_L_vector = zeros(size(F_z_R_vector));
F_y_R_vector = zeros(size(F_z_R_vector));
alpha_vector = zeros(size(F_z_R_vector));
for i = 1:numel(F_z_R_vector)
    [F_y_tot, F_y_L, F_y_R, M_x_L, M_x_R, alpha_val] = ...
    singleAxleCamberEvaluation(F_z_L, F_z_R_vector(i), gamma_L, gamma_R, tire);
    F_y_tot_vector(i) = F_y_tot;
    F_y_L_vector(i) = F_y_L;
    F_y_R_vector(i) = F_y_R;
    alpha_vector(i) = alpha_val;
end

plot(F_z_R_vector, F_y_tot_vector);
hold on;
%figure;
%plot(F_z_R_vector, alpha_vector);



%% test

gamma_vector = -1:0.1:1;
alpha_vector = -20:0.5:20;

F_y_matrix = zeros(numel(gamma_vector), numel(alpha_vector));

for i = 1:numel(gamma_vector)
    tire.gamma = gamma_vector(i);
    for j = i:numel(alpha_vector)
        F_y_matrix(i,j) = F_y(tire,alpha_vector(j),0,100);
    end
end

hold off;
for i = 1:numel(gamma_vector)
    plot(alpha_vector,F_y_matrix(i,:),'DisplayName', ['gamma = ' char(string(gamma_vector(i)))]);
    hold on;
end
legend();
xline(0);
yline(0);
%% test displacement code

%Vehicle Width in inches
C.front_width = 47;
C.rear_width = 47;

%Spring Roll Stiffness
C.front_spring_stiff = 2860; %numbers from LLTD Doc
C.rear_spring_stiff = 2911;

%ARB Roll Stiffness (Currently on Short)
C.front_ARB_stiff = 0;
C.rear_ARB_stiff = 699;

%Car Weight Distribution
C.weight_dist = 0.54; % percentage of weight in rear

%Total Car Weight
C.mass = 395; % not including driver (lb)

roll_angle = 1;

[normal_load_FL, dist_FL, normal_load_FR, dist_FR, normal_load_RL,dist_RL, normal_load_RR, dist_RR] = calcWheelForcesAndDisplacements(roll_angle, C);
clc
format compact;
normal_load_FL
dist_FL
normal_load_FR
dist_FR
normal_load_RL
dist_RL
normal_load_RR
dist_RR

