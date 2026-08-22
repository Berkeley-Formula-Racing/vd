%% setup

clc
close all
clear all

format compact;

%% taken from carConfig
load('Fx_combined_parameters_run38_30.mat'); % F_x combined magic formula parameters
Fx_parameters = cell2mat(Xbestcell);
load('Fy_combined_parameters_16x10x7.5_LC0.mat'); % F_y combined magic formula parameters
Fy_parameters = cell2mat(Xbestcell);
friction_scaling_factor = 1.05*0.55;
p_i = 12; % tire pressure

%initialize tire object, used for camber evaluation
tire = BDITire2(0,p_i,Fx_parameters,Fy_parameters,friction_scaling_factor);

%Vehicle Width in inches
C.front_width = 47;
C.rear_width = 47;

%Spring Roll Stiffness
C.front_spring_roll_stiffness = 3000;2860; %numbers from LLTD Doc
C.rear_spring_roll_stiffness = 0;2911;

%ARB Roll Stiffness (Currently on Short)
C.front_ARB_roll_stiffness = 0;
C.rear_ARB_roll_stiffness = 0;%699;986;1493;

%Car Weight Distribution
C.weight_dist = 0.54; % percentage of weight in rear

%Total Car Weight
C.mass = (395 + 150); % lbs


% [normal_load_vector, wheel_displacement_vector] = calcWheelForcesAndDisplacements(1, C);
% normal_load_vector
% load_transfer_front = abs(normal_load_vector(1) - normal_load_vector(2))/2
% load_transfer_rear = abs(normal_load_vector(3) - normal_load_vector(4))/2
% LLTD = load_transfer_front/(load_transfer_front+load_transfer_rear)

%F_y(left_tire, alpha, 0, F_z_L);


%% camber curve optimization
%roll angles and cambers to search through
roll_angle_vector = 0:0.25:1;
camber_vector = -5:0.1:5;

%arrays to store ideal camber curves 
wheel_displacement_matrix = zeros(4, numel(roll_angle_vector)); % 1-FL, 2-FR, 3-RL, 4-RR
normal_load_matrix = zeros(4, numel(roll_angle_vector)); % 1-FL, 2-FR, 3-RL, 4-RR
ideal_camber_matrix = zeros(4, numel(roll_angle_vector));

max_Fy_front_vector = zeros(3, numel(roll_angle_vector));
max_Fy_rear_vector = zeros(3, numel(roll_angle_vector));


%Loop through Roll Angles for front axle
for i = 1:numel(roll_angle_vector)
    
    %car roll angle
    roll_angle = roll_angle_vector(i);
    
    %calculate normal loads and wheel displacements at each wheel
    [normal_load_vector, wheel_displacement_vector] =...
        calcWheelForcesAndDisplacements(roll_angle, C);
    
    normal_load_vector
    
    normal_load_matrix(:,i) = normal_load_vector;
    wheel_displacement_matrix(:,i) = wheel_displacement_vector;
       
    %front axle, search through camber combinations for ideal combination
    for camber_FL = camber_vector
        for camber_FR = camber_vector
            [F_y_tot, F_y_L, F_y_R, M_x_L, M_x_R, alpha_val] = ....
                singleAxleCamberEvaluation(normal_load_vector(1), normal_load_vector(2), -camber_FL*0, camber_FR*0, tire);
            
            % if new configuration develops more lateral force, replace stored values
            if max_Fy_front_vector(1,i) < F_y_tot
                ideal_camber_matrix(1,i) = camber_FL;
                ideal_camber_matrix(2,i) = camber_FR;
                max_Fy_front_vector(1:3,i) = [F_y_tot;F_y_L;F_y_R];
            end
        end
    end

    %front axle, search through camber combinations for ideal combination
    for camber_RL = camber_vector
        for camber_RR = camber_vector
            [F_y_tot, F_y_L, F_y_R, M_x_L, M_x_R, alpha_val] = ...
                singleAxleCamberEvaluation(normal_load_vector(3), normal_load_vector(4), -camber_RL*0, camber_RR*0, tire);
            
            % if new configuration develops more lateral force, replace stored values 
            if max_Fy_rear_vector(1,i) < F_y_tot
                ideal_camber_matrix(3,i) = camber_RL;
                ideal_camber_matrix(4,i) = camber_RR;
                max_Fy_rear_vector(1:3,i) = [F_y_tot;F_y_L;F_y_R];
            end
        end
    end
end

%% plotting ideal results
subplot(4,1,1)
plot([-flip(roll_angle_vector) roll_angle_vector], [flip(ideal_camber_matrix(1,:)) ideal_camber_matrix(2,:)],'displayName', 'front camber ideal');
hold on;
title('Front Roll Camber');
xlabel('roll angle');
ylabel('camber (deg)');

subplot(4,1,2)
plot([-flip(roll_angle_vector) roll_angle_vector], [flip(ideal_camber_matrix(3,:)) ideal_camber_matrix(4,:)], 'displayName', 'rear camber ideal');
hold on;
title('Rear Roll Camber');
xlabel('roll angle');
ylabel('camber (deg)');
legend('Location','southeast');


% subplot(4,1,3)
% plot(roll_angle_vector, max_Fy_front_vector, 'displayName', 'front Fy ideal');
% hold on;
% plot(roll_angle_vector, max_Fy_rear_vector, 'displayName', 'rear Fy ideal');
% title('Lateral Force');
% xlabel('roll angle');
% ylabel('lateral force');
% legend('Location','southeast');
normal_load_matrix = 0*normal_load_matrix +1;

subplot(4,1,3);
plot(roll_angle_vector, max_Fy_front_vector(2,:)./normal_load_matrix(1,:), 'displayName', 'front left Fy');
hold on;
plot(roll_angle_vector, max_Fy_front_vector(3,:)./normal_load_matrix(2,:), 'displayName', 'front right Fy');
plot(roll_angle_vector, max_Fy_rear_vector(2,:)./normal_load_matrix(3,:), 'displayName', 'rear left Fy');
plot(roll_angle_vector, max_Fy_rear_vector(3,:)./normal_load_matrix(4,:), 'displayName', 'rear right Fy');
title('Lateral Force');
xlabel('roll angle');
ylabel('lateral force');
legend('Location','southeast');

subplot(4,1,4);
plot(roll_angle_vector, max_Fy_rear_vector(1,:)./(max_Fy_front_vector(1,:) + max_Fy_rear_vector(1,:))-0.53, 'displayName', 'ratio');
hold on;
%yline(C.weight_dist, 'displayName', '% of weight in the rear');
%ylim([0.4,0.6]);
title('% Lateral Force in Rear');
xlabel('roll angle');
ylabel('rear Fy / total Fy');
legend('Location','southeast');
