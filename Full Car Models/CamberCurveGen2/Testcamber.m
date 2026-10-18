%% setup

clc
close all
clear all

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
C.front_spring_roll_stiffness = 2860; %numbers from LLTD Doc
C.rear_spring_roll_stiffness = 2911;

%ARB Roll Stiffness (Currently on Short)
C.front_ARB_roll_stiffness = 0;
C.rear_ARB_roll_stiffness =986; %699, 1493

%Car Weight Distribution
C.weight_dist = 0.54; % percentage of weight in rear

%Total Car Weight
C.mass = 395 + 150; % lbs

%% camber curve optimization
%roll angles and cambers to search through
roll_angle_vector = 0:0.05:1;
camber_vector = -5:0.5:5;

%arrays to store ideal camber curves 
wheel_displacement_matrix = zeros(4, numel(roll_angle_vector)); % 1-FL, 2-FR, 3-RL, 4-RR
ideal_camber_matrix = zeros(4, numel(roll_angle_vector));

max_Fy_front_vector = zeros(1, numel(roll_angle_vector));
max_Fy_rear_vector = zeros(1, numel(roll_angle_vector));


%Loop through Roll Angles for front axle
for i = 1:numel(roll_angle_vector)
    
    %car roll angle
    roll_angle = roll_angle_vector(i)
    
    %calculate normal loads and wheel displacements at each wheel
    [normal_load_vector, wheel_displacement_vector] =...
        calcWheelForcesAndDisplacements(roll_angle, C);
    
    wheel_displacement_matrix(:,i) = wheel_displacement_vector;
       
    %front axle, search through camber combinations for ideal combination
    for camber_FL = camber_vector
        for camber_FR = camber_vector
            [F_y_tot, F_y_L, F_y_R, M_x_L, M_x_R, alpha_val] = ....
                singleAxleCamberEvaluation(normal_load_vector(1), normal_load_vector(2), -camber_FL, camber_FR, tire);
            
            % if new configuration develops more lateral force, replace stored values
            if max_Fy_front_vector(i) < F_y_tot
                ideal_camber_matrix(1,i) = camber_FL;
                ideal_camber_matrix(2,i) = camber_FR;
                max_Fy_front_vector(i) = F_y_tot;
            end
        end
    end

    %front axle, search through camber combinations for ideal combination
    for camber_RL = camber_vector
        for camber_RR = camber_vector
            [F_y_tot, F_y_L, F_y_R, M_x_L, M_x_R, alpha_val] = ...
                singleAxleCamberEvaluation(normal_load_vector(3), normal_load_vector(4), -camber_RL, camber_RR, tire);
            
            % if new configuration develops more lateral force, replace stored values 
            if max_Fy_rear_vector(i) < F_y_tot
                ideal_camber_matrix(3,i) = camber_RL;
                ideal_camber_matrix(4,i) = camber_RR;
                max_Fy_rear_vector(i) = F_y_tot;
            end
        end
    end
end

%% plotting ideal results
subplot(3,1,1)
plot([-flip(roll_angle_vector) roll_angle_vector], [flip(ideal_camber_matrix(1,:)) ideal_camber_matrix(2,:)],'displayName', 'Front Camber Ideal');
hold on;
title('Front Roll Camber');
xlabel('Roll Angle (deg)');
ylabel('Camber (deg)');

subplot(3,1,2)
plot([-flip(roll_angle_vector) roll_angle_vector], [flip(ideal_camber_matrix(3,:)) ideal_camber_matrix(4,:)], 'displayName', 'Rear Camber Ideal');
hold on;
title('Rear Roll Camber');
xlabel('Roll Angle (deg)');
ylabel('Camber (deg)');

subplot(3,1,3)
plot(roll_angle_vector, max_Fy_front_vector + max_Fy_rear_vector, 'displayName', 'Fy Ideal');
hold on;
%plot(roll_angle_vector, max_Fy_rear_vector, 'displayName', 'rear Fy ideal');
title('Max Lateral Force');
xlabel('Roll Angle (deg)');
ylabel('Lateral Force');

% subplot(4,1,4)
% plot(roll_angle_vector, max_Fy_rear_vector./(max_Fy_front_vector + max_Fy_rear_vector), 'displayName', 'ideal ratio');
% hold on;
% title('% Lateral Force in Rear');
% xlabel('roll angle');
% ylabel('rear Fy / total Fy');


%% compare to existing

existing_camber_matrix = modifiedRollCamberCurves(roll_angle_vector);
existing_Fy_front_vector = zeros(1,numel(roll_angle_vector));
existing_Fy_rear_vector = zeros(1,numel(roll_angle_vector));

for i = 1:numel(roll_angle_vector)
    %car roll angle
    roll_angle = roll_angle_vector(i);
    
    %calculate normal loads and wheel displacements at each wheel
    [normal_load_vector, wheel_displacement_vector] =...
        calcWheelForcesAndDisplacements(roll_angle, C);
    
    
    [existing_Fy_front_vector(i), ~, ~, ~, ~, ~] = ...
        singleAxleCamberEvaluation(normal_load_vector(1), normal_load_vector(2), -existing_camber_matrix(1,i), existing_camber_matrix(2,i), tire);
    
    [existing_Fy_rear_vector(i), ~, ~, ~, ~, ~] = ...
        singleAxleCamberEvaluation(normal_load_vector(3), normal_load_vector(4), -existing_camber_matrix(3,i), existing_camber_matrix(4,i), tire);
end

%% plotting comparison camber curves

% front camber curve
subplot(3,1,1);
plot([-flip(roll_angle_vector) roll_angle_vector], [flip(existing_camber_matrix(1,:)) existing_camber_matrix(2,:)], '--', 'displayName', 'Front Camber Modified');
legend('Location','southeast');

% rear camber curve
subplot(3,1,2);
plot([-flip(roll_angle_vector) roll_angle_vector], [flip(existing_camber_matrix(3,:)) existing_camber_matrix(4,:)], '--', 'displayName', 'Rear Camber Modified');
legend('Location','southeast');

% Fy plots
subplot(3,1,3)
plot(roll_angle_vector, existing_Fy_front_vector + existing_Fy_rear_vector,'--', 'displayName', 'Fy Modified');
%plot(roll_angle_vector, existing_Fy_rear_vector,'--', 'displayName', 'rear Fy existing');
%ylim([0,425]);
legend('Location','southeast');

% Car balance plot
% subplot(4,1,4)
% plot(roll_angle_vector, existing_Fy_rear_vector./(existing_Fy_rear_vector + existing_Fy_front_vector),'--', 'displayName', 'existing ratio');
% yline(C.weight_dist, 'displayName', '% of weight in the rear');
% ylim([0.4,0.6]);
% legend('Location','southeast');

%% compare to existing

existing_camber_matrix = existingRollCamberCurves(roll_angle_vector);
existing_Fy_front_vector = zeros(1,numel(roll_angle_vector));
existing_Fy_rear_vector = zeros(1,numel(roll_angle_vector));

for i = 1:numel(roll_angle_vector)
    %car roll angle
    roll_angle = roll_angle_vector(i);
    
    %calculate normal loads and wheel displacements at each wheel
    [normal_load_vector, wheel_displacement_vector] =...
        calcWheelForcesAndDisplacements(roll_angle, C);
    
    
    [existing_Fy_front_vector(i), ~, ~, ~, ~, ~] = ...
        singleAxleCamberEvaluation(normal_load_vector(1), normal_load_vector(2), -existing_camber_matrix(1,i), existing_camber_matrix(2,i), tire);
    
    [existing_Fy_rear_vector(i), ~, ~, ~, ~, ~] = ...
        singleAxleCamberEvaluation(normal_load_vector(3), normal_load_vector(4), -existing_camber_matrix(3,i), existing_camber_matrix(4,i), tire);
end

%% plotting comparison camber curves

% front camber curve
subplot(3,1,1);
plot([-flip(roll_angle_vector) roll_angle_vector], [flip(existing_camber_matrix(1,:)) existing_camber_matrix(2,:)], '--', 'displayName', 'Front Camber Existing');
legend('Location','southeast');

% rear camber curve
subplot(3,1,2);
plot([-flip(roll_angle_vector) roll_angle_vector], [flip(existing_camber_matrix(3,:)) existing_camber_matrix(4,:)], '--', 'displayName', 'Rear Camber Existing');
legend('Location','southeast');

% Fy plots
subplot(3,1,3)
plot(roll_angle_vector, existing_Fy_front_vector + existing_Fy_rear_vector,'--', 'displayName', 'Fy Existing');
%plot(roll_angle_vector, existing_Fy_rear_vector,'--', 'displayName', 'rear Fy existing');
%ylim([0,425]);
legend('Location','southeast');

% Car balance plot
% subplot(4,1,4)
% plot(roll_angle_vector, existing_Fy_rear_vector./(existing_Fy_rear_vector + existing_Fy_front_vector),'--', 'displayName', 'existing ratio');
% yline(C.weight_dist, 'displayName', '% of weight in the rear');
% ylim([0.4,0.6]);
% legend('Location','southeast')
