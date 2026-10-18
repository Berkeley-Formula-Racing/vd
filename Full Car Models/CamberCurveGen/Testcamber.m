%% setup

clc
close all
clear

%% taken from carConfig

%change Fx parameters currently here - Aeden.
load('Fx_combined_parameters_run38_30.mat'); % F_x combined magic formula parameters
Fx_parameters = cell2mat(Xbestcell);
load('Fy_pure_parameters_1965run24.mat'); % F_y combined magic formula parameters
%load('Lapsim_Fy_combined_parameters_1965run15.mat');
Fy_parameters = cell2mat(Xbestcell);

friction_scaling_factor = 1.05*0.55;
p_i = 12; % tire pressure

%initialize tire object, used for camber evaluation
tire = BDITire2(0,p_i,Fx_parameters,Fy_parameters,friction_scaling_factor);

%Vehicle Width in inches
C.front_width = 47;
C.rear_width = 47;
C.wheelbase = 62;

%Spring Roll Stiffness
C.front_spring_roll_stiffness = 2860; %numbers from LLTD Doc
C.rear_spring_roll_stiffness = 2911;
C.front_bump_stiffness = 2*191.66;%wheel rates
C.rear_bump_stiffness = 2*196.02; 

%ARB Roll Stiffness (Currently on Short)
C.front_ARB_roll_stiffness = 0;
C.rear_ARB_roll_stiffness =986; %699, 1493

%Car Weight Distribution
C.weight_dist = 0.54; % percentage of weight in rear

%Total Car Weight
C.mass = (395 + 150)*4.44; % lbs to N

%anti effects
front_anti_dive = .545;
rear_anti_lift = 1.02;
rear_anti_squat = .457;

%% camber curve optimization
%roll angles and cambers to search through
roll_angle_vector = 0:0.05:1;
%-.143431, .5
pitch_angle_vector = linspace(-.143431,.5,100);
camber_vector = -10:.5:10;


%arrays to store ideal camber curves 
wheel_displacement_matrix = zeros(4, numel(roll_angle_vector)); % 1-FL, 2-FR, 3-RL, 4-RR
ideal_camber_matrix = zeros(5, numel(roll_angle_vector));

max_Fy_front_vector = zeros(1, numel(roll_angle_vector));


max_Fy_rear_vector = zeros(1, numel(roll_angle_vector));
max_Fx_vector = zeros(1,numel(pitch_angle_vector));


%Loop through Roll Angles for front axle
for i = 1:numel(roll_angle_vector)
   
    %car roll angle
    roll_angle = roll_angle_vector(i);
    %pitch_angle = pitch_angle_vector(i);
    
    %calculate normal loads and wheel displacements at each wheel
    [normal_load_vector, wheel_displacement_vector] =...
        calcWheelForcesAndDisplacements(roll_angle, C);
    %[pitch_load_vector, pitch_displacement_vector] = calcWheelForcesAndDisplacements_pitch(pitch_angle,C);

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


    %rear axle, search through camber combinations for ideal combination
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
subplot(2,2,1)
plot([-flip(roll_angle_vector) roll_angle_vector], [flip(ideal_camber_matrix(1,:)) ideal_camber_matrix(2,:)],'displayName', 'Front Camber Ideal');
hold on;
title('Front Roll Camber');
xlabel('Roll Angle (deg)');
ylabel('Camber (deg)');

subplot(2,2,2)
plot([-flip(roll_angle_vector) roll_angle_vector], [flip(ideal_camber_matrix(3,:)) ideal_camber_matrix(4,:)], 'displayName', 'Rear Camber Ideal');
hold on;
title('Rear Roll Camber');
xlabel('Roll Angle (deg)');
ylabel('Camber (deg)');

subplot(2,2,3)
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
subplot(2,2,1);
plot([-flip(roll_angle_vector) roll_angle_vector], [flip(existing_camber_matrix(1,:)) existing_camber_matrix(2,:)], '--', 'displayName', 'Front Camber Modified');
legend('Location','southeast');

% rear camber curve
subplot(2,2,2);
plot([-flip(roll_angle_vector) roll_angle_vector], [flip(existing_camber_matrix(3,:)) existing_camber_matrix(4,:)], '--', 'displayName', 'Rear Camber Modified');
legend('Location','southeast');

% Fy plots
subplot(2,2,3)
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
    %pitch_angle = pitch_angle_vector(i);
    
    %calculate normal loads and wheel displacements at each wheel
    [normal_load_vector, wheel_displacement_vector] =...
        calcWheelForcesAndDisplacements(roll_angle, C);
    %[pitch_load_vector, wheel_disp_long] = calcWheelForcesAndDisplacements_pitch(pitch_angle,C));
    %[existing_Fx_vector(i)] = BicycleCamberEvaluation(pitch_load_vector(1),pitch_load_vector(3),)
    
    [existing_Fy_front_vector(i), ~, ~, ~, ~, ~] = ...
        singleAxleCamberEvaluation(normal_load_vector(1), normal_load_vector(2), -existing_camber_matrix(1,i), existing_camber_matrix(2,i), tire);
    
    [existing_Fy_rear_vector(i), ~, ~, ~, ~, ~] = ...
        singleAxleCamberEvaluation(normal_load_vector(3), normal_load_vector(4), -existing_camber_matrix(3,i), existing_camber_matrix(4,i), tire);
end
%% put stuff here
%get ideal camber vs disp arrays
idealCamber_vs_disp_front = [flip(ideal_camber_matrix(1,:)), ideal_camber_matrix(2,:); ...
                              flip(wheel_displacement_matrix(1,:)),wheel_displacement_matrix(2,:)];
idealCamber_vs_disp_rear = [flip(ideal_camber_matrix(3,:)), ideal_camber_matrix(4,:); ...
                              flip(wheel_displacement_matrix(3,:)),wheel_displacement_matrix(4,:)];
idealCamber_Fx_vector = zeros(1,numel(pitch_angle_vector));
[no_repeat,ia] = unique(idealCamber_vs_disp_front(2,:),'stable');
idealCamber_vs_disp_front= [idealCamber_vs_disp_front(1,ia);no_repeat];

[no_repeat,ia] = unique(idealCamber_vs_disp_rear(2,:),'stable');
idealCamber_vs_disp_rear= [idealCamber_vs_disp_rear(1,ia);no_repeat];

%get existing camber vs displacement arrays
existingCamber_vs_disp_front = [flip(existing_camber_matrix(1,:)), existing_camber_matrix(2,:); ...
                              flip(wheel_displacement_matrix(1,:)),wheel_displacement_matrix(2,:)];
existingCamber_vs_disp_rear = [flip(existing_camber_matrix(3,:)), existing_camber_matrix(4,:); ...
                              flip(wheel_displacement_matrix(3,:)),wheel_displacement_matrix(4,:)];
existingCamber_Fx_vector = zeros(1,numel(pitch_angle_vector));

[no_repeat,ia] = unique(existingCamber_vs_disp_front(2,:),'stable');
existingCamber_vs_disp_front= [existingCamber_vs_disp_front(1,ia);no_repeat];

[no_repeat,ia] = unique(existingCamber_vs_disp_rear(2,:),'stable');
existingCamber_vs_disp_rear= [existingCamber_vs_disp_rear(1,ia);no_repeat];

TLV = zeros(4,numel(pitch_angle_vector));

for i = 1:numel(pitch_angle_vector)
    pitch_angle = pitch_angle_vector(i);

    if pitch_angle < 0
        state = "brake";
    else
        state = 'accel';
    end
    [pitch_load_vector, pitch_displacement_vector] = calcWheelForcesAndDisplacements_pitch(pitch_angle,C);
    index_f = 1;
    index_r = 1;
    minDiff_r = 100;
    minDiff_f = 100; %just so it gets trumped fast
   %{ 
for j = 1:length(idealCamber_vs_disp_rear)
        min_temp_f = min(abs(idealCamber_vs_disp_front(2,j)-pitch_displacement_vector(1)),minDiff_f);

        min_temp_r = min(abs(idealCamber_vs_disp_rear(2,j)-pitch_displacement_vector(3)),minDiff_r);
        if min_temp_f ~= minDiff_f
            minDiff_f = min_temp_f;
            index_f = j;  
        elseif min_temp_r ~= minDiff_r
            minDiff_r = min_temp_r;
            index_r = j;
        end
end 
   %}

    %interpolation method can struggle for pitch displacement that is
    %greater than disp due to roll.
    existing_f = interp1(existingCamber_vs_disp_front(2,:),existingCamber_vs_disp_front(1,:),pitch_displacement_vector(1));
    existing_r = interp1(existingCamber_vs_disp_rear(2,:),existingCamber_vs_disp_rear(1,:),pitch_displacement_vector(3));

    ideal_f = interp1(idealCamber_vs_disp_front(2,:),idealCamber_vs_disp_front(1,:),pitch_displacement_vector(1));
    ideal_r = interp1(idealCamber_vs_disp_rear(2,:),idealCamber_vs_disp_rear(1,:),pitch_displacement_vector(3));


    %calc long accel w ideal camber
    [Fx_tot, ~, ~, ~, ~, ~] = BicycleCamberEvaluation(pitch_load_vector(1),pitch_load_vector(3), ...
        ideal_f,ideal_r,tire,state);
    idealCamber_Fx_vector(i) = Fx_tot;
    TLV(1:4,i)= pitch_load_vector; %so greg can look at the loads on each tire
    
    %calc long accel w 0 camber (ideal for long accel)
    [Fx_tot, Fx_front, Fx_rear, M_x_L, M_x_R, kappa_val] = BicycleCamberEvaluation(pitch_load_vector(1),pitch_load_vector(3),0,0,tire,state);
    max_Fx_vector(i) = Fx_tot;

    %calc long accel w existing camber
    [Fx_tot, ~, ~, ~, ~, ~] = BicycleCamberEvaluation(pitch_load_vector(1),pitch_load_vector(3), ...
        existing_f,existing_r,tire,state);
    existingCamber_Fx_vector(i) = Fx_tot;
    
end



%% plotting comparison camber curves

% front camber curve (right tire)
subplot(2,2,1);
plot([-flip(roll_angle_vector) roll_angle_vector], [flip(existing_camber_matrix(1,:)) existing_camber_matrix(2,:)], '--', 'displayName', 'Front Camber Existing');
legend('Location','southeast');

% rear camber curve (right tire)
subplot(2,2,2);
plot([-flip(roll_angle_vector) roll_angle_vector], [flip(existing_camber_matrix(3,:)) existing_camber_matrix(4,:)], '--', 'displayName', 'Rear Camber Existing');
legend('Location','southeast');

% Fy plots
subplot(2,2,3)
plot(roll_angle_vector, existing_Fy_front_vector + existing_Fy_rear_vector,'--', 'displayName', 'Fy Existing');
%plot(roll_angle_vector, existing_Fy_rear_vector,'--', 'displayName', 'rear Fy existing');
%ylim([0,425]);
legend('Location','southeast');

subplot(2,2,4)
plot(pitch_angle_vector, max_Fx_vector,pitch_angle_vector,idealCamber_Fx_vector)%,pitch_angle_vector,zeros(length(idealCamber_Fx_vector)))

% Car balance plot
% subplot(4,1,4)
% plot(roll_angle_vector, existing_Fy_rear_vector./(existing_Fy_rear_vector + existing_Fy_front_vector),'--', 'displayName', 'existing ratio');
% yline(C.weight_dist, 'displayName', '% of weight in the rear');
% ylim([0.4,0.6]);
% legend('Location','southeast')
