function [normal_load_vector, wheel_displacement_vector] = calcWheelForcesAndDisplacements_pitch(pitch_angle, C)
    % negative pitch = nose down
    
    % 1-FL, 2-FR, 3-RL, 4-RR
    wheel_displacement_vector = zeros(4,1);
    normal_load_vector = zeros(4,1);

    front_anti_lift = 0;
    front_anti_dive = .545;
    rear_anti_lift = 1.02;
    rear_anti_squat = .457;

    if pitch_angle<0
        %disp_f/disp_r
        disp_ratio = C.rear_bump_stiffness*(1-front_anti_dive)/(C.front_bump_stiffness*(1-rear_anti_lift));
    else
        disp_ratio = -C.rear_bump_stiffness*(1-front_anti_lift)/(C.front_bump_stiffness*(1-rear_anti_squat));
    end
    

    %The distance the rear wheels move as a function of pitch angle
    wheel_displacement_vector(3) = C.wheelbase*tand(pitch_angle)/(1-disp_ratio);
    wheel_displacement_vector(4) = C.wheelbase*tand(pitch_angle)/(1-disp_ratio);

    %The distance the front wheels move as a function of pitch angle
    wheel_displacement_vector(1) = wheel_displacement_vector(3)*disp_ratio; %calculated by wheelRateFront/wheelRateTotal
    wheel_displacement_vector(2) = wheel_displacement_vector(3)*disp_ratio;

    

    if pitch_angle<0
        %Normal Load on front wheels as a function of pitch angle
        normal_load_vector(1) = (1 - C.weight_dist) * (C.mass/2) - ...
            wheel_displacement_vector(1)/(1-front_anti_dive)*C.front_bump_stiffness;
        normal_load_vector(2) = (1 - C.weight_dist) * (C.mass/2) - ...
            wheel_displacement_vector(2)/(1-front_anti_dive)*C.front_bump_stiffness;
        
        %Normal Load on rear wheels as a function of pitch angle
        normal_load_vector(3) = (C.weight_dist) * (C.mass/2) + ...
            wheel_displacement_vector(3)/(1-rear_anti_lift)*C.rear_bump_stiffness;
        normal_load_vector(4) = (C.weight_dist) * (C.mass/2) + ...
            wheel_displacement_vector(4)/(1-rear_anti_lift)*C.rear_bump_stiffness;

    else
        %Normal Load on front wheels as a function of pitch angle
        normal_load_vector(1) = (1 - C.weight_dist) * (C.mass/2) + ...
            wheel_displacement_vector(1)/(1-front_anti_lift)*C.front_bump_stiffness;
        normal_load_vector(2) = (1 - C.weight_dist) * (C.mass/2) + ...
            wheel_displacement_vector(2)/(1-front_anti_lift)*C.front_bump_stiffness;
        
        %Normal Load on rear wheels as a function of pitch angle
        normal_load_vector(3) = (C.weight_dist) * (C.mass/2) + ...
            wheel_displacement_vector(3)/(1-rear_anti_squat)*C.rear_bump_stiffness;
        normal_load_vector(4) = (C.weight_dist) * (C.mass/2) + ...
            wheel_displacement_vector(4)/(1-rear_anti_squat)*C.rear_bump_stiffness;


    end
    %disp(normal_load_vector(1)+normal_load_vector(3))
    %disp(wheel_displacement_vector(1))
    %disp(wheel_displacement_vector(3))
    %rear_load_transfer = pitch_angle / C.rear_width * ( C.rear_spring_roll_stiffness + C.rear_ARB_roll_stiffness )
end