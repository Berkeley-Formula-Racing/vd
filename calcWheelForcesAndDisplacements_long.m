function [normal_load_vector, wheel_displacement_vector] = calcWheelForcesAndDisplacements_long(pitch_angle, C)
    % negative pitch = nose down
    
    % 1-FL, 2-FR, 3-RL, 4-RR
    wheel_displacement_vector = zeros(4,1);
    normal_load_vector = zeros(4,1);
    
    %The distance the front wheels move as a function of roll angle
    wheel_displacement_vector(1) = C.wheelbase*ratio*sind(pitch_angle);
    wheel_displacement_vector(2) = C.wheelbase*ratio*sind(pitch_angle);

    %The distance the rear wheels move as a function of roll angle
    wheel_displacement_vector(3) = -(C.wheelbase*(1-ratio)*sind(pitch_angle));
    wheel_displacement_vector(4) = -(C.wheelbase*(1-ratio)*sind(pitch_angle));
    

    %Normal Load on front wheels as a function of roll angle
    normal_load_vector(1) = (1 - C.weight_dist) * (C.mass/2) + ...
        pitch_angle / C.front_width * ( C.front_spring_roll_stiffness + C.front_ARB_roll_stiffness );
    normal_load_vector(2) = (1 - C.weight_dist) * (C.mass/2) - ...
        pitch_angle / C.front_width * ( C.front_spring_roll_stiffness + C.front_ARB_roll_stiffness );
    
    %Normal Load on rear wheels as a function of roll angle
    normal_load_vector(3) = (C.weight_dist) * (C.mass/2) + ...
        pitch_angle / C.rear_width * ( C.rear_spring_roll_stiffness + C.rear_ARB_roll_stiffness );
    normal_load_vector(4) = (C.weight_dist) * (C.mass/2) - ...
        pitch_angle / C.rear_width * ( C.rear_spring_roll_stiffness + C.rear_ARB_roll_stiffness );
    
    %rear_load_transfer = pitch_angle / C.rear_width * ( C.rear_spring_roll_stiffness + C.rear_ARB_roll_stiffness )
end