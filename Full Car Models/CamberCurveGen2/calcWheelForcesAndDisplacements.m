function [normal_load_vector, wheel_displacement_vector] = calcWheelForcesAndDisplacements(roll_angle, C)
    % + roll angle => car rolls to the left
    
    % 1-FL, 2-FR, 3-RL, 4-RR
    wheel_displacement_vector = zeros(4,1);
    normal_load_vector = zeros(4,1);
    
    %The distance the front wheels move as a function of roll angle
    wheel_displacement_vector(1) = (C.front_width*sind(roll_angle))/2;
    wheel_displacement_vector(2) = -(C.front_width*sind(roll_angle))/2;

    %The distance the rear wheels move as a function of roll angle
    wheel_displacement_vector(3) = (C.rear_width*sind(roll_angle))/2;
    wheel_displacement_vector(4) = -(C.rear_width*sind(roll_angle))/2;
    

    %Normal Load on front wheels as a function of roll angle
    normal_load_vector(1) = (1 - C.weight_dist) * (C.mass/2) + ...
        roll_angle / C.front_width * ( C.front_spring_roll_stiffness + C.front_ARB_roll_stiffness );
    normal_load_vector(2) = (1 - C.weight_dist) * (C.mass/2) - ...
        roll_angle / C.front_width * ( C.front_spring_roll_stiffness + C.front_ARB_roll_stiffness );
    
    %Normal Load on rear wheels as a function of roll angle
    normal_load_vector(3) = (C.weight_dist) * (C.mass/2) + ...
        roll_angle / C.rear_width * ( C.rear_spring_roll_stiffness + C.rear_ARB_roll_stiffness );
    normal_load_vector(4) = (C.weight_dist) * (C.mass/2) - ...
        roll_angle / C.rear_width * ( C.rear_spring_roll_stiffness + C.rear_ARB_roll_stiffness );
    
    %rear_load_transfer = roll_angle / C.rear_width * ( C.rear_spring_roll_stiffness + C.rear_ARB_roll_stiffness )
end







