function [lat_accel_vec,yaw_accel_vec] = generate_constant_beta_line(car,velocity,beta)

counter = 1;

[x_table_ss,lat_accel,yaw_accel,lat_accel_guess] = max_constant_beta(velocity,beta,car,-1);
steer_angle1 = x_table_ss{1,'steer_angle'};
x0 = lat_accel_guess;
[x_table_ss,lat_accel,yaw_accel,lat_accel_guess] = max_constant_beta(velocity,beta,car,1);
steer_angle2 = x_table_ss{1,'steer_angle'};

steer_angles = linspace(-20,20,30);
lat_accel_vec = zeros(size(steer_angles));
yaw_accel_vec = zeros(size(steer_angles));

for i = steer_angles  
    lat_vel = velocity*tand(beta);
    [x_table_ss,lat_accel,yaw_accel,lat_accel_guess] = constant_steer(velocity,i,lat_vel,car); %,x0);
    
%     exitflag = x_table_ss{1,'exitflag'};
%     if exitflag ~= 1 
%         break
%     end
    
    x0 = lat_accel_guess;
    yaw_accel_vec(counter) = yaw_accel; 
    lat_accel_vec(counter) = lat_accel;
    counter = counter+1;
end

end
