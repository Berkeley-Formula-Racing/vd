function [lat_accel_vec,yaw_accel_vec] = generate_constant_steer_line(car,velocity,steer_angle)

counter = 1;

[x_table_ss,lat_accel,yaw_accel,lat_accel_guess] = max_constant_steer(velocity,steer_angle,car,-1);
lat_vel1 = x_table_ss{1,'lat_vel'};
x0 = lat_accel_guess;
[x_table_ss,lat_accel,yaw_accel,lat_accel_guess] = max_constant_steer(velocity,steer_angle,car,1);
lat_vel2 = x_table_ss{1,'lat_vel'};

lat_velocities = linspace(lat_vel1,lat_vel2,15);
lat_accel_vec = zeros(size(lat_velocities));
yaw_accel_vec = zeros(size(lat_velocities));

for i = lat_velocities  
    [x_table_ss,lat_accel,yaw_accel,lat_accel_guess] = constant_steer(velocity,steer_angle,i,car,x0);
    
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
