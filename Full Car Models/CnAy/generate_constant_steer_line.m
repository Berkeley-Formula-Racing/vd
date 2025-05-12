function [lat_accel_vec, yaw_accel_vec, lat_velocities, x_table_all] = generate_constant_steer_line(car, velocity, steer_angle, dataPointStability)

betaVal = deg2rad(dataPointStability(1));
counter = 1;

[x_table_ss, lat_accel, yaw_accel, lat_accel_guess] = max_constant_steer(velocity, steer_angle, car, -1);
lat_vel1 = x_table_ss{1, 'lat_vel'};
x0 = lat_accel_guess;
[x_table_ss, lat_accel, yaw_accel, lat_accel_guess] = max_constant_steer(velocity, steer_angle, car, 1);
lat_vel2 = x_table_ss{1, 'lat_vel'};

lat_velocities = linspace(lat_vel1, lat_vel2, 30);
lat_velocities = [lat_velocities, 20*tan(betaVal), 20*tan(betaVal+deg2rad(1))];

lat_accel_vec = zeros(size(lat_velocities));
yaw_accel_vec = zeros(size(lat_velocities));
x_table_all = cell(size(lat_velocities));

for i = lat_velocities
    [x_table_ss, lat_accel, yaw_accel, lat_accel_guess] = constant_steer(velocity, steer_angle, i, car, x0);
    x0 = lat_accel_guess;

    yaw_accel_vec(counter) = yaw_accel;
    lat_accel_vec(counter) = lat_accel;
    x_table_all{counter} = x_table_ss;
    counter = counter + 1;
end

end
