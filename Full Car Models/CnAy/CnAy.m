%% Parameters
clear
setup_paths
carCell = carConfig(); % generate all cars to sim over
numCars = size(carCell,1);
car = carCell{1,1};

%%
velocity = 20; % m/s

for steer_angle = [-16 -8 -4 -2 0 2 4 8 16]
    [lat_accel,yaw_accel] = generate_constant_steer_line(car,velocity,steer_angle);
    Cn = yaw_accel*car.I_zz/(car.M*9.81*car.W_b);
    plot(lat_accel/9.81,Cn,'k')
    hold on
end

for beta = [-8 -4 -2 -1 0 1 2 4 8]
    [lat_accel,yaw_accel] = generate_constant_beta_line(car,velocity,beta);
    Cn = yaw_accel*car.I_zz/(car.M*9.81*car.W_b);
    plot(lat_accel/9.81,Cn,'k')
    hold on
end

xlabel('Ay')
ylabel('Cn')
grid on
