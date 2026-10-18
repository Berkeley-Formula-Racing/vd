%% Parameters
clear; clc; close all
setup_paths
carCell = carConfig(); % generate all cars to sim over
numCars = size(carCell,1);

car = carCell{1,1};

%%

% x: long_vel, lat_vel, yaw_rate
%x0 = [10 0.5 1]';
steer_angle = 7;

long_vel = 15;
lat_vel = 0;
for yaw_rate = linspace(-1.5,1.5,5)
    for lat_vel = linspace(-1,1,5)
        x0 = [long_vel,lat_vel,yaw_rate]';
        f = @(t,x) car.PhasePlaneODE(x,steer_angle);
        tspan = [0 10];
        [t,x] = ode45(f,tspan,x0);
        
        dxdt = zeros(size(x));
        for i = 1:numel(t)
            dxdt(i,:) = car.PhasePlaneODE(x(i,:),steer_angle);
        end
        
        beta = tan(x(:,2)./x(:,1))*180/pi;
        r = x(:,3);
        
        
        Cn = dxdt(:,3)*car.I_zz/(car.M*9.81*car.W_b);
        Ay = dxdt(:,2)/9.81;
        %plot(Ay,Cn);
        
        plot(beta,r);
        hold on
    end
end

