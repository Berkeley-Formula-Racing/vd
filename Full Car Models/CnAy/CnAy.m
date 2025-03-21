clear
setup_paths
carCell = carConfig(); % generate all cars to sim over
numCars = size(carCell, 1);
dataPointControl= [0, 0]; %beta, delta
dataPointStability= [1,0];
control = [];
stability = []; 
steer_angles = linspace(-20,20, 41);

velocity = 20; % m/s
velocityList = linspace(12,20, 3); %see how car handles through different speed corners
for i = 1:numCars

%for i = 1:5
    %velocity = velocityList(i);
    figure
    colors = ['r' 'k' 'b' 'g' 'm'];
    car = carCell{i, 1};
    %car = carCell{1,1};
    max_lat_accel = -inf; 
    corresponding_Cn = 0;
    max_lat_accel_at_zero_Cn = -inf; 
    for steer_angle = [-16 -8 -4 -2 0 2 4 8 16]
        [lat_accel, yaw_accel, lat_velocities] = generate_constant_steer_line(car, velocity, steer_angle, dataPointStability);
        Cn = yaw_accel * car.I_zz;
        % Cn = yaw_accel * car.I_zz / (car.M * 9.81 * car.W_b);
        beta_vals = (atan(lat_velocities/20))*180/pi; %lat_vel to beta values in deg
        
        if (steer_angle == dataPointStability(2))
            stability(i) = (Cn(31)-Cn(32));
            plot(lat_accel(31) / 9.81, Cn(31), 'go'); 
            plot(lat_accel(32) / 9.81, Cn(32), 'bo'); 
        end
        
        lat_accel(32) = []; lat_accel(31) = []; Cn(32) = []; Cn(31) = [];%remove extra data from plot
        plot(lat_accel / 9.81, Cn, colors(i))
        hold on

        % Update max lateral acceleration and corresponding Cn
        
        if max(lat_accel) > max_lat_accel
            max_lat_accel = max(lat_accel);
            corresponding_Cn = Cn(lat_accel == max(lat_accel));
        end
        % Update max lateral acceleration at Cn approximately 0
        if abs(Cn(end)) < 40 && max(lat_accel) > max_lat_accel_at_zero_Cn
            max_lat_accel_at_zero_Cn = max(lat_accel);
        end
        
    end
    
    for beta = [-8 -4 -2 -1 0 1 2 4 8]
        [lat_accel, yaw_accel] = generate_constant_beta_line(car, velocity, beta);
        Cn = yaw_accel * car.I_zz;
        %Cn = yaw_accel * car.I_zz / (car.M * 9.81 * car.W_b);
        plot(lat_accel / 9.81, Cn, colors(i))
        hold on
        
        if (beta == dataPointControl(1))
            pointIndex = find(steer_angles == dataPointControl(2));
            control(i) = (Cn(pointIndex+1)-Cn(pointIndex));
            plot(lat_accel(pointIndex) / 9.81, Cn(pointIndex), 'mo'); 
            plot(lat_accel(pointIndex+1) / 9.81, Cn(pointIndex+1), 'ro'); 
        end
        % Update max lateral acceleration and corresponding Cn
        if max(lat_accel) > max_lat_accel
            max_lat_accel = max(lat_accel);
            corresponding_Cn = Cn(lat_accel == max(lat_accel));
        end

        % Update max lateral acceleration at Cn approximately 0
        if abs(Cn(end)) < 40 && max(lat_accel) > max_lat_accel_at_zero_Cn
            max_lat_accel_at_zero_Cn = max(lat_accel);
        end
        
    end

    xlabel('Ay')
    ylabel('Cn')
    grid on
    title(['Car ' num2str(i) ' - Max Ay: ' num2str(max_lat_accel) ' m/s^2, Cn: ' num2str(corresponding_Cn) ', Max Ay at Cn=0: ' num2str(max_lat_accel_at_zero_Cn) ' m/s^2'])
end

stability
control