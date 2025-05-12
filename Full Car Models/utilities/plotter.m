function [] = plotter(car,g_g_vel,plot_choice)
% plotting gg-diagram variations for car

max_vel = car.max_vel;

long_g_accel = car.longAccelLookup(:,1)'/9.81;
lat_g_accel = car.longAccelLookup(:,2)'/9.81;
vel_accel = car.longAccelLookup(:,3)';

long_g_braking = car.longDecelLookup(:,1)'/9.81;
lat_g_braking = car.longDecelLookup(:,2)'/9.81;
vel_braking = car.longDecelLookup(:,3)';
    

if plot_choice(1)
    % g-g diagram for different velocities, scatter plot    
    figure
    
    
    scatter3([lat_g_accel -lat_g_accel -lat_g_braking lat_g_braking],...
        [long_g_accel long_g_accel long_g_braking long_g_braking],...
        [vel_accel vel_accel vel_braking vel_braking],'b');
    

end

if plot_choice(2)
    % g-g diagram using crust (alpha shape-based surface)
        figure

    % Prepare data
    x = [lat_g_accel, -lat_g_accel, -lat_g_braking, lat_g_braking]';
    y = [long_g_accel, long_g_accel, long_g_braking, long_g_braking]';
    z = [vel_accel, vel_accel, vel_braking, vel_braking]';

    % Create alpha shape and crust
    shp = alphaShape(x, y, z, 1.5);
    [tri, pts] = boundaryFacets(shp);

    % Plot with edge lines
    trisurf(tri, pts(:,1), pts(:,2), pts(:,3), ...
            'EdgeColor', 'k', ...    % black grid lines
            'FaceColor', 'interp', ...
            'FaceAlpha', 1);         % adjust alpha if needed

    xlabel('Lateral G');
    ylabel('Longitudinal G');
    zlabel('Velocity (m/s)');
    title('Velocity Dependent G-G Diagram Surface');
    view(3);
end


if plot_choice(3)
    figure
    
    x = lat_g_accel;
    y = vel_accel;
    z = long_g_accel;

    F_accel = scatteredInterpolant([x' y'],z');

    [Xq,Yq] = meshgrid(-2:0.05:2, 5:0.05:max_vel);
    Vq = F_accel(Xq,Yq);
    mesh(Xq,Yq,Vq);
    title('Max Accel (Scattered Interpolant)','FontSize',18)        
    xlabel('Lat G')
    ylabel('Velocity')
    zlabel('Long G')

    hold on
    scatter3(x,y,z);
    xlim([0 2]);
    zlim([0 1.5]);
    hold off
end

if plot_choice(4)
    figure
    
    x = lat_g_braking;
    y = vel_braking;
    z = long_g_braking;

    F_braking = scatteredInterpolant([x' y'],z');

    [Xq,Yq] = meshgrid(-2:0.05:2, 5:0.05:max_vel);
    Vq = F_braking(Xq,Yq);
    mesh(Xq,Yq,Vq);
    title('Max Braking (Scattered Interpolant)','FontSize',18)    
    xlabel('Lat G')
    ylabel('Velocity')
    zlabel('Long G')

    hold on
    scatter3(x,y,z);
    xlim([0 2]);
    zlim([-1.5 0]);
    hold off
end

if plot_choice(5)
    % === Load real data ===
    data = readtable("autocross_2.csv");

    % Convert wheel speeds from mph to m/s
    mph_to_mps = 0.44704;
    real_vel = mean([data.WheelSpdFL, data.WheelSpdFR, ...
                     data.WheelSpdRL, data.WheelSpdRR], 2) * mph_to_mps;

    real_lat_g = data.AccelY;  % assumed in g
    real_long_g = data.AccelX; % assumed in g
    raw_lat_g = data.AccelY;
    raw_long_g = data.AccelX;
    raw_vel = mean([data.WheelSpdFL, data.WheelSpdFR, ...
                    data.WheelSpdRL, data.WheelSpdRR], 2) * 0.44704;
    
    % Filter out points beyond the 1st and 99th percentiles in each signal
    lat_bounds = prctile(raw_lat_g, [1 99]);
    long_bounds = prctile(raw_long_g, [1 99]);
    vel_bounds = prctile(raw_vel, [1 99]);
    
    valid = raw_lat_g >= lat_bounds(1) & raw_lat_g <= lat_bounds(2) & ...
            raw_long_g >= long_bounds(1) & raw_long_g <= long_bounds(2) & ...
            raw_vel >= vel_bounds(1) & raw_vel <= vel_bounds(2);
    
    % Cleaned signals
    real_lat_g = raw_lat_g(valid);
    real_long_g = raw_long_g(valid);
    real_vel = raw_vel(valid);
    % === Automatically define velocity slices ===
    bin_edges = 7:7:car.max_vel;             % e.g. slices from 5 to max in 5 m/s bins
    bin_centers = bin_edges(1:end-1) + 2.5;  % midpoint of each bin

    for i = 1:length(bin_centers)
        v_min = bin_edges(i);
        v_max = bin_edges(i+1);
        v_mid = bin_centers(i);

        figure

        % Simulated car data at this velocity
        idx_accel = vel_accel >= v_min & vel_accel < v_max;
        idx_brake = vel_braking >= v_min & vel_braking < v_max;

        scatter([lat_g_accel(idx_accel) lat_g_braking(idx_brake) ...
                 -lat_g_braking(idx_brake) -lat_g_accel(idx_accel)], ...
                [long_g_accel(idx_accel) long_g_braking(idx_brake) ...
                 long_g_braking(idx_brake) long_g_accel(idx_accel)], ...
                'DisplayName', sprintf('Sim @ %.1f m/s', v_mid));
        hold on        

        % Real data in same bin
        idx_real = real_vel >= v_min & real_vel < v_max;
        scatter(real_lat_g(idx_real), real_long_g(idx_real), 10, 'r', 'filled', ...
            'DisplayName', sprintf('Real @ ~%.1f m/s', v_mid));

        title(sprintf('G-G Diagram at %.1f m/s', v_mid), 'FontSize', 18)
        xlabel('Lat G', 'FontSize', 15)
        ylabel('Long G', 'FontSize', 15)

        legend('show')
        grid on
        axis equal
        xlim([-2.5 2.5])
        ylim([-2.5 2.5])
        hold off
    end
end
