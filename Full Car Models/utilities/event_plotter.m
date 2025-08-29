function [] = event_plotter(comp,plot_choice)

if plot_choice(1)
    figure
    hold on
    title('Autocross Track')
    plot(comp.autocross_track(1,:),comp.autocross_track(2,:))
    scatter(comp.autocross_track)
    xlabel('Distance (m)')
    ylabel('Curvature (1/m)')
end

if plot_choice(2)
    figure
    
    title('Endurance Track')
    plot(comp.endurance_track(1,:),comp.endurance_track(2,:))

    xlabel('Distance (m)')
    ylabel('Curvature (1/m)')
end

if plot_choice(3)
    figure
    plot(comp.interp_info.radius_vector,comp.interp_info.max_vel_corner_vector)
    xlabel('Radius (m)')
    ylabel('Maximum possible velocity')
end

if plot_choice(4)
    figure
    plot(comp.interp_info.long_vel_guess,comp.interp_info.long_accel_matrix/9.81)
    xlabel('Velocity (m/s)')
    ylabel('Maximum possible longitudinal acceleration (g)')
end

if plot_choice(5)
    %figure
    hold on
    plot(comp.accel.time_vec,comp.accel.long_vel_vector)
    title('Accel')
    xlabel('Time (s)')
    ylabel('Velocity (m/s)')
end

if plot_choice(6)
    figure
    plot(comp.accel.time_vec,comp.accel.long_accel_vector)
    title('Accel')
    xlabel('Time (s)')
    ylabel('Longitudinal acceleration (g)')
end

if plot_choice(7)
    figure
    plot(comp.autocross.time_vec,comp.autocross.long_vel);
    xlabel('Time (s)')
    ylabel('Longitudinal Velocity (m/s)')
    title('Autocross Longitudinal Velocity')
    
    hold on
    for i = 1:numel(comp.car.powertrain.switch_gear_velocities)
        plot(comp.autocross.time_vec,...
            comp.car.powertrain.switch_gear_velocities(i)*ones(size(comp.autocross.time_vec)),'b--');
        gear_shifts(i) = numel(find(diff(sign(comp.autocross.long_vel-comp.car.powertrain.switch_gear_velocities(i)))~=0));
    end
    
    legend(['Number of gear shifts: = ' num2str(gear_shifts)])

end


if plot_choice(8)
    figure
    plot(comp.endurance_track(1,:))
    xlabel('Distance (m)')
    ylabel('Curvature (1/m)')
end
if plot_choice(9)
    figure; scatter(comp.skidpad.ay, comp.skidpad.steer, 'b.'); hold on
    xlabel('Lateral Acceleration (g)'); ylabel('Steer Angle (deg)');

    % Sort
    [ay_sorted, idx] = sort(comp.skidpad.ay(:));
    steer_sorted = comp.skidpad.steer(idx);
    n = numel(ay_sorted);
    % --- Adaptive thresholds ---
    [ay_sorted, idx] = sort(comp.skidpad.ay(:));
    steer_sorted = comp.skidpad.steer(idx);
    n = numel(ay_sorted);
    
    minPts = max(10, ceil(0.25*n));           % more stable seed
    p0   = polyfit(ay_sorted(1:minPts), steer_sorted(1:minPts), 1);
    r0   = steer_sorted(1:minPts) - polyval(p0, ay_sorted(1:minPts));
    sigma = 1.4826*mad(r0,1);
    max_error_deg = max(0.15, 3*sigma);        % residual gate
    
    slope_ref = p0(1);
    slope_tol = max(0.03, 0.20*abs(slope_ref));% <-- key: allow at least 0.03 deg/g drift
    W = 6;                                      % window used for local slope
    
    % --- Grow the linear region ---
    cutoff = n;
    for i = minPts+W:n
        % residuals on the prefix
        p_tmp = polyfit(ay_sorted(1:i), steer_sorted(1:i), 1);
        r     = steer_sorted(1:i) - polyval(p_tmp, ay_sorted(1:i));
    
        % local slope over last W points vs initial slope
        p_win = polyfit(ay_sorted(i-W+1:i), steer_sorted(i-W+1:i), 1);
        slope_dev = abs(p_win(1) - slope_ref);
    
        % stop only when BOTH: residuals exceed noise AND slope has drifted
        if max(abs(r)) > max_error_deg && slope_dev > slope_tol
            cutoff = i - W; 
            break
        end
    end
    
    lin_mask = false(size(ay_sorted)); lin_mask(1:cutoff) = true;
    p = polyfit(ay_sorted(lin_mask), steer_sorted(lin_mask), 1);
    Ku = p(1);

    % Plot fit + cutoff marker
    plot(ay_sorted(lin_mask), polyval(p, ay_sorted(lin_mask)), 'r-', 'LineWidth', 1.5);
    xline(ay_sorted(cutoff), '--k', 'Cutoff','LabelVerticalAlignment','bottom');

    legend('Data', sprintf('Linear fit (K_u = %.3f deg/g)', Ku), 'Location','northwest');
end
