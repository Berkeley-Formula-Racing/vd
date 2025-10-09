function [] = tesla_plots(carCell)

accel_times = zeros(size(carCell, 1), 1);
autocross_times = zeros(size(carCell, 1), 1);
endurance_times = zeros(size(carCell, 1), 1);
peak_torque = zeros(size(carCell, 1), 1);
autocross_vel = [];
endurance_vel = [];

for i = 1:size(carCell, 1)
    accel_times(i, 1) = carCell{i, 1}.comp.times.accel;
    autocross_times(i, 1) = carCell{i, 1}.comp.times.autocross;
    endurance_times(i, 1) = carCell{i, 1}.comp.times.endurance;
    peak_torque(i, 1) = carCell{i, 1}.powertrain.torque_fn(2,2);
    autocross_vel{i} = carCell{i, 1}.comp.autocross.long_vel;
    endurance_vel{i} = carCell{i, 1}.comp.endurance.long_vel;
end
%{
smoothed_accel_times = smooth(accel_times);
smoothed_autocross_times = smooth(autocross_times);
smoothed_endurance_times = smooth(endurance_times);

f_auto = fit(peak_torque, autocross_times, "poly9");
f_accel = fit(peak_torque, accel_times, "poly9");
f_end = fit(peak_torque, endurance_times, "poly9");


auto_objective = @(auto_x_val) feval(f_auto, auto_x_val);

auto_x_min = fminbnd(auto_objective, min(peak_torque), max(peak_torque));
auto_min_y = feval(f_auto, auto_x_min);

figure
plot(f_auto, peak_torque, smoothed_autocross_times)
hold on
line([auto_x_min auto_x_min], [43 52], 'Color', 'r', 'LineStyle', '--');
autoLabelText = sprintf('Optimal Torque: %.2flb-ft; Time: %.2fs', auto_x_min, auto_min_y);
xlabel("Peak Torque (lb-ft)")
ylabel("Autocross Times (s)")
legend("Data", "Fit Line", autoLabelText)
title("Relationship Between Peak Torque and Autocross Event For FSAE Vehicle")
hold off


accel_objective = @(accel_x_val) feval(f_accel, accel_x_val);

accel_x_min = fminbnd(accel_objective, min(peak_torque), max(peak_torque));
accel_min_y = feval(f_accel, accel_x_min);

figure
plot(f_accel, peak_torque, smoothed_accel_times)
hold on
line([110.09 110.09], [2.5 6], 'Color', 'r', 'LineStyle', '--');
accelLabelText = 'Optimal Torque: 110.09lb-ft; Time: 2.95s';
xlabel("Peak Torque (lb-ft)")
ylabel("Acceleration Times (s)")
legend("Data", "Fit Line", accelLabelText)
title("Relationship Between Peak Torque and Acceleration Event For FSAE Vehicle")
hold off


end_objective = @(end_x_val) feval(f_end, end_x_val);

end_x_min = fminbnd(end_objective, min(peak_torque), max(peak_torque));
end_min_y = feval(f_end, end_x_min);

figure
plot(f_end, peak_torque, smoothed_endurance_times)
hold on
line([end_x_min end_x_min], [1250 1600], 'Color', 'r', 'LineStyle', '--');
endLabelText = sprintf('Optimal Torque: %.2flb-ft; Time: %.2fs', end_x_min, end_min_y);
xlabel("Peak Torque (lb-ft)")
ylabel("Endurance Times (s)")
legend("Data", "Fit Line", endLabelText)
title("Relationship Between Peak Torque and Endurance Event For FSAE Vehicle")
hold off
%}

%{
max_rpm_corner_vector = ((carCell{1,1}.comp.interp_info.max_vel_corner_vector) / ((2 * pi * 8) * 0.0254)) * 60 * 7.2918;

useable_auto_corner_velocities = [carCell{1,1}.comp.interp_info.auto_max_apex_velocities(1:77)];
auto_corner_rpm = ((useable_auto_corner_velocities) / ((2 * pi * 8) * 0.0254)) * 60 * 7.2918;

useable_end_corner_velocities = [carCell{1,1}.comp.interp_info.end_max_apex_velocities(1:128) carCell{1,1}.comp.interp_info.end_max_apex_velocities(618:663)];
end_corner_rpm = ((useable_end_corner_velocities) / ((2 * pi * 8) * 0.0254)) * 60 * 7.2918;

bin_edges = [2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000];

figure
histogram(auto_corner_rpm)
hold on
title("Wheel RPM with Final Drive of 7.2918 at Corner Apexes of Autocross Event for FSAE Vehicle")
xlabel("RPM")
ylabel("Count")
hold off

figure
histogram(end_corner_rpm)
hold on
title("Wheel RPM with Final Drive of 7.2918 at Corner Apexes of Endurance Event for FSAE Vehicle")
xlabel("RPM")
ylabel("Count")
hold off

figure
plot(carCell{1,1}.comp.interp_info.radius_vector, max_rpm_corner_vector)
hold on
title("Relationship Between Wheel RPM with Final Drive of 7.2918 and Corner Radii (m)")
xlabel("Corner Radius (m)")
ylabel("RPM")
hold off

figure
histogram(useable_auto_corner_velocities)
hold on
title("Maximum Velocities (m/s) at Corner Apexes of Autocross Event for FSAE Vehicle")
xlabel("Velocity (m/s)")
ylabel("Count")
hold off

figure
histogram(useable_end_corner_velocities)
hold on
title("Maximum Velocities (m/s) at Corner Apexes of Endurance Event for FSAE Vehicle")
xlabel("Velocity (m/s)")
ylabel("Count")
hold off

figure
plot(carCell{1,1}.comp.interp_info.radius_vector, carCell{1,1}.comp.interp_info.max_vel_corner_vector)
hold on
title("Relationship Between Maximum Velocities (m/s) and Increasing Corner Radii (m)")
xlabel("Corner Radius (m)")
ylabel("Velocity (m/s)")
hold off
%}

%{
max_torque_corner_vector = ((carCell{1,1}.M .* carCell{1,1}.comp.interp_info.max_vel_corner_vector.^2) ./ carCell{1,1}.comp.interp_info.radius_vector) .* 0.1956

max_frictional_force = 0.9 * carCell{1,1}.M * 9.81

figure
plot(carCell{1,1}.comp.interp_info.radius_vector, max_torque_corner_vector)
xlabel("Corner Radius")
ylabel("Maximum Predicted Torque")
%}

all_auto_torques = [];
all_auto_velocities = [];
all_end_torques = [];
all_end_velocities = [];

% Flatten the data into two vectors using curly braces for cell array access
for i = 1:length(peak_torque)
    all_auto_torques = [all_auto_torques; repmat(peak_torque(i), length(autocross_vel{i}), 1)];
    all_end_torques = [all_end_torques; repmat(peak_torque(i), length(endurance_vel{i}), 1)];
    all_auto_velocities = [all_auto_velocities; autocross_vel{i}'];
    all_end_velocities = [all_end_velocities; endurance_vel{i}'];
end

% Convert to column vectors (if needed)
all_auto_torques = all_auto_torques(:);
all_auto_velocities = all_auto_velocities(:);
all_end_torques = all_end_torques(:);
all_end_velocities = all_end_velocities(:);

all_auto_rpm = ((all_auto_velocities) / ((2 * pi * 8) * 0.0254)) * 60 * 7.2918;
all_end_rpm = ((all_end_velocities) / ((2 * pi * 8) * 0.0254)) * 60 * 7.2918;

% Define number of bins for torque and velocity
nbins_auto_torque = length(peak_torque); % Number of unique torque values
nbins_auto_velocity = 150; % Adjust as needed for resolution
nbins_auto_rpm = nthroot(length(unique(all_auto_rpm)), 3);
nbins_end_torque = length(peak_torque);
nbins_end_velocity = 150;
nbins_end_rpm = nthroot(length(unique(all_end_rpm)), 3);

% Ensure proper bin edges for torque and velocity
auto_torque_edges = linspace(min(all_auto_torques), max(all_auto_torques), nbins_auto_torque + 1);
auto_velocity_edges = linspace(min(all_auto_velocities), max(all_auto_velocities), nbins_auto_velocity + 1);
auto_rpm_edges = linspace(min(all_auto_rpm), max(all_auto_rpm), nbins_auto_rpm + 1);
end_torque_edges = linspace(min(all_end_torques), max(all_end_torques), nbins_end_torque + 1);
end_velocity_edges = linspace(min(all_end_velocities), max(all_end_velocities), nbins_end_velocity + 1);
end_rpm_edges = linspace(min(all_end_rpm), max(all_end_rpm), nbins_end_rpm + 1);

% Check if the bin edges have at least two elements
if length(auto_torque_edges) < 2 || length(auto_velocity_edges) < 2 || length(auto_rpm_edges) < 2 
    error('Autocross Torque, Velocity, or RPM bin edges must have at least two elements.');
end

if length(end_torque_edges) < 2 || length(end_velocity_edges) < 2|| length(end_rpm_edges) < 2 
    error('Endurance Torque, Velocity, RPM bin edges must have at least two elements.');
end

% Create a 2D histogram (counts) for the heatmap
[auto_vel_counts, ~, ~] = histcounts2(all_auto_torques, all_auto_velocities, auto_torque_edges, auto_velocity_edges);
[auto_rpm_counts, ~, ~] = histcounts2(all_auto_torques, all_auto_rpm, auto_torque_edges, auto_rpm_edges);
[end_vel_counts, ~, ~] = histcounts2(all_end_torques, all_end_velocities, end_torque_edges, end_velocity_edges);
[end_rpm_counts, ~, ~] = histcounts2(all_end_torques, all_end_rpm, end_torque_edges, end_rpm_edges);

%{
% Autocross Torque vs. Velocity
figure;
imagesc(auto_velocity_edges(1:end-1), auto_torque_edges(1:end-1), auto_vel_counts);
colormap('jet'); % Choose a colormap
colorbar; % Add a colorbar
xlabel('Longitudinal Velocity');
ylabel('Engine Torque');
title('Heatmap of Engine Torque vs Longitudinal Velocity for Autocross Event');
axis xy; % Ensure correct orientation

% Endurance Torque vs. Velocity
figure;
imagesc(end_velocity_edges(1:end-1), end_torque_edges(1:end-1), end_vel_counts);
colormap('jet'); % Choose a colormap
colorbar; % Add a colorbar
xlabel('Longitudinal Velocity');
ylabel('Engine Torque');
title('Heatmap of Engine Torque vs Longitudinal Velocity for Endurance Event');
axis xy; % Ensure correct orientation

%}


%{
% Autocross Torque vs. Wheel RPM
figure;
imagesc(auto_rpm_edges, auto_torque_edges, auto_rpm_counts);
colormap('jet'); % Choose a colormap
colorbar; % Add a colorbar
xlabel('RPM');
ylabel('Engine Torque (lb-ft)');
title('Heatmap of Engine Torque vs RPM for FSAE Autocross Event');
axis xy; % Ensure correct orientation

% Autocross Torque vs. Wheel RPM
figure;
imagesc(end_rpm_edges, end_torque_edges, end_rpm_counts);
colormap('jet'); % Choose a colormap
colorbar; % Add a colorbar
xlabel('RPM');
ylabel('Engine Torque (lb-ft)');
title('Heatmap of Engine Torque vs RPM for FSAE Endurance Event');
axis xy; % Ensure correct orientation

%}

auto_s = carCell{1,1}.comp.autocross_track(1, :);
auto_kappa = carCell{1,1}.comp.autocross_track(2, :);

end_s = carCell{1,1}.comp.endurance_track(1, :);
end_kappa = carCell{1,1}.comp.endurance_track(2, :);

auto_theta = cumtrapz(auto_s, auto_kappa);
end_theta = cumtrapz(end_s, end_kappa);

auto_x = cumtrapz(auto_s, cos(auto_theta));
auto_y = cumtrapz(auto_s, sin(auto_theta));

end_x = cumtrapz(end_s, cos(end_theta));
end_y = cumtrapz(end_s, sin(end_theta));

figure
plot(auto_x, auto_y, 'k-', 'LineWidth', 2)
hold on
scatter(auto_x, auto_y, 50, carCell{1,1}.comp.autocross.long_vel, 'filled');
colormap(jet);
colorbar;
title('2024 FSAE Michigan Autocross Track with Velocity (m/s) Heatmap');
axis equal
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
hold off

figure
plot(end_x, end_y, 'k-', 'LineWidth', 2)
hold on
scatter(end_x, end_y, 50, carCell{1,1}.comp.endurance.long_vel, 'filled');
colormap(jet);
colorbar;
title('2024 FSAE Michigan Endurance Track with Velocity (m/s) Heatmap');
axis equal
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
hold off

figure
plot(auto_x, auto_y, 'k-', 'LineWidth', 2)
hold on
scatter(auto_x, auto_y, 50, (carCell{1,1}.comp.autocross.long_vel) / ((2 * pi * 8) * 0.0254) * 60 * 7.2918, 'filled');
colormap(jet);
colorbar;
title('2024 FSAE Michigan Autocross Track with RPM Heatmap');
axis equal
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
hold off

figure
plot(end_x, end_y, 'k-', 'LineWidth', 2)
hold on
scatter(end_x, end_y, 50, (carCell{1,1}.comp.endurance.long_vel) / ((2 * pi * 8) * 0.0254) * 60 * 7.2918, 'filled');
colormap(jet);
colorbar;
title('2024 FSAE Michigan Endurance Track with RPM Heatmap');
axis equal
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
hold off

figure
plot(auto_x, auto_y, 'k-', 'LineWidth', 2)
hold on
scatter(auto_x, auto_y, 50, carCell{1,1}.comp.autocross.long_accel/9.81, 'filled');
colormap(jet);
colorbar;
title('2024 FSAE Michigan Autocross Track with Longitudinal Gs Heatmap');
axis equal
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
hold off

figure
plot(end_x, end_y, 'k-', 'LineWidth', 2)
hold on
scatter(end_x, end_y, 50, carCell{1,1}.comp.endurance.long_accel/9.81, 'filled');
colormap(jet);
colorbar;
title('2024 FSAE Michigan Endurance Track with Longitudinal Gs Heatmap');
axis equal
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
hold off

figure
plot(auto_x, auto_y, 'k-', 'LineWidth', 2)
hold on
scatter(auto_x, auto_y, 50, carCell{1,1}.comp.autocross.lat_accel/9.81, 'filled');
colormap(jet);
colorbar;
title('2024 FSAE Michigan Autocross Track with Lateral Gs Heatmap');
axis equal
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
hold off

figure
plot(end_x, end_y, 'k-', 'LineWidth', 2)
hold on
scatter(end_x, end_y, 50, carCell{1,1}.comp.endurance.lat_accel/9.81, 'filled');
colormap(jet);
colorbar;
title('2024 FSAE Michigan Endurance Track with Lateral Gs Heatmap');
axis equal
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
hold off