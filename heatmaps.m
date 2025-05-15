function [] = heatmaps(carCell, n)

autocross_vel = [];
endurance_vel = [];

for i = 1:size(carCell, 1)
    autocross_vel{i} = carCell{i, 1}.comp.autocross.long_vel;
    endurance_vel{i} = carCell{i, 1}.comp.endurance.long_vel;
end

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


auto_d = [0, cumsum(sqrt(diff(auto_x).^2 + diff(auto_y).^2))];
auto_dq = linspace(auto_d(1), auto_d(end), 100000);

auto_xq = interp1(auto_d, auto_x, auto_dq, 'spline');
auto_yq = interp1(auto_d, auto_y, auto_dq, 'spline');
auto_vq = interp1(auto_d, autocross_vel{1,1}, auto_dq, 'spline');


end_d = [0, cumsum(sqrt(diff(end_x).^2 + diff(end_y).^2))];
end_dq = linspace(end_d(1), end_d(end), 1000000);

end_xq = interp1(end_d, end_x, end_dq, 'spline');
end_yq = interp1(end_d, end_y, end_dq, 'spline');
end_vq = interp1(end_d, endurance_vel{1,1}, end_dq, 'spline');



figure
plot(auto_x, auto_y, 'k-', 'LineWidth', 2)
hold on
scatter(auto_xq, auto_yq, 10, auto_vq, 'filled');
colormap('jet');
colorbar;
title('2025 FSAE Michigan Autocross Track with Velocity (m/s) Heatmap');
axis equal
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
hold off

figure
plot(-end_yq, end_xq, 'k-', 'LineWidth', 2)
hold on
scatter(-end_yq, end_xq, 10, end_vq, 'filled');
colormap('jet');
colorbar;
title('2025 FSAE Michigan Endurance Track with Velocity (m/s) Heatmap');
axis equal
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
hold off




auto_s2 = carCell{2,1}.comp.autocross_track(1, :);
auto_kappa2 = carCell{2,1}.comp.autocross_track(2, :);

end_s2 = carCell{2,1}.comp.endurance_track(1, :);
end_kappa2 = carCell{2,1}.comp.endurance_track(2, :);

auto_theta2 = cumtrapz(auto_s2, auto_kappa2);
end_theta2 = cumtrapz(end_s2, end_kappa2);

auto_x2 = cumtrapz(auto_s2, cos(auto_theta2));
auto_y2 = cumtrapz(auto_s2, sin(auto_theta2));

end_x2 = cumtrapz(end_s2, cos(end_theta2));
end_y2 = cumtrapz(end_s2, sin(end_theta2));


auto_d2 = [0, cumsum(sqrt(diff(auto_x2).^2 + diff(auto_y2).^2))];
auto_dq2 = linspace(auto_d2(1), auto_d2(end), 100000);

auto_xq2 = interp1(auto_d2, auto_x2, auto_dq2, 'spline');
auto_yq2 = interp1(auto_d2, auto_y2, auto_dq2, 'spline');
auto_vq2 = interp1(auto_d2, autocross_vel{1,2}, auto_dq2, 'spline');


end_d2 = [0, cumsum(sqrt(diff(end_x2).^2 + diff(end_y2).^2))];
end_dq2 = linspace(end_d2(1), end_d2(end), 1000000);

end_xq2 = interp1(end_d2, end_x2, end_dq2, 'spline');
end_yq2 = interp1(end_d2, end_y2, end_dq2, 'spline');
end_vq2 = interp1(end_d2, endurance_vel{1,2}, end_dq2, 'spline');

figure
plot(auto_x, auto_y, 'k-', 'LineWidth', 2)
hold on
scatter(auto_xq, auto_yq, 10, auto_vq - auto_vq2, 'filled');
colormap('jet');
colorbar;
title('vel diff');
axis equal
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
hold off

figure
plot(end_x, end_y, 'k-', 'LineWidth', 2)
hold on
scatter(end_xq, end_yq, 10, end_vq - end_vq2, 'filled');
colormap('jet');
colorbar;
title('vel diff');
axis equal
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
hold off