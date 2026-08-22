%% Master RCVD Bicycle Model: Target Wheelbase Optimization Suite
clear; clc;
setup_paths;

%% 1. Initialize Baseline Vehicle & Lock Static Parameters
disp('Step 1: Loading baseline car and locking static weight distribution...');
[carCell,eventParams] = carConfig();
car = carCell{1,1};       
accelCar = carCell{1,2};  
numWorkers = 16; 

g = car.g;
track_width = car.t_f;    
rho = 1.162;              
ClA_base = car.aero.cla;       
aero_dist_f = 0.539; 
aero_dist_r = 1 - aero_dist_f;

% --- STATIC WEIGHT DISTRIBUTION SETUP (Locked as Requested) ---
M = car.M;                        % Total vehicle mass with driver from scales
weight_dist_r = 0.512;            % Locked at 51.2% rear distribution
weight_dist_f = 1 - weight_dist_r;
total_weight = M * g;

% --- UNDERTRAY DOWNFORCE SCALING PARAMETERS ---
ut_fraction = 230 / 830;          % Undertray provides 27.7% of baseline downforce
ut_sensitivity = 1.0;             % 1.0 = Linear loss with length; >1.0 = Aggressive loss; <1.0 = Conservative
baseline_L_in = 62.0;             % Baseline chassis length in inches

fprintf('----------------------------------------------------\n');
fprintf('  STATIC CONFIGURATION LOCKED\n');
fprintf('----------------------------------------------------\n');
fprintf('Total Vehicle Mass        : %6.2f kg (%5.1f lbs)\n', M, M * 2.20462);
fprintf('Static Front Distribution : %6.2f %%\n', weight_dist_f * 100);
fprintf('Static Rear Distribution  : %6.2f %%\n', weight_dist_r * 100);
fprintf('----------------------------------------------------\n\n');

%% 2. Run Run-Once Track Solver & Extract 3 Low / 3 High Speed Corners
disp('Step 2: Processing track solver to cluster velocity conditions...');
paramsArr = gg2(car, numWorkers); 
car = makeGG(paramsArr, car);
comp = Events2(car, accelCar, eventParams);
comp.calcTimes(); 

load('michigantrack2024.mat'); 
V_profile  = comp.autocross.long_vel;
ay_profile = comp.autocross.lat_accel;

% Extract 3 Low-Speed Hairpins & 3 High-Speed Sweepers
[~, locs_low] = findpeaks(abs(curvature), 'SortStr', 'descend', 'MinPeakDistance', 50);
idx_low_3 = locs_low(1:3);
R_low = 1 ./ abs(curvature(idx_low_3));
V_low = V_profile(idx_low_3);
ay_low = abs(ay_profile(idx_low_3));

sweep_candidates = find(abs(curvature) > 0.015 & abs(curvature) < 0.035);
[~, sweep_locs] = findpeaks(V_profile(sweep_candidates), 'SortStr', 'descend', 'MinPeakDistance', 30);
idx_high_3 = sweep_candidates(sweep_locs(1:3));
R_high = 1 ./ abs(curvature(idx_high_3));
V_high = V_profile(idx_high_3);
ay_high = abs(ay_profile(idx_high_3));

%% 3. Dynamic Empirical Tire Load Sensitivity (Pacejka Sourced)
disp('Step 3: Calculating empirical Pacejka load sensitivity slopes...');
Fz_f_stat = (total_weight * weight_dist_f) / 2;
Fz_r_stat = (total_weight * weight_dist_r) / 2;
alpha_deg = 0.1; 
alpha_rad = alpha_deg * (pi / 180);

Fy_f1 = car.tire.F_y(alpha_deg, 0, Fz_f_stat, 0);
Fy_f2 = car.tire.F_y(alpha_deg, 0, Fz_f_stat + 500, 0);
Fy_r1 = car.tire.F_y(alpha_deg, 0, Fz_r_stat, 0);
Fy_r2 = car.tire.F_y(alpha_deg, 0, Fz_r_stat + 500, 0);

Cf_base = abs(2 * (Fy_f1 / alpha_rad));
Cr_base = abs(2 * (Fy_r1 / alpha_rad));
Cf_loaded = abs(2 * (Fy_f2 / alpha_rad));
Cr_loaded = abs(2 * (Fy_r2 / alpha_rad));

dCf_dFz = (Cf_base - Cf_loaded) / 1000; 
dCr_dFz = (Cr_base - Cr_loaded) / 1000;
fprintf('Empirical Load Sensitivity: dCf/dFz = %.4f | dCr/dFz = %.4f\n\n', dCf_dFz, dCr_dFz);

%% 4. Wheelbase Parameter Sweep Loop across Multi-Corner Cluster
disp('Step 4: Sweeping geometric parameters and scaling undertray downforce...');
wheelbase_sweep_in = linspace(60, 62.5, 6); 
numSteps = length(wheelbase_sweep_in);

steer_low  = zeros(numSteps, 3);
steer_high = zeros(numSteps, 3);
SI_low     = zeros(numSteps, 3);
SI_high    = zeros(numSteps, 3);
vel_sweep_mps = linspace(5, 35, 10); 
sens_matrix   = zeros(numSteps, length(vel_sweep_mps));
dynamic_ClA_log = zeros(1, numSteps);

% Pre-allocate arrays for target optimization metrics
launch_grip_gain = zeros(1, numSteps);
low_speed_agility = zeros(1, numSteps);
high_speed_aero = zeros(1, numSteps);
high_speed_stability = zeros(1, numSteps);

for i = 1:numSteps
    L = wheelbase_sweep_in(i) * 0.0254; % meters
    l_f = L * weight_dist_r;            
    l_r = L * weight_dist_f;            
    
    % --- DYNAMIC UNDERTRAY DOWNFORCE DEGRADATION ---
    length_ratio = wheelbase_sweep_in(i) / baseline_L_in;
    dynamic_ClA = ClA_base * ((1 - ut_fraction) + ut_fraction * (length_ratio ^ ut_sensitivity));
    dynamic_ClA_log(i) = dynamic_ClA;
    car.aero.cla = dynamic_ClA;
    % -----------------------------------------------
    
    % --- Process 3 Low-Speed Corners ---
    for c = 1:3
        Fy_f = (M * ay_low(c) * l_r) / L;
        Fy_r = (M * ay_low(c) * l_f) / L;
        Fz_f = ((M * g * l_r) / L) + (0.5 * rho * (V_low(c)^2) * dynamic_ClA * aero_dist_f);
        Fz_r = ((M * g * l_f) / L) + (0.5 * rho * (V_low(c)^2) * dynamic_ClA * aero_dist_r);
        
        Cf = Cf_base - (dCf_dFz * (Fz_f - (total_weight*weight_dist_f)));
        Cr = Cr_base - (dCr_dFz * (Fz_r - (total_weight*weight_dist_r)));
        
        balance = (Fy_f / Cf) - (Fy_r / Cr);
        steer_low(i, c) = rad2deg(L / R_low(c)) + rad2deg(balance);
        
        SM = (l_f * Cf - l_r * Cr) / (L * (Cf + Cr));
        C0 = (Cf * Cr * L^2) / (2 * M * (Cf + Cr));
        SI_low(i, c) = -SM + (C0 * (R_low(c) * g) / (V_low(c)^2)) / 1000; 
    end
    
    % --- Process 3 High-Speed Corners ---
    for c = 1:3
        Fy_f = (M * ay_high(c) * l_r) / L;
        Fy_r = (M * ay_high(c) * l_f) / L;
        Fz_f = ((M * g * l_r) / L) + (0.5 * rho * (V_high(c)^2) * dynamic_ClA * aero_dist_f);
        Fz_r = ((M * g * l_f) / L) + (0.5 * rho * (V_high(c)^2) * dynamic_ClA * aero_dist_r);
        
        Cf = Cf_base - (dCf_dFz * (Fz_f - (total_weight*weight_dist_f)));
        Cr = Cr_base - (dCr_dFz * (Fz_r - (total_weight*weight_dist_r)));
        
        balance = (Fy_f / Cf) - (Fy_r / Cr);
        steer_high(i, c) = rad2deg(L / R_high(c)) + rad2deg(balance);
        
        SM = (l_f * Cf - l_r * Cr) / (L * (Cf + Cr));
        C0 = (Cf * Cr * L^2) / (2 * M * (Cf + Cr));
        SI_high(i, c) = -SM + (C0 * (R_high(c) * g) / (V_high(c)^2)) / 1000;
    end
    
    % --- Calculate RCVD Fig 5.52 Steering Sensitivity across 10 Speeds in m/s ---
    for v = 1:length(vel_sweep_mps)
        V_val = vel_sweep_mps(v);
        Fz_f_val = ((M * g * l_r) / L) + (0.5 * rho * (V_val^2) * dynamic_ClA * aero_dist_f);
        Fz_r_val = ((M * g * l_f) / L) + (0.5 * rho * (V_val^2) * dynamic_ClA * aero_dist_r);
        
        Cf_val = Cf_base - (dCf_dFz * (Fz_f_val - (total_weight*weight_dist_f)));
        Cr_val = Cr_base - (dCr_dFz * (Fz_r_val - (total_weight*weight_dist_r)));
        
        UG_deg_per_g = rad2deg((total_weight * weight_dist_f) / Cf_val - (total_weight * weight_dist_r) / Cr_val);
        ack_grad_deg_per_g = rad2deg(L * g / (V_val^2));
        sens_matrix(i, v) = 1 / (UG_deg_per_g + ack_grad_deg_per_g);
    end
    
    % --- EXTRACT TARGET OPTIMIZATION METRICS ---
    % 1. Launch Traction (Rear load under 1.2G launch)
    launch_grip_gain(i) = (total_weight * weight_dist_r) + ((M * (1.2 * g) * car.h_g) / L);
    % 2. Low-Speed Agility (Average total steering angle in hairpins)
    low_speed_agility(i) = mean(steer_low(i, :));
    % 3. High-Speed Aero Grip (Dynamic ClA)
    high_speed_aero(i) = dynamic_ClA;
    % 4. High-Speed Drivability Risk (Average positive SI in sweepers)
    high_speed_stability(i) = mean(SI_high(i, :));
end

%% 5. Plot RCVD Multi-Corner Diagnostics across Separate Tabs
disp('Step 5: Rendering dark-mode diagnostic tabs with sloping Ackermann baselines...');
colors = ['b', 'c', 'm']; 

% TAB 1: LOW-SPEED STEERING DEMAND
figure('Name', 'Low-Speed Steering', 'Position', [100, 100, 750, 500]);
hold on;
for c = 1:3
    plot(wheelbase_sweep_in, steer_low(:, c), [colors(c) '-o'], 'LineWidth', 2.5, ...
        'DisplayName', sprintf('Corner %d (R=%.1fm, V=%.1fm/s)', c, R_low(c), V_low(c)));
    ack_base_low = rad2deg((wheelbase_sweep_in * 0.0254) / R_low(c));
    plot(wheelbase_sweep_in, ack_base_low, [colors(c) '--'], 'LineWidth', 1.5, 'HandleVisibility', 'off');
end
grid on; box on;
xlim([wheelbase_sweep_in(1) wheelbase_sweep_in(end)]);
xlabel('Wheelbase Length (inches)', 'FontWeight', 'bold', 'FontSize', 11);
ylabel('Total Steering Angle (degrees)', 'FontWeight', 'bold', 'FontSize', 11);
title('LOW-SPEED CLUSTER: Steering Demand vs. Sloping Ackermann Geometry', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best');

% TAB 2: HIGH-SPEED STEERING DEMAND
figure('Name', 'High-Speed Steering', 'Position', [150, 150, 750, 500]);
hold on;
for c = 1:3
    plot(wheelbase_sweep_in, steer_high(:, c), [colors(c) '-o'], 'LineWidth', 2.5, ...
        'DisplayName', sprintf('Sweeper %d (R=%.1fm, V=%.1fm/s)', c, R_high(c), V_high(c)));
    ack_base_high = rad2deg((wheelbase_sweep_in * 0.0254) / R_high(c));
    plot(wheelbase_sweep_in, ack_base_high, [colors(c) '--'], 'LineWidth', 1.5, 'HandleVisibility', 'off');
end
grid on; box on;
xlim([wheelbase_sweep_in(1) wheelbase_sweep_in(end)]);
xlabel('Wheelbase Length (inches)', 'FontWeight', 'bold', 'FontSize', 11);
ylabel('Total Steering Angle (degrees)', 'FontWeight', 'bold', 'FontSize', 11);
title('HIGH-SPEED CLUSTER: Steering Demand vs. Sloping Ackermann Geometry', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best');

% TAB 3: RCVD STABILITY INDEX (SI)
figure('Name', 'RCVD Stability Index', 'Position', [200, 200, 750, 500]);
hold on;
for c = 1:3
    plot(wheelbase_sweep_in, SI_low(:, c), [colors(c) '-o'], 'LineWidth', 2, ...
        'DisplayName', sprintf('Low-Speed %d (R=%.1fm)', c, R_low(c)));
    plot(wheelbase_sweep_in, SI_high(:, c), [colors(c) '--s'], 'LineWidth', 2, ...
        'DisplayName', sprintf('High-Speed %d (R=%.1fm)', c, R_high(c)));
end
yline(0, 'g--', 'LineWidth', 2, 'DisplayName', 'Neutral Steer / Critical Speed Line');
grid on; box on;
xlim([wheelbase_sweep_in(1) wheelbase_sweep_in(end)]);
xlabel('Wheelbase Length (inches)', 'FontWeight', 'bold', 'FontSize', 11);
ylabel('Stability Index SI = \partial C_N / \partial A_y', 'FontWeight', 'bold', 'FontSize', 11);
title('RCVD Eq 5.68: Directional Stability Index vs Wheelbase', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9);
yl3 = ylim; ymid3 = mean(yl3);
text(wheelbase_sweep_in(2), ymid3 + (yl3(2)-ymid3)*0.4, '\leftarrow Unstable Directional Zone (+)', 'FontWeight', 'bold', 'Color', [0.7 0.2 0.2]);
text(wheelbase_sweep_in(2), ymid3 - (ymid3-yl3(1))*0.4, '\leftarrow Stable Weathercock Zone (-)', 'FontWeight', 'bold', 'Color', [0.2 0.5 0.2]);

% TAB 4: RCVD FIG 5.52 STEERING SENSITIVITY
figure('Name', 'Steering Sensitivity', 'Position', [250, 250, 750, 500]);
hold on;
colors_sens = cool(numSteps); 
for i = 1:numSteps
    if abs(wheelbase_sweep_in(i) - 62) < 0.1
        plot(vel_sweep_mps, sens_matrix(i, :), '--o', 'Color', colors_sens(i, :), 'LineWidth', 3, ...
            'MarkerSize', 6, 'MarkerFaceColor', colors_sens(i, :), ...
            'DisplayName', sprintf('%.1f in Wheelbase (Current Baseline)', wheelbase_sweep_in(i)));
    else
        plot(vel_sweep_mps, sens_matrix(i, :), '-o', 'Color', colors_sens(i, :), 'LineWidth', 2, ...
            'MarkerSize', 5, 'MarkerFaceColor', colors_sens(i, :), ...
            'DisplayName', sprintf('%.1f in Wheelbase', wheelbase_sweep_in(i)));
    end
end
grid on; box on;
xlim([vel_sweep_mps(1) vel_sweep_mps(end)]);
xlabel('Velocity (m/s)', 'FontWeight', 'bold', 'FontSize', 11);
ylabel('Steering Sensitivity dA_y / d\delta (g per deg)', 'FontWeight', 'bold', 'FontSize', 11);
title('RCVD Fig 5.52: Steering Sensitivity across 10 Speeds & All 6 Wheelbases', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'northwest', 'FontSize', 9);
yl4 = ylim;
text(vel_sweep_mps(4), yl4(2)*0.85, 'Shorter wheelbases sit on higher, steeper sensitivity curves \rightarrow', ...
    'FontWeight', 'bold', 'Color', colors_sens(1, :), 'HorizontalAlignment', 'left');

%% 6. Target Wheelbase Optimization & Trade-Study Matrix
disp('Step 6: Executing objective trade-study optimization to find target wheelbase...');

% Normalize all 4 metrics from 0 to 1 (1 = Best Performance/Lowest Risk)
norm_launch = (launch_grip_gain - min(launch_grip_gain)) / (max(launch_grip_gain) - min(launch_grip_gain));
norm_agility = (max(low_speed_agility) - low_speed_agility) / (max(low_speed_agility) - min(low_speed_agility));
norm_aero = (high_speed_aero - min(high_speed_aero)) / (max(high_speed_aero) - min(high_speed_aero));
norm_stability = (max(high_speed_stability) - high_speed_stability) / (max(high_speed_stability) - min(high_speed_stability));

% Apply FSAE Dynamic Event Weightings (Total = 100 Points)
% 30% High-Speed Aero Grip (AutoX / Endurance cornering speeds)
% 25% High-Speed Stability (Endurance driver consistency / spin prevention)
% 25% Launch Traction (75m Acceleration event sprint)
% 20% Low-Speed Agility (Slalom and hairpin steering reduction)
total_score = (30 * norm_aero) + (25 * norm_stability) + (25 * norm_launch) + (20 * norm_agility);

% Identify Optimal Target Wheelbase
[max_score, target_idx] = max(total_score);
L_target = wheelbase_sweep_in(target_idx);

% Print Decision Matrix to Command Window
fprintf('\n=========================================================================================\n');
fprintf('                 FSAE OBJECTIVE WHEELBASE TRADE-STUDY MATRIX                             \n');
fprintf('=========================================================================================\n');
fprintf(' Wheelbase | Launch Grip | Hairpin Steer | Undertray ClA | Stability Risk |  TOTAL SCORE  \n');
fprintf('   (in)    |   Gain (N)  |  Demand (deg) |   (Aero Grip) |  (SI at 45mph) |   (out of 100)\n');
fprintf('-----------------------------------------------------------------------------------------\n');
for i = 1:numSteps
    if i == target_idx
        marker = '<--- OPTIMAL TARGET';
    else
        marker = '';
    end
    fprintf('  %5.1f    |   %6.1f    |     %5.2f     |     %5.3f     |    +%5.4f     |    %5.1f  %s\n', ...
        wheelbase_sweep_in(i), launch_grip_gain(i), low_speed_agility(i), ...
        high_speed_aero(i), high_speed_stability(i), total_score(i), marker);
end
fprintf('=========================================================================================\n\n');

% TAB 5: OBJECTIVE TARGET WHEELBASE OPTIMIZATION CURVE
figure('Name', 'Target Wheelbase Optimizer', 'Position', [300, 300, 750, 450]);
plot(wheelbase_sweep_in, total_score, 'g-o', 'LineWidth', 3, 'MarkerSize', 8, ...
    'MarkerFaceColor', 'g', 'DisplayName', 'Composite FSAE Performance Score');
hold on;
plot(L_target, max_score, 'wp', 'MarkerSize', 16, 'MarkerFaceColor', 'm', 'LineWidth', 2, ...
    'DisplayName', sprintf('Optimal Target: %.1f inches', L_target));
grid on; box on;
xlim([wheelbase_sweep_in(1) wheelbase_sweep_in(end)]);
xlabel('Wheelbase Length (inches)', 'FontWeight', 'bold', 'FontSize', 11);
ylabel('FSAE Performance Score (out of 100)', 'FontWeight', 'bold', 'FontSize', 11);
title('OBJECTIVE TARGET NUMBER: Composite Trade-Study Optimization', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 10);
yl5 = ylim; ymid5 = mean(yl5);
text(L_target, max_score - (yl5(2)-yl5(1))*0.15, sprintf('\\leftarrow Recommended Target: %.1f"', L_target), ...
    'FontWeight', 'bold', 'Color', 'w', 'FontSize', 11);

disp('All 5 diagnostic and optimization tabs rendered successfully.');