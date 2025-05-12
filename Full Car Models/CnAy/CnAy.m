clear
setup_paths
carCell = carConfig(); % generate all cars to sim over
numCars = size(carCell, 1);

% Initialize results struct
results(numCars) = struct();

steer_angles = linspace(-20,20, 41);
betas = [-8 -4 -2 -1 0 1 2 4 8];

velocity = 20; % m/s

for i = 1:numCars
    car = carCell{i,1};
    colors = ['r' 'k' 'b' 'g' 'm'];
    figure

    % Initialize storage variables
    max_lat_accel = -inf;
    max_Cn = -inf;
    Cn_at_max_Ay = 0;
    Ay_at_max_Cn = 0;
    max_Ay_at_Cn0 = -inf;

    steer_angle_store = [];
    beta_store = [];
    Cn_store = [];
    Ay_store = [];
    SAf_store = [];
    SAr_store = [];
    Delta_store = [];
    Fyf_store = [];
    Fyr_store = [];

    all_delta = [];
    all_Cn = [];
    all_beta = [];
    all_Cn_beta = [];

    %% Sweep over steer angles
    for steer_angle = [-16 -8 -4 -2 0 2 4 8 16]
        [lat_accel, yaw_accel, lat_velocities, x_table_all] = generate_constant_steer_line(car, velocity, steer_angle, [0,0]);
        Cn = yaw_accel * car.I_zz;
        beta_vals = (atan(lat_velocities / velocity)) * 180 / pi;

        lat_accel(32) = []; lat_accel(31) = []; Cn(32) = []; Cn(31) = [];
        x_table_all(32) = []; x_table_all(31) = [];

        plot(lat_accel / 9.81, Cn, colors(i))
        hold on

        % Store for slope estimation
        all_delta = [all_delta, steer_angle * ones(1, length(Cn))];
        all_Cn = [all_Cn, Cn];

        [Aymax, idx] = max(lat_accel);
        if Aymax > max_lat_accel
            max_lat_accel = Aymax;
            Cn_at_max_Ay = Cn(idx);

            % Interpolate all quantities at Ay = max_lat_accel
            SAf_interp = NaN(size(x_table_all));
            SAr_interp = NaN(size(x_table_all));
            Fyf_interp = NaN(size(x_table_all));
            Fyr_interp = NaN(size(x_table_all));
            beta_interp = NaN(size(x_table_all));

            for j = 1:length(x_table_all)
                x_row = x_table_all{j};
                alpha = table2array(x_row(1, contains(x_row.Properties.VariableNames, 'alpha')));
                T = table2array(x_row(1, contains(x_row.Properties.VariableNames, 'T')));
                SAf_interp(j) = mean(alpha(1:2));
                SAr_interp(j) = mean(alpha(3:4));
                Fyf_interp(j) = mean(T(1:2));
                Fyr_interp(j) = mean(T(3:4));
                beta_interp(j) = beta_vals(j);
            end

            SAf = interp1(lat_accel, SAf_interp, max_lat_accel, 'linear', 'extrap');
            SAr = interp1(lat_accel, SAr_interp, max_lat_accel, 'linear', 'extrap');
            Fyf = interp1(lat_accel, Fyf_interp, max_lat_accel, 'linear', 'extrap');
            Fyr = interp1(lat_accel, Fyr_interp, max_lat_accel, 'linear', 'extrap');
            Beta = interp1(lat_accel, beta_interp, max_lat_accel, 'linear', 'extrap');

            SAf_store(i) = SAf;
            SAr_store(i) = SAr;
            Fyf_store(i) = Fyf;
            Fyr_store(i) = Fyr;

            Delta_store(i) = steer_angle;
            beta_store(i) = Beta;
        end

        [Cnmax, Cnidx] = max(abs(Cn));
        if Cn(Cnidx) > max_Cn
            max_Cn = Cn(Cnidx);
            Ay_at_max_Cn = lat_accel(Cnidx);
        end

        idx_near_zero_Cn = find(abs(Cn) < 40);
        if ~isempty(idx_near_zero_Cn)
            max_Ay_near_Cn0 = max(lat_accel(idx_near_zero_Cn));
            if max_Ay_near_Cn0 > max_Ay_at_Cn0
                max_Ay_at_Cn0 = max_Ay_near_Cn0;
            end
        end
    end

    %% Sweep over beta values
    for beta = betas
        [lat_accel, yaw_accel] = generate_constant_beta_line(car, velocity, beta);
        Cn = yaw_accel * car.I_zz;
        plot(lat_accel / 9.81, Cn, colors(i))
        hold on
        all_beta = [all_beta, beta * ones(1, length(Cn))];
        all_Cn_beta = [all_Cn_beta, Cn];
    end

    %% Revert to basic finite difference slope calculation at β=0 and δ=0
    idx_dpos = find(steer_angles == 2);
    idx_dneg = find(steer_angles == -2);
    if ~isempty(idx_dpos) && ~isempty(idx_dneg)
        dN_dDelta = (all_Cn(idx_dpos) - all_Cn(idx_dneg)) / deg2rad(4);
    else
        dN_dDelta = NaN;
    end

    idx_bpos = find(betas == 2);
    idx_bneg = find(betas == -2);
    if ~isempty(idx_bpos) && ~isempty(idx_bneg)
        dN_dBeta = (all_Cn_beta(idx_bpos) - all_Cn_beta(idx_bneg)) / deg2rad(4);
    else
        dN_dBeta = NaN;
    end

    %% Store results
    results(i).maxAy = max_lat_accel;
    results(i).Cn_at_maxAy = Cn_at_max_Ay;
    results(i).SAf = SAf_store(i);
    results(i).SAr = SAr_store(i);
    results(i).Beta = beta_store(i);
    results(i).Delta = Delta_store(i);
    results(i).Fyf = Fyf_store(i);
    results(i).Fyr = Fyr_store(i);
    results(i).maxCn = max_Cn;
    results(i).Ay_at_maxCn = Ay_at_max_Cn;
    results(i).maxAy_at_Cn0 = max_Ay_at_Cn0;
    results(i).dN_dDelta = dN_dDelta;
    results(i).dN_dBeta = dN_dBeta;

    xlabel('Ay (g)')
    ylabel('Cn')
    grid on
    title(['Car ' num2str(i)])
end

% Display table of results
T = struct2table(results);
disp(T);