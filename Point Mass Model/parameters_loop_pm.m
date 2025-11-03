function carCell = parameters_loop_pm(carParams, aeroParams, eParams, DTparams, Bparams, muParams)
    trk = track_pm.loadFromMichiganMat('michigantrack2024.mat');
    sim = lapsim_pm();
    vmax = 33;

    % get torque source once
    [rpm_vec_src, tq_lbft_vec_src] = getTorqueSource(eParams.engineSource);

    % lengths
    L.mass      = numel(carParams.mass);
    L.drv       = numel(carParams.driver_weight);
    L.wd        = numel(carParams.weight_dist);
    L.wr        = numel(carParams.wheel_radius);
    L.croll     = numel(carParams.c_roll);

    L.cda       = numel(aeroParams.cda);
    L.cla       = numel(aeroParams.cla);
    L.cop       = numel(aeroParams.distribution);

    L.fd        = numel(DTparams.final_drive);
    L.eta       = numel(DTparams.drivetrain_efficiency);

    L.brktq     = numel(Bparams.max_braking_torque);
    L.brkTL     = numel(Bparams.tire_limited_brake);

    L.mu_scale  = numel(muParams.track_friction_scale);
    L.mu_coeff  = numel(muParams.mu_model_coeff);
    L.mu_exp    = numel(muParams.mu_model_exp);

    [i_mass, i_drv, i_wd, i_wr, i_croll, ...
     i_cda, i_cla, i_cop, ...
     i_fd, i_eta, ...
     i_brktq, i_brkTL, ...
     i_mu_scale, i_mu_coeff, i_mu_exp] = ndgrid( ...
        1:L.mass, 1:L.drv, 1:L.wd, 1:L.wr, 1:L.croll, ...
        1:L.cda, 1:L.cla, 1:L.cop, ...
        1:L.fd, 1:L.eta, ...
        1:L.brktq, 1:L.brkTL, ...
        1:L.mu_scale, 1:L.mu_coeff, 1:L.mu_exp);

    N = numel(i_mass);
    carCell = cell(N,1);

    for n = 1:N
        spec = struct();

        % car
        spec.m            = carParams.mass(i_mass(n)) + carParams.driver_weight(i_drv(n));
        spec.weight_r     = carParams.weight_dist(i_wd(n));
        spec.wheel_radius = carParams.wheel_radius(i_wr(n));
        spec.c_roll       = carParams.c_roll(i_croll(n));

        % aero
        spec.cda          = aeroParams.cda(i_cda(n));
        spec.cla          = aeroParams.cla(i_cla(n));
        spec.cop          = aeroParams.distribution(i_cop(n));

        % env
        spec.rho = 1.2;
        spec.g   = 9.81;

        % μ
        spec.track_friction_scale = muParams.track_friction_scale(i_mu_scale(n));
        spec.mu_model_coeff       = muParams.mu_model_coeff(i_mu_coeff(n));
        spec.mu_model_exp         = muParams.mu_model_exp(i_mu_exp(n));

        % powertrain
        spec.gears         = eParams.gears;  % fixed
        spec.final_drive   = DTparams.final_drive(i_fd(n));
        spec.eta_driveline = DTparams.drivetrain_efficiency(i_eta(n));
        spec.rpm_redline   = eParams.redline;
        spec.rpm_idle      = eParams.rpm_idle;
        spec.rpm_vec       = rpm_vec_src;
        spec.tq_lbft_vec   = tq_lbft_vec_src;

        % brakes
        spec.max_brake_torque   = Bparams.max_braking_torque(i_brktq(n));
        spec.tire_limited_brake = Bparams.tire_limited_brake(i_brkTL(n));

        % build + run
        carObj = buildCarFromSpec_pm(spec);
        [lapT, vprof, vlim] = sim.run(carObj, trk, vmax);

        out = struct();
        out.carObj    = carObj;
        out.spec      = spec;
        out.lap_time  = lapT;
        out.v_profile = vprof;
        out.v_latlim  = vlim;

        carCell{n} = out;

        fprintf('combo %d/%d: lap time = %.2f s\n', n, N, lapT);
    end
end

function [rpm_vec, tq_lbft_vec] = getTorqueSource(engineSource)
    switch string(engineSource)
        case "KTM450"
            KTM = KTM450();
            rpm_vec     = KTM(1,:).';
            tq_lbft_vec = KTM(2,:).';
        case "flat30"
            rpm_vec     = (3000:500:11000).';
            tq_lbft_vec = 30*ones(size(rpm_vec));
        otherwise
            error("Unknown engineSource '%s'", engineSource);
    end
end
