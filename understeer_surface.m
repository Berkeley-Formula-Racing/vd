function [] = understeer_surface(carCell)
    % Use first column if passed as an N×M cell matrix
    carCell = carCell(:,1);

    num_cars = numel(carCell);

    % Preallocations
    Ku   = nan(1, num_cars);   % understeer gradient (deg/g)
    Cf   = nan(1, num_cars);   % camber compliance front
    Cr   = nan(1, num_cars);   % camber compliance rear
    Mass = nan(1, num_cars);   % mass if available

    for i = 1:num_cars
        car = carCell{i};
        comp  = car.comp;
        ay    = comp.skidpad.ay(:);
        steer = comp.skidpad.steer(:);

        % --- sort ---
        [ay, k]  = sort(ay);
        steer    = steer(k);
        n        = numel(ay);
        if n < 8, warning('Car %d: not enough points. Skipping.', i); continue; end

        % --- adaptive thresholds from low-g region ---
        minPts = max(10, ceil(0.25*n));                 % seed length
        p0   = polyfit(ay(1:minPts), steer(1:minPts), 1);
        r0   = steer(1:minPts) - polyval(p0, ay(1:minPts));
        sigma = 1.4826*mad(r0,1);                       % robust std
        max_error_deg = max(0.15, 3*sigma);             % residual gate

        slope_ref = p0(1);
        slope_tol = max(0.03, 0.20*abs(slope_ref));     % allow >=0.03 deg/g drift
        W = max(6, min(10, floor(0.2*n)));              % trailing window for local slope

        % --- grow linear region until BOTH residual + slope gates fail ---
        cutoff = n;
        for j = (minPts+W):n
            p_prefix = polyfit(ay(1:j), steer(1:j), 1);
            r        = steer(1:j) - polyval(p_prefix, ay(1:j));

            p_win    = polyfit(ay(j-W+1:j), steer(j-W+1:j), 1);
            slope_dev = abs(p_win(1) - slope_ref);

            if max(abs(r)) > max_error_deg && slope_dev > slope_tol
                cutoff = j - W;            % back off by the window length
                break
            end
        end

        lin_mask = false(size(ay)); lin_mask(1:cutoff) = true;
        p = polyfit(ay(lin_mask), steer(lin_mask), 1);
        Ku(i) = p(1);  % deg/g

        % --- extract car-level metadata safely ---
        hasF = (isstruct(car) && isfield(car,'camber_compliance_f')) || ...
               (isobject(car) && isprop(car,'camber_compliance_f'));
        hasR = (isstruct(car) && isfield(car,'camber_compliance_r')) || ...
               (isobject(car) && isprop(car,'camber_compliance_r'));
        hasM = (isstruct(car) && isfield(car,'M')) || ...
               (isobject(car) && isprop(car,'M'));

        if hasF, Cf(i) = car.camber_compliance_f; end
        if hasR, Cr(i) = car.camber_compliance_r; end
        if hasM, Mass(i) = car.M; end
    end

    % --- filter valid points ---
    valid = ~isnan(Ku) & ~isnan(Cf) & ~isnan(Cr);
    if ~any(valid)
        warning('No valid points to plot.');
        return;
    end
    Cf = Cf(valid);  Cr = Cr(valid);  Ku = Ku(valid);

    % =========================
    %   SURFACE FITS / MODELS
    % =========================

    % 1) Planar LS fit: Ku = b0 + b1*Cf + b2*Cr
    X1 = [ones(numel(Cf),1), Cf(:), Cr(:)];
    b_lin = X1 \ Ku(:);

    % 2) Quadratic LS fit (optional): Ku = c0 + c1*Cf + c2*Cr + c3*Cf^2 + c4*Cr^2 + c5*Cf*Cr
    X2 = [ones(numel(Cf),1), Cf(:), Cr(:), Cf(:).^2, Cr(:).^2, Cf(:).*Cr(:)];
    c_quad = X2 \ Ku(:);

    % 3) Smooth interpolant (no extrapolation outside convex hull)
    F = scatteredInterpolant(Cf(:), Cr(:), Ku(:), 'natural', 'none');

    % --- mesh grid for visualization (expand bounds slightly) ---
    pad = 0.05;
    cf_min = min(Cf); cf_max = max(Cf);
    cr_min = min(Cr); cr_max = max(Cr);
    dcf = max(1e-9, cf_max - cf_min);
    dcr = max(1e-9, cr_max - cr_min);

    cfv = linspace(cf_min - pad*dcf, cf_max + pad*dcf, 60);
    crv = linspace(cr_min - pad*dcr, cr_max + pad*dcr, 60);
    [CF, CR] = meshgrid(cfv, crv);

    % Predictions on grid
    Ku_lin  = b_lin(1) + b_lin(2)*CF + b_lin(3)*CR;
    Ku_quad = c_quad(1) + c_quad(2)*CF + c_quad(3)*CR + c_quad(4)*CF.^2 + c_quad(5)*CR.^2 + c_quad(6).*CF.*CR;
    Ku_nat  = F(CF, CR);   % NaN outside data hull

    % =========================
    %        PLOTTING
    % =========================
    figure('Color','w'); hold on; grid on; box on;

    % Interpolated surface (semi-transparent)
    s1 = surf(CF, CR, Ku_nat, 'EdgeColor','black', 'FaceAlpha',0.55);
   
    % Raw points
    sc = scatter3(Cf, Cr, Ku, 70, Ku, 'filled', 'MarkerEdgeColor',[0 0 0], 'LineWidth',0.5);

    % Labels
    cb = colorbar; ylabel(cb, 'Understeer Gradient (deg/g)');
    xlabel('Front camber compliance (deg/g or deg/kN)');
    ylabel('Rear camber compliance (deg/g or deg/kN)');
    zlabel('Understeer gradient, K_u (deg/g)');
    title('Understeer Gradient Surface vs. Front/Rear Camber Compliance');
    view(135, 25);

    % Legend / annotations (avoid too many handles)
    legend([sc s1], ...
        {'Data (K_u)', 'Interpolant (natural)'}, ...
        'Location','bestoutside');

    % =========================
    %     FIT QUALITY PRINT
    % =========================
    Ku_hat_lin  = X1*b_lin;
    Ku_hat_quad = X2*c_quad;

    SSE_lin  = sum((Ku(:) - Ku_hat_lin(:)).^2);
    SSE_quad = sum((Ku(:) - Ku_hat_quad(:)).^2);
    SST      = sum( (Ku(:) - mean(Ku(:))).^2 );
    R2_lin   = 1 - SSE_lin/SST;
    R2_quad  = 1 - SSE_quad/SST;

    fprintf('\n=== Understeer surface fits ===\n');
    fprintf('Plane fit:   Ku = %.4f + %.4f*Cf + %.4f*Cr   (R^2 = %.4f)\n', b_lin(1), b_lin(2), b_lin(3), R2_lin);
    fprintf('Quad  fit:   Ku = %.4f + %.4f*Cf + %.4f*Cr + %.4f*Cf^2 + %.4f*Cr^2 + %.4f*Cf*Cr   (R^2 = %.4f)\n', ...
        c_quad(1), c_quad(2), c_quad(3), c_quad(4), c_quad(5), c_quad(6), R2_quad);

    % Optional: return models via assignin for quick reuse
    assignin('base','understeer_plane_coeffs', b_lin);
    assignin('base','understeer_quad_coeffs',  c_quad);
    assignin('base','understeer_interpolant',  F);
end
