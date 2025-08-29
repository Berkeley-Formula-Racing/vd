function [] = understeer_plots(carCell)
    % Use first column if passed as an N×M cell matrix
    carCell = carCell(:,1);

    num_cars = numel(carCell);

    vals = nan(1, num_cars);  % understeer gradient (deg/g)

    for i = 1:num_cars
        comp  = carCell{i}.comp;
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
        vals(i) = p(1);                     % Ku (deg/g)
    end

    % --- plot summary across cars ---
    figure; hold on; grid on; box on;
    x = 1:num_cars;
    scatter(x, vals, 40, 'filled');
    xlabel('Vehicle index'); ylabel('Understeer gradient (deg/g)');
    title('Understeer Gradient (linear region)');
    xlim([0, num_cars+1]);

    p = polyfit(x, vals, 1);
    disp(p(1));

    % --- labels: camber_compliance_f / camber_compliance_r ---
    for i = 1:num_cars
        if isnan(vals(i)), continue; end
        car = carCell{i};
        % handle struct vs object safely
        hasF = (isstruct(car) && isfield(car,'camber_compliance_f')) || ...
               (isobject(car) && isprop(car,'camber_compliance_f'));
        hasR = (isstruct(car) && isfield(car,'camber_compliance_r')) || ...
               (isobject(car) && isprop(car,'camber_compliance_r'));
        cf = NaN; cr = NaN;
        if hasF, cf = car.camber_compliance_f; end
        if hasR, cr = car.camber_compliance_r; end

        lbl = sprintf('Cf=%.3f, Cr=%.3f', cf, cr);
        text(x(i)+0.05, vals(i), lbl, 'FontSize',8, ...
             'HorizontalAlignment','left','VerticalAlignment','bottom');
    end
    disp(vals);
end
