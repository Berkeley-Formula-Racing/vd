function [] = understeer_plots(carCell)
    % Use first column if passed as an N×M cell matrix
    carCell = carCell(:,1);

    num_cars = numel(carCell);

    % --- prealloc for outputs ---
    vals  = nan(num_cars,1);   % understeer gradient (deg/g)
    Cf    = nan(num_cars,1);   % camber compliance front
    Cr    = nan(num_cars,1);   % camber compliance rear
    Mass  = nan(num_cars,1);   % mass if available
    Names = strings(num_cars,1);

    for i = 1:num_cars
        car   = carCell{i};
        comp  = car.comp;
        ay    = comp.skidpad.ay(:);
        steer = comp.skidpad.steer(:);

        % --- sort ---
        [ay, k]  = sort(ay);
        steer    = steer(k);
        n        = numel(ay);
        if n < 8
            warning('Car %d: not enough points. Skipping.', i);
            continue;
        end

        % --- adaptive thresholds from low-g region ---
        minPts = max(10, ceil(0.25*n));                 % seed length
        p0     = polyfit(ay(1:minPts), steer(1:minPts), 1);
        r0     = steer(1:minPts) - polyval(p0, ay(1:minPts));
        sigma  = 1.4826*mad(r0,1);                       % robust std
        max_error_deg = max(0.15, 3*sigma);              % residual gate

        slope_ref = p0(1);
        slope_tol = max(0.03, 0.20*abs(slope_ref));      % allow >=0.03 deg/g drift
        W = max(6, min(10, floor(0.2*n)));               % trailing window for local slope

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

        % --- metadata for table/labels ---
        hasF = (isstruct(car) && isfield(car,'camber_compliance_f')) || ...
               (isobject(car) && isprop(car,'camber_compliance_f'));
        hasR = (isstruct(car) && isfield(car,'camber_compliance_r')) || ...
               (isobject(car) && isprop(car,'camber_compliance_r'));

        if hasF, Cf(i) = car.camber_compliance_f; end
        if hasR, Cr(i) = car.camber_compliance_r; end

        if (isstruct(car) && isfield(car,'M')) || (isobject(car) && isprop(car,'M'))
            Mass(i) = car.M;
        end

        % Attempt to grab a name/id if present
        if (isstruct(car) && isfield(car,'name')) || (isobject(car) && isprop(car,'name'))
            Names(i) = string(car.name);
        elseif (isstruct(car) && isfield(car,'Name')) || (isobject(car) && isprop(car,'Name'))
            Names(i) = string(car.Name);
        elseif (isstruct(car) && isfield(car,'id')) || (isobject(car) && isprop(car,'id'))
            Names(i) = string(car.id);
        elseif (isstruct(car) && isfield(car,'ID')) || (isobject(car) && isprop(car,'ID'))
            Names(i) = string(car.ID);
        else
            Names(i) = "Car " + i;
        end
    end

    % --- plot summary across cars (colored by Ku) ---
    figure('Color','w'); hold on; grid on; box on;
    x = 1:num_cars;
    s = scatter(x, vals, 40, vals, 'filled'); %#ok<NASGU>
    xlabel('Vehicle index'); ylabel('Understeer gradient (deg/g)');
    title('Understeer Gradient (linear region)');
    xlim([0, num_cars+1]);
    cb = colorbar; cb.Label.String = 'Ku (deg/g)';

    pfit = polyfit(x(isfinite(vals)), vals(isfinite(vals)), 1);
    disp('Slope of Ku vs. index:'), disp(pfit(1));

    % --- light label overlay (optional) ---
    for i = 1:num_cars
        if ~isfinite(vals(i)), continue; end
        lbl = sprintf('Cf=%.3f, Cr=%.3f, m=%g', Cf(i), Cr(i), Mass(i));
        text(x(i)+0.05, vals(i), lbl, 'FontSize',8, ...
             'HorizontalAlignment','left','VerticalAlignment','bottom');
    end

    % --- build & save the table ---
    T = table((1:num_cars)', Names, vals, Cf, Cr, Mass, ...
        'VariableNames', {'Index','Name','Ku_deg_per_g','Cf_deg_per_g','Cr_deg_per_g','Mass'});

    % expose to base workspace and save files
    assignin('base','understeer_table', T);
    try
        writetable(T, 'understeer_gradients.csv');
        save('understeer_gradients.mat', 'T');
    catch ME
        warning('Could not save table to disk: %s', ME.message);
    end
end
