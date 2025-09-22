function [models, tables] = fit_points_understeer(carCell, opts)
%FIT_POINTS_UNDERSTEER  Fit Points and Understeer Gradient vs. Mass, ccF, ccR.
%
% Usage:
%   [models, tables] = fit_points_understeer(carCell)
%   [models, tables] = fit_points_understeer(carCell, struct('useInteractions',true,'useQuadratic',false,'standardize',false))
%
% Inputs:
%   carCell : {N×1} cell array of car structs/objects with fields:
%       - camber_compliance_f (deg/g)
%       - camber_compliance_r (deg/g)
%       - M (kg) or Mass (kg)   [optional but recommended]
%       - comp.skidpad.ay  (vector, g)
%       - comp.skidpad.steer (vector, deg)
%       - comp.points.total (scalar)
%
%   opts (optional): struct with fields
%       .useInteractions : include pairwise interaction terms (default=true)
%       .useQuadratic    : include squared terms (default=false)
%       .standardize     : z-score predictors before fitting (default=false)
%
% Outputs:
%   models: struct with fields for each response (Points, Ku), each containing:
%       .coef, .se, .t, .p, .R2, .adjR2, .RMSE, .dof, .predictFcn, .spec
%   tables: struct with underlying data table and per-model diagnostics tables
%       .data : table of predictors+responses per car
%       .diag_Points, .diag_Ku : tables with fitted vs residuals
%
% Notes:
%   - Understeer gradient Ku is estimated using an adaptive linear-region
%     detector identical to your understeer_surface() routine.
%   - If opts.standardize is true, coefficients are on z-scales; the returned
%     predictFcn handles raw (unscaled) inputs transparently.
%   - You can toggle interactions/quadratic terms without changing call sites.


if nargin < 2 || isempty(opts)
    opts = struct; 
end
opts = set_default(opts, 'useInteractions', true);
opts = set_default(opts, 'useQuadratic',    false);
opts = set_default(opts, 'standardize',     false);

% Ensure column cell
carCell = carCell(:,1);
N = numel(carCell);

% -----------------------------
% Collect predictors & responses
% -----------------------------
M   = nan(N,1);
ccF = nan(N,1);
ccR = nan(N,1);
Pts = nan(N,1);
Ku  = nan(N,1);
Name = strings(N,1);

for i = 1:N
    car = carCell{i};

    % mass (try multiple common names)
    if hasFieldOrProp(car,'M'),        M(i)   = getFieldOrProp(car,'M');        end
    if isnan(M(i)) && hasFieldOrProp(car,'Mass'), M(i) = getFieldOrProp(car,'Mass'); end

    if hasFieldOrProp(car,'camber_compliance_f'), ccF(i) = getFieldOrProp(car,'camber_compliance_f'); end
    if hasFieldOrProp(car,'camber_compliance_r'), ccR(i) = getFieldOrProp(car,'camber_compliance_r'); end

    % points
    try
        Pts(i) = car.comp.points.total;
    catch
        % leave NaN
    end

    % name
    if     hasFieldOrProp(car,'name'), Name(i) = string(getFieldOrProp(car,'name'));
    elseif hasFieldOrProp(car,'Name'), Name(i) = string(getFieldOrProp(car,'Name'));
    else,                              Name(i) = "Car " + i; end

    % understeer gradient from skidpad (deg/g)
    try
        ay    = car.comp.skidpad.ay(:);
        steer = car.comp.skidpad.steer(:);
        [Ku(i),~] = estimate_understeer(ay, steer);
    catch
        % leave NaN
    end
end

T = table(Name, M, ccF, ccR, Pts, Ku);

% Valid rows per response
keep_pts = isfinite(T.M) & isfinite(T.ccF) & isfinite(T.ccR) & isfinite(T.Pts);
keep_ku  = isfinite(T.M) & isfinite(T.ccF) & isfinite(T.ccR) & isfinite(T.Ku);

if nnz(keep_pts) < 3 && nnz(keep_ku) < 3
    error('Not enough valid rows for either response.');
end

% -----------------------------
% Design spec and matrices
% -----------------------------
[spec, X_P, y_P, mu, sg] = build_design(T(keep_pts,:), 'Pts', opts);
[~,    X_K, y_K, ~,  ~ ] = build_design(T(keep_ku ,:), 'Ku',  opts, spec, mu, sg);

% -----------------------------
% Fit OLS for both responses
% -----------------------------
models = struct; tables = struct; tables.data = T;
if ~isempty(X_P)
    models.Points = do_ols_fit(X_P, y_P, spec, mu, sg, opts);
    tables.diag_Points = make_diag_table(T(keep_pts,:), y_P, X_P*models.Points.coef);
end
if ~isempty(X_K)
    models.Ku = do_ols_fit(X_K, y_K, spec, mu, sg, opts);
    tables.diag_Ku = make_diag_table(T(keep_ku,:), y_K, X_K*models.Ku.coef);
end

% Pretty print
fprintf('\n=== Model spec ===\n'); disp(spec);
if isfield(models,'Points');  print_model('Points', models.Points, spec);  end
if isfield(models,'Ku');      print_model('Ku',     models.Ku,     spec);  end

% Quick visualization at median mass slice
try
    do_quick_plots(T, models, spec);
catch ME
    warning('Plotting failed: %s', ME.message);
end

end

%=====================================================================
% Estimator for understeer gradient using adaptive linear region growth
%=====================================================================
function [Ku, meta] = estimate_understeer(ay, steer)
    Ku = NaN; meta = struct('cutoff',NaN,'n',numel(ay));
    if numel(ay) < 8 || numel(steer) < 8
        return
    end

    [ay, k] = sort(ay(:)); steer = steer(:); steer = steer(k);
    n = numel(ay);

    minPts = max(10, ceil(0.25*n));
    p0   = polyfit(ay(1:minPts), steer(1:minPts), 1);
    r0   = steer(1:minPts) - polyval(p0, ay(1:minPts));
    sigma = 1.4826*mad(r0,1);
    max_error_deg = max(0.15, 3*sigma);

    slope_ref = p0(1);
    slope_tol = max(0.03, 0.20*abs(slope_ref));
    W = max(6, min(10, floor(0.2*n)));

    cutoff = n;
    for j = (minPts+W):n
        p_prefix = polyfit(ay(1:j), steer(1:j), 1);
        r        = steer(1:j) - polyval(p_prefix, ay(1:j));

        p_win    = polyfit(ay(j-W+1:j), steer(j-W+1:j), 1);
        slope_dev = abs(p_win(1) - slope_ref);

        if max(abs(r)) > max_error_deg && slope_dev > slope_tol
            cutoff = j - W; break
        end
    end

    idx = false(n,1); idx(1:cutoff) = true;
    p = polyfit(ay(idx), steer(idx), 1);
    Ku = p(1);
    meta.cutoff = cutoff;
end

%=====================================================================
% Build design matrix X given options. Supports reuse of spec/mu/sg.
%=====================================================================
function [spec, X, y, mu, sg] = build_design(T, yname, opts, spec_in, mu_in, sg_in)
    % Predictors (raw)
    P = [T.M, T.ccF, T.ccR];
    pred_names = {'M','ccF','ccR'};

    % Standardize if requested
    if nargin >= 5 && ~isempty(mu_in)
        mu = mu_in; sg = sg_in;
    else
        mu = mean(P,1,'omitnan');
        sg = std(P,0,1,'omitnan'); sg(sg==0) = 1;
    end
    if opts.standardize
        Z = (P - mu)./sg;
        base_cols = Z;
        base_names = strcat('z_', pred_names);
    else
        base_cols = P;
        base_names = pred_names;
    end

    % Start with intercept + main effects
    X = [ones(size(T,1),1), base_cols];
    names = [{'1'}, base_names];

    if opts.useInteractions
        % pairwise interactions
        k = size(base_cols,2);
        for a = 1:k
            for b = (a+1):k
                X = [X, base_cols(:,a).*base_cols(:,b)]; %#ok<AGROW>
                names{end+1} = sprintf('%s:%s', base_names{a}, base_names{b}); %#ok<AGROW>
            end
        end
    end

    if opts.useQuadratic
        for a = 1:size(base_cols,2)
            X = [X, base_cols(:,a).^2]; %#ok<AGROW>
            names{end+1} = sprintf('%s^2', base_names{a}); %#ok<AGROW>
        end
    end

    y = T.(yname);

    if nargin >= 4 && ~isempty(spec_in)
        spec = spec_in; % reuse existing spec so columns align across responses
        % If spec_in provided, rebuild X to match it exactly
        [~, X] = realize_spec(base_cols, spec);
    else
        spec = struct('useInteractions',opts.useInteractions, ...
                      'useQuadratic',opts.useQuadratic, ...
                      'standardize',opts.standardize, ...
                      'mu',mu,'sg',sg,'pred_names',{pred_names}, ...
                      'term_names',{names});
        % Freeze a canonical realization function
        [spec, X] = realize_spec(base_cols, spec);
    end
end

% Canonicalize term order and provide a function for future predictions
function [spec, X] = realize_spec(B, spec)
    % B: base predictor matrix after (optional) standardization
    pnames = spec.pred_names; %#ok<NASGU>
    names = {'1','x1','x2','x3'}; % map: 1 -> intercept; x1->M, x2->ccF, x3->ccR

    % base
    X = [ones(size(B,1),1), B(:,1), B(:,2), B(:,3)];

    % interactions (pairwise)
    if spec.useInteractions
        X = [X, B(:,1).*B(:,2), B(:,1).*B(:,3), B(:,2).*B(:,3)];
        names = [names, {'x1:x2','x1:x3','x2:x3'}]; 
    end

    % quadratic
    if spec.useQuadratic
        X = [X, B(:,1).^2, B(:,2).^2, B(:,3).^2];
        names = [names, {'x1^2','x2^2','x3^2'}];
    end

    spec.term_names = names;
end

%=====================================================================
% OLS core + wrapped predictor for raw inputs
%=====================================================================
function model = do_ols_fit(X, y, spec, mu, sg, opts)
    n = size(X,1); k = size(X,2);
    coef = X \ y;
    res  = y - X*coef;
    dof  = n - k;
    s2   = (res' * res) / max(dof,1);
    COVb = s2 * inv(X' * X);
    se   = sqrt(diag(COVb));
    tval = coef ./ se;
    pval = 2 * (1 - tcdf(abs(tval), max(dof,1)));
    SS_tot = sum((y - mean(y)).^2);
    SS_res = sum(res.^2);
    R2     = 1 - SS_res/SS_tot;
    adjR2  = 1 - (1-R2)*(n-1)/max(dof,1);
    RMSE   = sqrt(s2);

    % prediction wrapper expects RAW M, ccF, ccR
    function yhat = predict_raw(Mr, ccFr, ccRr)
        B = [Mr, ccFr, ccRr];
        if opts.standardize
            B = (B - mu)./sg;
        end
        % realize spec for a single row (or vectorized)
        [~, Xr] = realize_spec(B, spec);
        yhat = Xr * coef;
    end

    model = struct('coef',coef,'se',se,'t',tval,'p',pval,'R2',R2, ...
                   'adjR2',adjR2,'RMSE',RMSE,'dof',dof, ...
                   'predictFcn',@predict_raw,'spec',spec);
end

function D = make_diag_table(T, y, yhat)
    D = table(T.Name, T.M, T.ccF, T.ccR, y, yhat, y - yhat, ...
        'VariableNames', {'Name','M','ccF','ccR','Observed','Predicted','Residual'});
end

function print_model(label, model, spec)
    fprintf('\n--- %s model ---\n', label);
    fprintf('R^2 = %.4f (adj %.4f), RMSE = %.4f, dof = %d\n', model.R2, model.adjR2, model.RMSE, model.dof);
    for i = 1:numel(model.coef)
        fprintf('  %-8s  %+ .6g   (se=%.3g, t=%.2f, p=%.3g)\n', spec.term_names{i}, model.coef(i), model.se(i), model.t(i), model.p(i));
    end
end

function s = set_default(s, f, v)
    if ~isfield(s,f) || isempty(s.(f)), s.(f) = v; end
end

%=====================================================================
% Quick slices over ccF/ccR at median Mass for visualization
%=====================================================================
function do_quick_plots(T, models, spec)
    medM = median(T.M, 'omitnan');
    cfv = linspace(min(T.ccF), max(T.ccF), 41);
    crv = linspace(min(T.ccR), max(T.ccR), 41);
    [CF, CR] = meshgrid(cfv, crv);

    if isfield(models,'Points'); do_plot(models.Points, 'Pts'); end
    if isfield(models,'Ku');     do_plot(models.Ku,     'Ku');  end

    % ---- nested helper ----
    function do_plot(mdl, ttl)
        figure('Color','w'); hold on; grid on; box on;
        Z = mdl.predictFcn(repmat(medM, numel(CF), 1), CF(:), CR(:));
        Z = reshape(Z, size(CF));
        surf(CF, CR, Z, 'EdgeAlpha', 0.15, 'FaceAlpha', 0.85);

        scatter3(T.ccF, T.ccR, T.(ttl), 28, T.(ttl), ...
            'filled', 'MarkerEdgeColor',[0 0 0.35]);  % data points

        xlabel('ccF (deg/g)'); ylabel('ccR (deg/g)'); zlabel(ttl);
        title(sprintf('%s vs ccF, ccR @ M=%.1f (spec: int=%d, quad=%d, std=%d)', ...
            ttl, medM, spec.useInteractions, spec.useQuadratic, spec.standardize));
        view(40,28);
        colorbar;
    end
end

% ----------------- helpers used above -----------------
function tf = hasFieldOrProp(S, name)
    tf = (isstruct(S) && isfield(S,name)) || (isobject(S) && isprop(S,name));
end

function v = getFieldOrProp(S, name)
    if isstruct(S), v = S.(name); else, v = S.(name); end
end
