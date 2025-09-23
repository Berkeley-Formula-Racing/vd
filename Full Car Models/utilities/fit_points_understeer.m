function [models, tables] = fit_points_understeer(carCell, opts)
%FIT_POINTS_UNDERSTEER  Fit Points vs. Mass, ccF, ccR.
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
%       - comp.points.total (scalar)
%
%   opts (optional): struct with fields
%       .useInteractions : include pairwise interaction terms (default=true)
%       .useQuadratic    : include squared terms (default=false)
%       .standardize     : z-score predictors before fitting (default=false)
%
% Outputs:
%   models: struct with field:
%       .Points: .coef, .se, .t, .p, .R2, .adjR2, .RMSE, .dof, .predictFcn, .spec
%   tables: struct with:
%       .data : table of predictors+response per car
%       .diag_Points : table with fitted vs residuals
%
% Notes:
%   - If opts.standardize is true, coefficients are on z-scales; the returned
%     predictFcn handles raw (unscaled) inputs transparently.
%   - You can toggle interactions/quadratic terms without changing call sites.

if nargin < 2 || isempty(opts), opts = struct; end
opts = set_default(opts, 'useInteractions', true);
opts = set_default(opts, 'useQuadratic',    false);
opts = set_default(opts, 'standardize',     false);

% Ensure column cell
carCell = carCell(:,1);
N = numel(carCell);

% -----------------------------
% Collect predictors & response
% -----------------------------
M   = nan(N,1);
ccF = nan(N,1);
ccR = nan(N,1);
Pts = nan(N,1);
Name = strings(N,1);

for i = 1:N
    car = carCell{i};

    M(i)  = car.M;
    ccF(i) = car.camber_compliance_f;
    ccR(i) = car.camber_compliance_r;
    Pts(i) = car.comp.points.total;

    Name(i) = "Car " + i;
end

T = table(Name, M, ccF, ccR, Pts);

% Valid rows
keep_pts = isfinite(T.M) & isfinite(T.ccF) & isfinite(T.ccR) & isfinite(T.Pts);
if nnz(keep_pts) < 3
    error('Not enough valid rows for Points model (need >= 3).');
end

% -----------------------------
% Design spec and matrices
% -----------------------------
[spec, X_P, y_P, mu, sg] = build_design(T(keep_pts,:), 'Pts', opts);

% -----------------------------
% Fit OLS for Points
% -----------------------------
models = struct; tables = struct; tables.data = T;

models.Points = do_ols_fit(X_P, y_P, spec, mu, sg, opts);
tables.diag_Points = make_diag_table(T(keep_pts,:), y_P, X_P*models.Points.coef);

% Pretty print
fprintf('\n=== Model spec ===\n'); disp(spec);
print_model('Points', models.Points, spec);

% Quick visualization at median mass slice
try
    do_quick_plots(T, models, spec);
catch ME
    warning('Plotting failed: %s', ME.message);
end

end

%=====================================================================
% Build design matrix X given options.
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
        spec = spec_in; % reuse existing spec
        [~, X] = realize_spec(base_cols, spec);
    else
        spec = struct('useInteractions',opts.useInteractions, ...
                      'useQuadratic',opts.useQuadratic, ...
                      'standardize',opts.standardize, ...
                      'mu',mu,'sg',sg,'pred_names',{pred_names}, ...
                      'term_names',{names});
        [spec, X] = realize_spec(base_cols, spec);
    end
end

% Canonicalize term order and provide a function for future predictions
function [spec, X] = realize_spec(B, spec)
    % B: base predictor matrix after (optional) standardization
    % Map base predictors in fixed order: x1->M, x2->ccF, x3->ccR
    names = {'1','x1','x2','x3'}; 

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
% Quick slice over ccF/ccR at median Mass for visualization
%=====================================================================
function do_quick_plots(T, models, spec)
    medM = median(T.M, 'omitnan');
    cfv = linspace(min(T.ccF), max(T.ccF), 41);
    crv = linspace(min(T.ccR), max(T.ccR), 41);
    [CF, CR] = meshgrid(cfv, crv);

    % Points surface
    mdl = models.Points;
    figure('Color','w'); hold on; grid on; box on;
    Z = mdl.predictFcn(repmat(medM, numel(CF), 1), CF(:), CR(:));
    Z = reshape(Z, size(CF));
    surf(CF, CR, Z, 'EdgeAlpha', 0.15, 'FaceAlpha', 0.85);

    scatter3(T.ccF, T.ccR, T.Pts, 28, T.Pts, ...
        'filled', 'MarkerEdgeColor',[0 0 0.35]);  % data points

    xlabel('ccF (deg/g)'); ylabel('ccR (deg/g)'); zlabel('Points');
    title(sprintf('Points vs ccF, ccR @ M=%.1f (spec: int=%d, quad=%d, std=%d)', ...
        medM, spec.useInteractions, spec.useQuadratic, spec.standardize));
    view(40,28);
    colorbar;
end
