function [model, out] = fit_points_generic(carCell, predictorFields, opts)
%FIT_POINTS_GENERIC  OLS of Points vs arbitrary predictors.
%
% Usage:
%   [model, out] = fit_points_generic(carCell, {'M','camber_compliance_f','camber_compliance_r'});
%   [model, out] = fit_points_generic(carCell, {'M','susp.roll_grad','aero.CL'}, ...
%                    struct('useInteractions',true,'useQuadratic',false,'standardize',true));
%
% Inputs:
%   carCell         : {N×1} cell array of car structs/objects
%   predictorFields : 1×P cellstr of field paths to use as predictors.
%                     Supports nested paths using dots, e.g. 'susp.roll_grad'
%   opts (optional) : struct with fields
%       .responseField   : default 'comp.points.total'
%       .useInteractions : default false
%       .useQuadratic    : default false
%       .standardize     : default false
%       .verbose         : default true
%
% Outputs:
%   model: struct with fields
%       .coef, .se, .t, .p, .R2, .adjR2, .RMSE, .dof, .term_names
%       .mu, .sg, .spec, .predictFcn
%   out: struct with fields
%       .data      : table of predictors + response + Name
%       .keep_mask : logical keep rows used
%       .diag      : table of Observed/Predicted/Residual by row
%
% Notes:
%   - Only numeric predictors are supported (categorical encoding not included).
%   - Missing/non-finite rows are dropped consistently across all variables.

if nargin < 3, opts = struct; end
opts = set_default(opts, 'responseField',   'comp.points.total');
opts = set_default(opts, 'useInteractions', false);
opts = set_default(opts, 'useQuadratic',    false);
opts = set_default(opts, 'standardize',     false);
opts = set_default(opts, 'verbose',         true);

carCell = carCell(:,1);
N = numel(carCell);
P = numel(predictorFields);

% -----------------------------
% Extract predictors & response
% -----------------------------
Xraw = nan(N, P);
Name = strings(N,1);
for i = 1:N
    c = carCell{i};
    for j = 1:P
        Xraw(i,j) = getpath(c, predictorFields{j});
    end
    Name(i) = "Car " + i;
end
y = nan(N,1);
for i = 1:N
    y(i) = getpath(carCell{i}, opts.responseField);
end

% Assemble table (for user visibility)
T = array2table(Xraw, 'VariableNames', sanitize_names(predictorFields));
T.Points = y;
T.Name = Name;

% Valid rows = all finite predictors & response
keep = all(isfinite(Xraw),2) & isfinite(y);
if nnz(keep) < max(P+1, 3)
    error('Not enough valid rows after filtering. Have %d, need at least %d.', nnz(keep), max(P+1,3));
end

% -----------------------------
% Standardize if requested
% -----------------------------
mu = mean(Xraw(keep,:),1);
sg = std(Xraw(keep,:),0,1);  sg(sg==0) = 1;
if opts.standardize
    B = (Xraw - mu)./sg;
    base_names = strcat('z_', sanitize_names(predictorFields));
else
    B = Xraw;
    base_names = sanitize_names(predictorFields);
end

% -----------------------------
% Build design matrix
% -----------------------------
[X, term_names] = build_X(B(keep,:), opts.useInteractions, opts.useQuadratic, base_names);

% -----------------------------
% OLS fit
% -----------------------------
yfit = y(keep);
n = size(X,1); k = size(X,2);
coef = X \ yfit;
res  = yfit - X*coef;
dof  = n - k;
s2   = (res' * res) / max(dof,1);
COVb = s2 * inv(X' * X);
se   = sqrt(diag(COVb));
tval = coef ./ se;
pval = 2 * (1 - tcdf(abs(tval), max(dof,1)));
SStot = sum((yfit - mean(yfit)).^2);
SSres = sum(res.^2);
R2    = 1 - SSres/SStot;
adjR2 = 1 - (1-R2)*(n-1)/max(dof,1);
RMSE  = sqrt(s2);

% -----------------------------
% Predictor wrapper(s)
% -----------------------------
spec = struct('useInteractions',opts.useInteractions, ...
              'useQuadratic',opts.useQuadratic, ...
              'standardize',opts.standardize, ...
              'mu',mu,'sg',sg,'predictorFields',{predictorFields}, ...
              'base_names',{base_names}, ...
              'term_names',{term_names});

predictFcn = @(varargin) predict_any(varargin{:}, coef, spec);

% Diagnostics table (on kept rows)
yhat_keep = X*coef;
diagTbl = table(T.Name(keep), yfit, yhat_keep, yfit - yhat_keep, ...
    'VariableNames', {'Name','Observed','Predicted','Residual'});

% Pack outputs
model = struct('coef',coef,'se',se,'t',tval,'p',pval,'R2',R2,'adjR2',adjR2, ...
               'RMSE',RMSE,'dof',dof,'term_names',{term_names}, ...
               'mu',mu,'sg',sg,'spec',spec,'predictFcn',predictFcn);

out = struct('data',T,'keep_mask',keep,'diag',diagTbl);

% Print
if opts.verbose
    fprintf('\n--- Generic Points model ---\n');
    fprintf('Predictors: %s\n', strjoin(predictorFields, ', '));
    fprintf('Response  : %s\n', opts.responseField);
    fprintf('R^2 = %.4f (adj %.4f), RMSE = %.4f, dof = %d\n', R2, adjR2, RMSE, dof);
    for i = 1:numel(coef)
        fprintf('  %-12s %+ .6g   (se=%.3g, t=%.2f, p=%.3g)\n', term_names{i}, coef(i), se(i), tval(i), pval(i));
    end
end

end

% =======================
% Helpers
% =======================
function yhat = predict_any(arg, coef, spec)
% Accepts:
%   - numeric matrix (N×P) in the order of spec.predictorFields
%   - struct scalar/array with those fields
%   - table with those variable names (or columns matching sanitize_names)
%
% Returns:
%   yhat (N×1)
    if istable(arg)
        Xraw = table2array(arg(:, sanitize_names(spec.predictorFields)));
    elseif isnumeric(arg)
        Xraw = arg;
    elseif isstruct(arg)
        Xraw = extract_from_struct_array(arg, spec.predictorFields);
    else
        error('Unsupported input to predictFcn. Use numeric matrix, table, or struct array.');
    end
    if size(Xraw,2) ~= numel(spec.predictorFields)
        error('Predict input has %d columns, expected %d.', size(Xraw,2), numel(spec.predictorFields));
    end
    % standardize if needed
    if spec.standardize
        B = (Xraw - spec.mu)./spec.sg;
    else
        B = Xraw;
    end
    % build design and evaluate
    [X, ~] = build_X(B, spec.useInteractions, spec.useQuadratic, spec.base_names);
    yhat = X * coef;
end

function [X, names] = build_X(B, useInter, useQuad, base_names)
    n = size(B,1); p = size(B,2);
    X = [ones(n,1), B];
    names = [{'1'}, base_names];
    % interactions
    if useInter
        for i = 1:p
            for j = (i+1):p
                X = [X, B(:,i).*B(:,j)]; %#ok<AGROW>
                names{end+1} = sprintf('%s:%s', base_names{i}, base_names{j}); %#ok<AGROW>
            end
        end
    end
    % quadratics
    if useQuad
        for i = 1:p
            X = [X, B(:,i).^2]; %#ok<AGROW>
            names{end+1} = sprintf('%s^2', base_names{i}); %#ok<AGROW>
        end
    end
end

function a = extract_from_struct_array(S, fields)
    n = numel(S); p = numel(fields);
    a = nan(n, p);
    for i = 1:n
        for j = 1:p
            a(i,j) = getpath(S(i), fields{j});
        end
    end
end

function v = getpath(s, path)
% safely get nested field via dot path; returns NaN if missing/non-numeric
    try
        parts = split(string(path), '.');
        v = s;
        for k = 1:numel(parts)
            pk = strtrim(parts{k});
            if isstruct(v) && isfield(v, pk)
                v = v.(pk);
            elseif isobject(v) && isprop(v, pk)
                v = v.(pk);
            else
                v = NaN; return;
            end
        end
        if ~isscalar(v) || ~isnumeric(v) || ~isfinite(v)
            if isnumeric(v) && isscalar(v)
                v = double(v);
            else
                v = NaN;
            end
        end
    catch
        v = NaN;
    end
end

function tf = hasfield(s, f)
    tf = (isstruct(s) && isfield(s, f)) || (isobject(s) && isprop(s, f));
end

function names = sanitize_names(fields)
    names = regexprep(fields, '[^A-Za-z0-9_]', '_');
end

function s = set_default(s, f, v)
    if ~isfield(s,f) || isempty(s.(f)), s.(f) = v; end
end
