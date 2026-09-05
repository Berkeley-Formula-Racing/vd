function fit = lc0NDFitCouplingExponent(combined,longFit,latFit,opts)
%LC0NDFITCOUPLINGEXPONENT Fit a bounded friction-ellipse exponent to donor data.
%   Pure normalized donor curves define the requested Fx/Fy demands.  Each
%   exponent is scored against measured combined force, with no force-level
%   transfer to the target LC0 tyre implied by this fit.

requiredCombined = {'slipRatio','slipAngle_deg','Fx_N','Fy_N','Fz_N', ...
    'pressure_psi','camber_deg'};
requiredOpts = {'load_N','loadHalfWidth_N','pressure_psi', ...
    'pressureHalfWidth_psi','camber_deg','camberHalfWidth_deg', ...
    'minAbsSlipRatio','minAbsSlipAngle_deg','pGrid','minPoints'};
for i = 1:numel(requiredCombined)
    if ~isfield(combined,requiredCombined{i})
        error('lc0NDFitCouplingExponent:missingField', ...
            'COMBINED must contain %s.',requiredCombined{i});
    end
end
for i = 1:numel(requiredOpts)
    if ~isfield(opts,requiredOpts{i})
        error('lc0NDFitCouplingExponent:missingOption', ...
            'OPTS must contain %s.',requiredOpts{i});
    end
end

kappa = combined.slipRatio(:); alpha = combined.slipAngle_deg(:);
fx = combined.Fx_N(:); fy = combined.Fy_N(:); fz = combined.Fz_N(:);
pressure = combined.pressure_psi(:); camber = combined.camber_deg(:);
n = numel(kappa);
if any([numel(alpha),numel(fx),numel(fy),numel(fz),numel(pressure),numel(camber)] ~= n)
    error('lc0NDFitCouplingExponent:sizeMismatch', ...
        'All COMBINED channels must have equal length.');
end
if any(fz <= 0)
    error('lc0NDFitCouplingExponent:nonpositiveLoad','COMBINED Fz_N must be positive.');
end

xCurve = longFit.curve(longFit.curve.is_qualified,:);
yCurve = latFit.curve(latFit.curve.is_qualified,:);
if height(xCurve) < 2 || height(yCurve) < 2
    error('lc0NDFitCouplingExponent:insufficientPureCurve', ...
        'Both qualified pure reference curves require at least two points.');
end
muX0 = interp1(xCurve.slip_ratio,xCurve.mu_x,kappa,'linear',NaN);
muY0 = interp1(yCurve.slip_angle_deg,yCurve.mu_y,alpha,'linear',NaN);
use = isfinite(kappa) & isfinite(alpha) & isfinite(fx) & isfinite(fy) & ...
    isfinite(fz) & isfinite(pressure) & isfinite(camber) & ...
    abs(kappa) >= opts.minAbsSlipRatio & ...
    abs(alpha) >= opts.minAbsSlipAngle_deg & ...
    abs(fz-opts.load_N) <= opts.loadHalfWidth_N & ...
    abs(pressure-opts.pressure_psi) <= opts.pressureHalfWidth_psi & ...
    abs(camber-opts.camber_deg) <= opts.camberHalfWidth_deg & ...
    isfinite(muX0) & isfinite(muY0);
muXMeasured = fx(use)./fz(use);
muYMeasured = fy(use)./fz(use);
muX0 = muX0(use); muY0 = muY0(use);
capX = max(abs(xCurve.mu_x)); capY = max(abs(yCurve.mu_y));
pGrid = opts.pGrid(:);
rmse = nan(numel(pGrid),1);
for i = 1:numel(pGrid)
    p = pGrid(i);
    if ~isfinite(p) || p < 1, continue, end
    utilization = (abs(muX0/capX).^p + abs(muY0/capY).^p).^(1/p);
    scale = min(1,1./max(utilization,eps));
    errorVector = [scale.*muX0-muXMeasured; scale.*muY0-muYMeasured];
    rmse(i) = sqrt(mean(errorVector.^2));
end
qualified = sum(use) >= opts.minPoints;
if qualified && any(isfinite(rmse))
    [bestRmse,index] = min(rmse);
    bestExponent = pGrid(index);
else
    bestRmse = NaN;
    bestExponent = NaN;
end
fit = struct('options',opts,'n_points',sum(use),'is_qualified',qualified, ...
    'mu_x_capacity',capX,'mu_y_capacity',capY, ...
    'score',table(pGrid,rmse,'VariableNames',{'exponent','rmse_mu'}), ...
    'best_exponent',bestExponent,'rmse_mu_at_best',bestRmse);
end
