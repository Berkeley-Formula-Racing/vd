function fit = lc0NDFitLongitudinal(donor,opts)
%LC0NDFITLONGITUDINAL Build a load-normalized donor Fx reference curve.
%   The curve uses only near-zero-slip-angle samples at an explicit
%   reference load, pressure, and camber.  It is a donor characterization,
%   not a claim that the donor force level is the target LC0 force level.

requiredData = {'slipRatio','slipAngle_deg','Fx_N','Fz_N', ...
    'pressure_psi','camber_deg'};
requiredOpts = {'load_N','loadHalfWidth_N','pressure_psi', ...
    'pressureHalfWidth_psi','camber_deg','camberHalfWidth_deg', ...
    'maxAbsSlipAngle_deg','kappaGrid','kappaHalfWidth','minSamplesPerBin'};
for i = 1:numel(requiredData)
    if ~isfield(donor,requiredData{i})
        error('lc0NDFitLongitudinal:missingField', ...
            'DONOR must contain %s.',requiredData{i});
    end
end
for i = 1:numel(requiredOpts)
    if ~isfield(opts,requiredOpts{i})
        error('lc0NDFitLongitudinal:missingOption', ...
            'OPTS must contain %s.',requiredOpts{i});
    end
end

kappa = donor.slipRatio(:);
alpha = donor.slipAngle_deg(:);
fx = donor.Fx_N(:);
fz = donor.Fz_N(:);
pressure = donor.pressure_psi(:);
camber = donor.camber_deg(:);
n = numel(kappa);
if any([numel(alpha),numel(fx),numel(fz),numel(pressure),numel(camber)] ~= n)
    error('lc0NDFitLongitudinal:sizeMismatch', ...
        'All DONOR channels must have equal length.');
end
if any(fz <= 0)
    error('lc0NDFitLongitudinal:nonpositiveLoad', ...
        'DONOR Fz_N must be positive.');
end

use = isfinite(kappa) & isfinite(alpha) & isfinite(fx) & isfinite(fz) & ...
    isfinite(pressure) & isfinite(camber) & ...
    abs(alpha) <= opts.maxAbsSlipAngle_deg & ...
    abs(fz-opts.load_N) <= opts.loadHalfWidth_N & ...
    abs(pressure-opts.pressure_psi) <= opts.pressureHalfWidth_psi & ...
    abs(camber-opts.camber_deg) <= opts.camberHalfWidth_deg;
muX = fx./fz;
grid = opts.kappaGrid(:);
mu = nan(numel(grid),1);
sigma = nan(numel(grid),1);
count = zeros(numel(grid),1);
for i = 1:numel(grid)
    inBin = use & abs(kappa-grid(i)) <= opts.kappaHalfWidth;
    count(i) = sum(inBin);
    if count(i) > 0
        mu(i) = mean(muX(inBin));
        sigma(i) = std(muX(inBin));
    end
end
qualified = count >= opts.minSamplesPerBin;
curve = table(grid,mu,sigma,count,qualified,'VariableNames', ...
    {'slip_ratio','mu_x','mu_x_std','n_samples','is_qualified'});
fit = struct('options',opts,'n_selected',sum(use),'curve',curve, ...
    'peak_drive_mu',peak(mu,qualified,1),'peak_brake_mu',peak(mu,qualified,-1));
end

function value = peak(mu,qualified,sign)
values = mu(qualified);
if isempty(values)
    value = NaN;
elseif sign > 0
    value = max(values);
else
    value = -min(values);
end
end
