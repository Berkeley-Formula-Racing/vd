function fit = lc0NDFitLateral(donor,opts)
%LC0NDFITLATERAL Build a load-normalized pure lateral donor reference.

requiredData = {'slipRatio','slipAngle_deg','Fy_N','Fz_N', ...
    'pressure_psi','camber_deg'};
requiredOpts = {'load_N','loadHalfWidth_N','pressure_psi', ...
    'pressureHalfWidth_psi','camber_deg','camberHalfWidth_deg', ...
    'maxAbsSlipRatio','alphaGrid','alphaHalfWidth_deg','minSamplesPerBin'};
for i = 1:numel(requiredData)
    if ~isfield(donor,requiredData{i})
        error('lc0NDFitLateral:missingField','DONOR must contain %s.',requiredData{i});
    end
end
for i = 1:numel(requiredOpts)
    if ~isfield(opts,requiredOpts{i})
        error('lc0NDFitLateral:missingOption','OPTS must contain %s.',requiredOpts{i});
    end
end

kappa = donor.slipRatio(:);
alpha = donor.slipAngle_deg(:);
fy = donor.Fy_N(:);
fz = donor.Fz_N(:);
pressure = donor.pressure_psi(:);
camber = donor.camber_deg(:);
n = numel(alpha);
if any([numel(kappa),numel(fy),numel(fz),numel(pressure),numel(camber)] ~= n)
    error('lc0NDFitLateral:sizeMismatch','All DONOR channels must have equal length.');
end
if any(fz <= 0)
    error('lc0NDFitLateral:nonpositiveLoad','DONOR Fz_N must be positive.');
end
use = isfinite(kappa) & isfinite(alpha) & isfinite(fy) & isfinite(fz) & ...
    isfinite(pressure) & isfinite(camber) & ...
    abs(kappa) <= opts.maxAbsSlipRatio & ...
    abs(fz-opts.load_N) <= opts.loadHalfWidth_N & ...
    abs(pressure-opts.pressure_psi) <= opts.pressureHalfWidth_psi & ...
    abs(camber-opts.camber_deg) <= opts.camberHalfWidth_deg;
muY = fy./fz;
grid = opts.alphaGrid(:);
mu = nan(numel(grid),1); sigma = nan(numel(grid),1); count = zeros(numel(grid),1);
for i = 1:numel(grid)
    inBin = use & abs(alpha-grid(i)) <= opts.alphaHalfWidth_deg;
    count(i) = sum(inBin);
    if count(i) > 0
        mu(i) = mean(muY(inBin));
        sigma(i) = std(muY(inBin));
    end
end
qualified = count >= opts.minSamplesPerBin;
curve = table(grid,mu,sigma,count,qualified,'VariableNames', ...
    {'slip_angle_deg','mu_y','mu_y_std','n_samples','is_qualified'});
fit = struct('options',opts,'n_selected',sum(use),'curve',curve, ...
    'peak_abs_mu_y',max(abs(mu(qualified)),[],'omitnan'));
end
