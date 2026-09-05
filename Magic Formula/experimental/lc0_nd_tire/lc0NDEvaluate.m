function [Fx,Fy,info] = lc0NDEvaluate(model,alpha_deg,kappa,Fz_N)
%LC0NDEVALUATE Evaluate the isolated target/donor LC0 force prototype.
%   Out-of-reference alpha or kappa points return NaN rather than silently
%   extrapolating beyond the measured target or provisional donor curves.

if any(Fz_N(:) <= 0)
    error('lc0NDEvaluate:nonpositiveLoad','Fz_N must be positive.');
end
yCurve = model.target_lateral.curve(model.target_lateral.curve.is_qualified,:);
xCurve = model.donor_longitudinal.curve(model.donor_longitudinal.curve.is_qualified,:);
alpha = alpha_deg(:);
kappa = kappa(:);
fz = Fz_N(:);
if ~(numel(alpha) == numel(kappa) && numel(alpha) == numel(fz))
    error('lc0NDEvaluate:sizeMismatch','alpha_deg, kappa, and Fz_N must match in size.');
end
muY = interp1(yCurve.slip_angle_deg,yCurve.mu_y,alpha,'linear',NaN);
muX = interp1(xCurve.slip_ratio,xCurve.mu_x,kappa*model.rho_stiff,'linear',NaN);
Fx0 = fz.*model.longitudinal_mu_scale.*muX;
Fy0 = fz.*muY;
muXCapacity = model.longitudinal_mu_scale*max(abs(xCurve.mu_x));
muYCapacity = model.target_lateral_capacity;
valid = isfinite(Fx0) & isfinite(Fy0);
Fx = nan(size(fz)); Fy = nan(size(fz)); utilization = nan(size(fz));
if any(valid)
    [Fx(valid),Fy(valid),utilization(valid)] = lc0NDCombinedForce( ...
        Fx0(valid),Fy0(valid),muXCapacity,muYCapacity,fz(valid), ...
        model.coupling_exponent);
end
info = struct('Fx0_N',Fx0,'Fy0_N',Fy0,'mu_x',muX,'mu_y',muY, ...
    'utilization',utilization,'is_supported',valid, ...
    'longitudinal_mu_scale',model.longitudinal_mu_scale);
end
