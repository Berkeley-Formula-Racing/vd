function [Fx,Fy,utilization] = lc0NDCombinedForce(Fx0,Fy0,muX,muY,Fz,p)
%LC0NDCOMBINEDFORCE Apply a bounded, anisotropic combined-force envelope.
%   This is a transparent baseline coupling law for the unmeasured LC0
%   longitudinal/combined-slip region.  It preserves feasible pure-slip
%   demands exactly and radially caps overloaded combined demands.  It is
%   not a substitute for combined-slip TTC measurements.

if nargin < 6
    error('lc0NDCombinedForce:notEnoughInputs', ...
        'Fx0, Fy0, muX, muY, Fz and p are required.');
end
if any(Fz(:) <= 0) || any(muX(:) <= 0) || any(muY(:) <= 0) || ...
        ~isscalar(p) || ~isfinite(p) || p < 1
    error('lc0NDCombinedForce:invalidEnvelope', ...
        'Fz, muX and muY must be positive and p must be finite and >= 1.');
end

qx = abs(Fx0)./(muX.*Fz);
qy = abs(Fy0)./(muY.*Fz);
utilization = (qx.^p + qy.^p).^(1/p);
scale = min(1,1./max(utilization,eps));
Fx = scale.*Fx0;
Fy = scale.*Fy0;
end
