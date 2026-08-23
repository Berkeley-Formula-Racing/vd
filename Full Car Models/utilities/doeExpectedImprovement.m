function ei = doeExpectedImprovement(mu,sigma,best)
%DOEEXPECTEDIMPROVEMENT Expected improvement for a maximization objective.

validateattributes(mu,{'numeric'},{'real','finite'})
validateattributes(sigma,{'numeric'},{'real','finite','nonnegative'})
validateattributes(best,{'numeric'},{'real','finite','scalar'})
if ~isequal(size(mu),size(sigma))
    error('doeExpectedImprovement:sizeMismatch', ...
        'mu and sigma must have the same size.')
end

z = (mu-best) ./ max(sigma,eps);
ei = (mu-best) .* normcdf(z) + sigma .* normpdf(z);
zeroSigma = sigma <= eps;
ei(zeroSigma) = max(mu(zeroSigma)-best,0);
end
