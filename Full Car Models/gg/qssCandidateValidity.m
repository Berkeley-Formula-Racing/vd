function valid = qssCandidateValidity(exitflag,residual,tolerance)
%QSSCANDIDATEVALIDITY Canonical g-g optimizer acceptance contract.
if nargin < 3 || isempty(tolerance), tolerance = 1e-4; end
valid = isscalar(exitflag) && isscalar(residual) && ...
    ismember(exitflag,[1,2]) && isfinite(residual) && ...
    residual <= tolerance;
end
