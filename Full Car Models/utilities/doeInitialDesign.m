function [U,T] = doeInitialDesign(resolved,n,seed,method)
%DOEINITIALDESIGN Generate a repeatable initial DOE with its baseline case.

if nargin < 4 || isempty(method), method = "LHS"; end
validateattributes(n,{'numeric'},{'scalar','integer','positive'})
validateattributes(seed,{'numeric'},{'scalar','integer','nonnegative'})
method = lower(string(method));
if ~ismember(method,["lhs","random"])
    error('doeInitialDesign:badMethod','method must be LHS or Random.')
end
if ~isstruct(resolved) || ~isfield(resolved,'parameters') || ...
        ~isfield(resolved,'toPhysical') || ~isa(resolved.toPhysical,'function_handle')
    error('doeInitialDesign:badResolved', ...
        'resolved must contain parameters and a toPhysical function handle.')
end

P = resolved.parameters;
required = {'baseline','lowerPhysical','upperPhysical'};
if ~istable(P) || ~all(ismember(required,P.Properties.VariableNames))
    error('doeInitialDesign:badResolved', ...
        'resolved.parameters must provide baseline and physical bounds.')
end
baselineU = (P.baseline-P.lowerPhysical) ./ (P.upperPhysical-P.lowerPhysical);
if any(~isfinite(baselineU) | baselineU < 0 | baselineU > 1)
    error('doeInitialDesign:baselineOutsideBounds', ...
        'Each baseline value must lie within its resolved physical bounds.')
end

previousRng = rng;
restoreRng = onCleanup(@() rng(previousRng)); %#ok<NASGU>
rng(seed,'twister');
nVar = height(P);
U = rand(n,nVar);
if method == "lhs"
    for j = 1:nVar
        U(:,j) = (randperm(n)' - U(:,j))/n;
    end
end

if ~any(all(abs(U-baselineU') <= 1e-12,2))
    [~,nearest] = min(sum((U-baselineU').^2,2));
    U(nearest,:) = baselineU';
end
T = resolved.toPhysical(U);
end
