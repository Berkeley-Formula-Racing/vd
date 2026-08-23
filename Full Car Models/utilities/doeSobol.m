function S = doeSobol(analysis,response,nSamples,seed)
%DOESOBOL Estimate surrogate Sobol first- and total-order sensitivities.
if nargin < 3 || isempty(nSamples), nSamples = 4096; end
if nargin < 4 || isempty(seed), seed = 0; end
if ~isscalar(nSamples) || nSamples < 2 || nSamples ~= floor(nSamples)
    error('doeSobol:badSampleCount','nSamples must be an integer of at least two.');
end
if ~isscalar(seed) || ~isfinite(seed)
    error('doeSobol:badSeed','seed must be a finite scalar.');
end
if ~isfield(analysis,'predictorBounds') || ~istable(analysis.predictorBounds)
    error('doeSobol:missingBounds','analysis.predictorBounds is required.');
end

response = char(string(response));
bounds = analysis.predictorBounds;
required = {'parameter','lower','upper'};
if ~all(ismember(required,bounds.Properties.VariableNames))
    error('doeSobol:badBounds','Predictor bounds require parameter, lower, and upper columns.');
end
parameter = string(bounds.parameter(:));
lower = bounds.lower(:);
upper = bounds.upper(:);
if isempty(parameter) || any(~isfinite(lower) | ~isfinite(upper) | upper <= lower)
    error('doeSobol:badBounds','Every predictor requires finite increasing bounds.');
end

originalRng = rng;
cleanup = onCleanup(@() rng(originalRng)); %#ok<NASGU>
rng(double(seed),'twister');
p = numel(parameter);
A = rand(nSamples,p);
B = rand(nSamples,p);
fA = predictPhysical(analysis,response,A,parameter,lower,upper);
fB = predictPhysical(analysis,response,B,parameter,lower,upper);
fAB = zeros(nSamples,p);
for i = 1:p
    ABi = A;
    ABi(:,i) = B(:,i);
    fAB(:,i) = predictPhysical(analysis,response,ABi,parameter,lower,upper);
end

[firstOrder,totalOrder,materialExcursion] = sobolIndices(fA,fB,fAB);
nBootstrap = bootstrapCount(analysis);
firstBootstrap = zeros(nBootstrap,p);
totalBootstrap = zeros(nBootstrap,p);
for b = 1:nBootstrap
    rows = randi(nSamples,nSamples,1);
    [firstBootstrap(b,:),totalBootstrap(b,:)] = sobolIndices( ...
        fA(rows),fB(rows),fAB(rows,:));
end

S = table(parameter,firstOrder,totalOrder, ...
    prctile(firstBootstrap,2.5,1)',prctile(firstBootstrap,97.5,1)', ...
    prctile(totalBootstrap,2.5,1)',prctile(totalBootstrap,97.5,1)', ...
    materialExcursion, ...
    'VariableNames',{'parameter','firstOrder','totalOrder', ...
    'firstOrderLow','firstOrderHigh','totalOrderLow','totalOrderHigh', ...
    'materialExcursion'});
end

function y = predictPhysical(analysis,response,U,parameter,lower,upper)
X = lower' + U .* (upper-lower)';
T = array2table(X,'VariableNames',cellstr(parameter));
y = predictDOEModel(analysis,response,T);
if any(~isfinite(y))
    error('doeSobol:nonfinitePrediction', ...
        'Selected surrogate returned a nonfinite prediction for %s.',response);
end
y = y(:);
end

function [firstOrder,totalOrder,material] = sobolIndices(fA,fB,fAB)
V = var([fA;fB],1);
if ~isfinite(V) || V <= eps(max(1,max(abs([fA;fB]))))
    error('doeSobol:constantResponse', ...
        'Selected surrogate has insufficient prediction variance for Sobol indices.');
end
firstRaw = 1 - mean((fB-fAB).^2,1)'/(2*V);
totalRaw = mean((fA-fAB).^2,1)'/(2*V);
tiny = 1e-10;
material = firstRaw < -tiny | firstRaw > 1+tiny | ...
    totalRaw < -tiny | totalRaw > 1+tiny;
firstOrder = clampTiny(firstRaw,tiny);
totalOrder = clampTiny(totalRaw,tiny);
end

function value = clampTiny(value,tiny)
value(value < 0 & value >= -tiny) = 0;
value(value > 1 & value <= 1+tiny) = 1;
end

function n = bootstrapCount(analysis)
n = 200;
if isfield(analysis,'settings') && isfield(analysis.settings,'sobolBootstrapSamples')
    n = analysis.settings.sobolBootstrapSamples;
end
if ~isscalar(n) || n < 1 || n ~= floor(n)
    error('doeSobol:badBootstrapCount', ...
        'sobolBootstrapSamples must be a positive integer.');
end
end
