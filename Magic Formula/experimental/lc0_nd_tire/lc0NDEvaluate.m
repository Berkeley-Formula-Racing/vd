function [Fx,Fy,info] = lc0NDEvaluate(model,alpha_deg,kappa,Fz_N,options)
%LC0NDEVALUATE Evaluate the target/donor nondimensional tire model.
%   OPTIONS.outOfRange is 'nan' by default. Use 'clamp' only at a vehicle
%   boundary where finite forces are required; INFO then reports every
%   clamped query so unsupported TTC regions cannot be hidden.

if nargin < 5 || isempty(options), options = struct(); end
if ~isstruct(options) || ~isscalar(options)
    error('lc0NDEvaluate:badOptions','OPTIONS must be a scalar struct.');
end
policy = getOption(options,'outOfRange',"nan");
policy = lower(string(policy));
if ~ismember(policy,["nan","clamp"])
    error('lc0NDEvaluate:badOptions', ...
        'OPTIONS.outOfRange must be ''nan'' or ''clamp''.');
end
scenarioName = getOption(options,'scenario', ...
    getOption(model,'active_scenario',"nominal"));
[scenario,scenarioName] = resolveScenario(model,scenarioName);

alpha = alpha_deg(:);
kappa = kappa(:);
fz = Fz_N(:);
if ~(numel(alpha) == numel(kappa) && numel(kappa) == numel(fz))
    error('lc0NDEvaluate:sizeMismatch','alpha_deg, kappa, and Fz_N must match in size.');
end
if any(~isfinite(fz)) || any(fz < 0)
    error('lc0NDEvaluate:negativeLoad','Fz_N must be finite and nonnegative.');
end

yCurve = qualifiedCurve(model.target_lateral.curve, ...
    'slip_angle_deg','mu_y','target lateral');
xCurve = qualifiedCurve(model.donor_longitudinal.curve, ...
    'slip_ratio','mu_x','donor longitudinal');
pressurePsi = getOption(options,'pressurePsi', ...
    getScalingReference(model,'pressure_psi',12));
camberDeg = getOption(options,'camberDeg', ...
    getScalingReference(model,'camber_deg',0));
if isscalar(pressurePsi), pressurePsi = repmat(pressurePsi,size(fz)); end
if isscalar(camberDeg), camberDeg = repmat(camberDeg,size(fz)); end
if numel(pressurePsi) ~= numel(fz) || numel(camberDeg) ~= numel(fz)
    error('lc0NDEvaluate:sizeMismatch', ...
        'pressurePsi and camberDeg must be scalar or match Fz_N.');
end
pressurePsi = pressurePsi(:);
camberDeg = camberDeg(:);
[muScale,stiffnessScale,scaleExtrapolated] = targetScales( ...
    model,fz,pressurePsi,camberDeg,policy);
[muY,clampY] = interpolateCurve(yCurve.x,yCurve.y, ...
    alpha.*stiffnessScale,policy);
[muX,clampX] = interpolateCurve(xCurve.x,xCurve.y, ...
    kappa.*scenario.rhoStiff,policy);

baseRhoMu = getOption(model,'rho_mu',1);
baseScale = model.longitudinal_mu_scale;
longitudinalScale = baseScale.*scenario.rhoMu./baseRhoMu;
Fx0 = fz.*longitudinalScale.*muX;
Fy0 = fz.*muScale.*muY;
muXCapacity = ones(size(fz)).*longitudinalScale.*max(abs(xCurve.y));
muYCapacity = model.target_lateral_capacity.*muScale;
valid = isfinite(Fx0) & isfinite(Fy0) & fz > 0;
Fx = zeros(size(fz)); Fy = zeros(size(fz)); utilization = zeros(size(fz));
if any(valid)
    [Fx(valid),Fy(valid),utilization(valid)] = lc0NDCombinedForce( ...
        Fx0(valid),Fy0(valid),muXCapacity(valid),muYCapacity(valid),fz(valid), ...
        scenario.couplingExponent);
end
Fx(~valid & fz > 0) = NaN;
Fy(~valid & fz > 0) = NaN;
utilization(~valid & fz > 0) = NaN;

isExtrapolated = (clampX | clampY | scaleExtrapolated) & fz > 0;
info = struct('Fx0_N',Fx0,'Fy0_N',Fy0,'mu_x',muX,'mu_y',muY, ...
    'utilization',utilization,'is_supported',valid | fz == 0, ...
    'is_extrapolated',isExtrapolated,'clamp_count',sum(isExtrapolated), ...
    'query_count',numel(fz),'scenario',string(scenarioName), ...
    'rho_mu',scenario.rhoMu,'rho_stiff',scenario.rhoStiff, ...
    'coupling_exponent',scenario.couplingExponent, ...
    'longitudinal_mu_scale',longitudinalScale, ...
    'lateral_mu_scale',muScale,'lateral_stiffness_scale',stiffnessScale, ...
    'pressure_psi',pressurePsi,'camber_deg',camberDeg);

Fx = reshape(Fx,size(alpha_deg));
Fy = reshape(Fy,size(alpha_deg));
end

function value = getOption(options,name,defaultValue)
if isstruct(options) && isfield(options,name) && ~isempty(options.(name))
    value = options.(name);
else
    value = defaultValue;
end
end

function [scenario,name] = resolveScenario(model,name)
name = string(name);
if isfield(model,'uncertainty_scenarios') && ~isempty(model.uncertainty_scenarios)
    names = string({model.uncertainty_scenarios.name});
    index = find(names == name,1);
    if isempty(index)
        error('lc0NDEvaluate:unknownScenario', ...
            'Unknown uncertainty scenario %s.',name);
    end
    scenario = model.uncertainty_scenarios(index);
else
    scenario = struct('name',name,'rhoMu',getOption(model,'rho_mu',1), ...
        'rhoStiff',getOption(model,'rho_stiff',1), ...
        'couplingExponent',model.coupling_exponent);
end
end

function curve = qualifiedCurve(tableData,xName,yName,label)
required = {xName,yName,'is_qualified'};
if ~all(ismember(required,tableData.Properties.VariableNames))
    error('lc0NDEvaluate:badCurve','%s curve is missing required fields.',label);
end
use = tableData.is_qualified & isfinite(tableData.(xName)) & ...
    isfinite(tableData.(yName));
x = tableData.(xName)(use); y = tableData.(yName)(use);
[x,order] = sort(x); y = y(order);
[x,uniqueIndex] = unique(x,'stable'); y = y(uniqueIndex);
if numel(x) < 2
    error('lc0NDEvaluate:insufficientCurve', ...
        '%s requires at least two finite qualified points.',label);
end
curve = struct('x',x,'y',y);
end

function [value,isClamped] = interpolateCurve(x,y,query,policy)
isClamped = query < x(1) | query > x(end);
if policy == "clamp"
    evaluationQuery = min(max(query,x(1)),x(end));
    value = interp1(x,y,evaluationQuery,'linear');
else
    value = interp1(x,y,query,'linear',NaN);
end
end

function value = getScalingReference(model,name,defaultValue)
value = defaultValue;
if isfield(model,'target_scaling') && isfield(model.target_scaling,'reference') && ...
        isfield(model.target_scaling.reference,name)
    value = model.target_scaling.reference.(name);
end
end

function [muScale,stiffnessScale,isExtrapolated] = targetScales( ...
        model,fz,pressure,camber,policy)
muScale = ones(size(fz));
stiffnessScale = ones(size(fz));
isExtrapolated = false(size(fz));
if ~isfield(model,'target_scaling') || isempty(model.target_scaling) || ...
        ~isfield(model.target_scaling,'table')
    return
end
tableData = model.target_scaling.table;
required = {'load_center_N','pressure_center_psi','camber_center_deg', ...
    'mu_scale','stiffness_scale','is_qualified'};
if ~all(ismember(required,tableData.Properties.VariableNames))
    error('lc0NDEvaluate:badTargetScaling', ...
        'Target scaling table is missing required fields.');
end
use = tableData.is_qualified & isfinite(tableData.load_center_N) & ...
    isfinite(tableData.pressure_center_psi) & ...
    isfinite(tableData.camber_center_deg) & isfinite(tableData.mu_scale) & ...
    isfinite(tableData.stiffness_scale) & tableData.mu_scale > 0 & ...
    tableData.stiffness_scale > 0;
if ~any(use)
    error('lc0NDEvaluate:badTargetScaling', ...
        'Target scaling table contains no qualified finite bins.');
end
loadValues = tableData.load_center_N(use);
pressureValues = tableData.pressure_center_psi(use);
camberValues = tableData.camber_center_deg(use);
muValues = tableData.mu_scale(use);
stiffValues = tableData.stiffness_scale(use);
loadMin = min(loadValues); loadMax = max(loadValues);
pressureMin = min(pressureValues); pressureMax = max(pressureValues);
camberMin = min(camberValues); camberMax = max(camberValues);
scale = [max(loadMax-loadMin,eps),max(pressureMax-pressureMin,eps), ...
    max(camberMax-camberMin,eps)];
for i = 1:numel(fz)
    outside = fz(i) < loadMin || fz(i) > loadMax || ...
        pressure(i) < pressureMin || pressure(i) > pressureMax || ...
        camber(i) < camberMin || camber(i) > camberMax;
    isExtrapolated(i) = outside;
    qLoad = min(max(fz(i),loadMin),loadMax);
    qPressure = min(max(pressure(i),pressureMin),pressureMax);
    qCamber = min(max(camber(i),camberMin),camberMax);
    distance = ((loadValues-qLoad)./scale(1)).^2 + ...
        ((pressureValues-qPressure)./scale(2)).^2 + ...
        ((camberValues-qCamber)./scale(3)).^2;
    [~,index] = min(distance);
    muScale(i) = muValues(index);
    stiffnessScale(i) = stiffValues(index);
    if policy == "nan" && outside
        muScale(i) = NaN;
        stiffnessScale(i) = NaN;
    end
end
end
