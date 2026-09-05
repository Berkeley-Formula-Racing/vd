function result = run_lc0_nd_pacejka_comparison(cfg,legacyEvaluator,makePlots)
%RUN_LC0_ND_PACEJKA_COMPARISON Compare target/donor prototype with Tire2.
%   Both models are evaluated at the same alpha, kappa, load, and camber.
%   This is a force-level comparison only; neither tyre is connected to a
%   vehicle simulation by this runner.

if nargin < 1 || isempty(cfg), cfg = lc0NDConfig(); end
if nargin < 2 || isempty(legacyEvaluator)
    legacyEvaluator = configuredTire2(cfg.reference.pressure_psi);
end
if nargin < 3 || isempty(makePlots), makePlots = true; end
if ~isfolder(cfg.outputDirectory), mkdir(cfg.outputDirectory); end

[targetData,targetManifest] = lc0NDLoadFreeRolling(cfg);
targetFit = lc0NDFitLateral(targetData,cfg.targetForceFit);
donorLong = lc0NDLoadDonor(cfg.donor);
longFit = lc0NDFitLongitudinal(donorLong,cfg.donorFit);
donorLat = lc0NDLoadDonor(cfg.donorLateral);
donorLatFit = lc0NDFitLateral(donorLat,cfg.couplingLateralFit);
model = lc0NDBuildModel(targetFit,longFit,donorLatFit,cfg.calibration);
cases = comparisonCases(targetFit,longFit,cfg.reference);
comparison = lc0NDCompareForceModels(model,cases,legacyEvaluator);
result = struct('config',cfg,'target_manifest',targetManifest, ...
    'target_fit',targetFit,'donor_longitudinal_fit',longFit, ...
    'donor_lateral_fit',donorLatFit,'model',model, ...
    'comparison',comparison);
save(fullfile(cfg.outputDirectory,'lc0_nd_pacejka_comparison.mat'), ...
    'cfg','targetManifest','targetFit','longFit','donorLatFit','model', ...
    'comparison','-v7.3');
if makePlots, savePlots(comparison,cfg.outputDirectory); end
end

function cases = comparisonCases(targetFit,longFit,reference)
y = targetFit.curve(targetFit.curve.is_qualified,:);
x = longFit.curve(longFit.curve.is_qualified,:);
lat = table(repmat("pure lateral",height(y),1),y.slip_angle_deg, ...
    zeros(height(y),1),repmat(reference.load_N,height(y),1), ...
    repmat(reference.camber_deg,height(y),1), 'VariableNames', ...
    {'case_type','alpha_deg','slip_ratio','Fz_N','camber_deg'});
lon = table(repmat("pure longitudinal",height(x),1),zeros(height(x),1), ...
    x.slip_ratio,repmat(reference.load_N,height(x),1), ...
    repmat(reference.camber_deg,height(x),1), 'VariableNames', ...
    lat.Properties.VariableNames);
alpha = selectNonzero(y.slip_angle_deg,4);
kappa = selectNonzero(x.slip_ratio,4);
[alphaGrid,kappaGrid] = ndgrid(alpha,kappa);
n = numel(alphaGrid);
combo = table(repmat("combined",n,1),alphaGrid(:),kappaGrid(:), ...
    repmat(reference.load_N,n,1),repmat(reference.camber_deg,n,1), ...
    'VariableNames',lat.Properties.VariableNames);
cases = [lat;lon;combo];
end

function values = selectNonzero(values,maxCount)
values = values(values ~= 0);
if numel(values) > maxCount
    values = values(round(linspace(1,numel(values),maxCount)));
end
end

function savePlots(comparison,outputDirectory)
f = figure('Name','LC0 experimental vs Pacejka','Color','w');
tiledlayout(1,3,'TileSpacing','compact');
plotMode(comparison,"pure lateral",'alpha_deg','Fy','Slip angle [deg]');
plotMode(comparison,"pure longitudinal",'slip_ratio','Fx','Slip ratio');
nexttile;
use = comparison.case_type == "combined" & comparison.experimental_supported;
scatter(comparison.experimental_Fx_N(use)./comparison.Fz_N(use), ...
    comparison.experimental_Fy_N(use)./comparison.Fz_N(use),45,'filled', ...
    'DisplayName','Experimental'); hold on;
scatter(comparison.legacy_Fx_N(use)./comparison.Fz_N(use), ...
    comparison.legacy_Fy_N(use)./comparison.Fz_N(use),35,'o', ...
    'DisplayName','Current Pacejka');
xlabel('F_x/F_z'); ylabel('F_y/F_z'); title('Combined-slip locus'); grid on; legend;
exportgraphics(f,fullfile(outputDirectory,'lc0_nd_vs_pacejka.png'),'Resolution',180);
savefig(f,fullfile(outputDirectory,'lc0_nd_vs_pacejka.fig'));
close(f);
end

function plotMode(comparison,mode,xName,forceName,xLabel)
nexttile;
use = comparison.case_type == mode & comparison.experimental_supported;
x = comparison.(xName)(use);
if forceName == "Fx"
    expForce = comparison.experimental_Fx_N(use)./comparison.Fz_N(use);
    legacyForce = comparison.legacy_Fx_N(use)./comparison.Fz_N(use);
else
    expForce = comparison.experimental_Fy_N(use)./comparison.Fz_N(use);
    legacyForce = comparison.legacy_Fy_N(use)./comparison.Fz_N(use);
end
plot(x,expForce,'LineWidth',1.6,'DisplayName','Experimental'); hold on;
plot(x,legacyForce,'--','LineWidth',1.4,'DisplayName','Current Pacejka');
xlabel(xLabel); ylabel([char(forceName) '/F_z']); title(char(mode)); grid on; legend;
end

function evaluator = configuredTire2(pressurePsi)
packageRoot = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(fileparts(packageRoot)));
addpath(fullfile(repoRoot,'Full Car Models'));
setup_paths;
carCell = carConfig();
tire = carCell{1,1}.tire;
tire.p_i = pressurePsi;
evaluator = @(alpha,kappa,Fz,gamma) tireForce(tire,alpha,kappa,Fz,gamma);
end

function [Fx,Fy] = tireForce(tire,alpha,kappa,Fz,gamma)
Fx = tire.F_x(alpha,kappa,Fz,gamma);
Fy = tire.F_y(alpha,kappa,Fz,gamma);
end
