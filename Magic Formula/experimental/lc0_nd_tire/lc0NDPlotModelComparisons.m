function artifacts = lc0NDPlotModelComparisons(resultOrFile,outputDirectory,options)
%LC0NDPLOTMODELCOMPARISONS Plot matched tire-model force comparisons.
%   Consumes the saved result from run_lc0_nd_pacejka_comparison (or the
%   result struct itself). It does not run a vehicle, DOE, or TTC fitting
%   simulation.

if nargin < 1 || isempty(resultOrFile)
    error('lc0NDPlotModelComparisons:missingResult', ...
        'Pass a comparison result struct or a saved MAT file.');
end
if nargin < 2, outputDirectory = []; end
if nargin < 3 || isempty(options), options = struct(); end
if ~isstruct(options) || ~isscalar(options)
    error('lc0NDPlotModelComparisons:badOptions', ...
        'OPTIONS must be a scalar struct.');
end

result = loadResult(resultOrFile);
validateResult(result);
result.comparison = normalizeComparison(result.comparison);
if isempty(outputDirectory)
    outputDirectory = inferOutputDirectory(result);
end
if ~(ischar(outputDirectory) || ...
        (isstring(outputDirectory) && isscalar(outputDirectory)))
    error('lc0NDPlotModelComparisons:badOutputDirectory', ...
        'OUTPUTDIRECTORY must be a path string.');
end
outputDirectory = char(outputDirectory);
if ~isfolder(outputDirectory), mkdir(outputDirectory); end

fileStem = char(getOption(options,'fileStem','lc0_nd_model_comparisons'));
visible = char(getOption(options,'visible','off'));
closeFigure = logical(getOption(options,'closeFigure',true));

comparison = result.comparison;
caseType = string(comparison.case_type);
supported = logical(comparison.experimental_supported);
finiteExperimental = isfinite(comparison.experimental_Fx_N) & ...
    isfinite(comparison.experimental_Fy_N) & isfinite(comparison.Fz_N) & ...
    comparison.Fz_N > 0;
finiteLegacy = isfinite(comparison.legacy_Fx_N) & ...
    isfinite(comparison.legacy_Fy_N);
pureLateral = caseType == "pure lateral" & supported & finiteExperimental & ...
    finiteLegacy & isfinite(comparison.alpha_deg);
pureLongitudinal = caseType == "pure longitudinal" & supported & ...
    finiteExperimental & finiteLegacy & isfinite(comparison.slip_ratio);
combined = caseType == "combined" & finiteExperimental & finiteLegacy;

metrics = buildMetrics(comparison,pureLateral,pureLongitudinal,combined);
summary = buildSummary(comparison,caseType,pureLateral,pureLongitudinal,combined);

f = figure('Name','LC0 tire-model comparisons','Color','w','Visible',visible);
layout = tiledlayout(f,2,4,'TileSpacing','compact','Padding','compact');
plotPureLateral(layout,comparison,pureLateral);
plotPureLongitudinal(layout,comparison,pureLongitudinal);
envelopes = plotCombinedEnvelope(layout,result,comparison,combined);
plotUtilization(layout,comparison,combined);
plotLateralError(layout,comparison,pureLateral);
plotLongitudinalError(layout,comparison,pureLongitudinal);
plotCombinedError(layout,comparison,combined);
plotMetrics(layout,metrics,summary);
title(layout,'Matched tire-model comparison; normalized forces use F_z');

pngFile = fullfile(outputDirectory,[fileStem '.png']);
figFile = fullfile(outputDirectory,[fileStem '.fig']);
csvFile = fullfile(outputDirectory,[fileStem '_metrics.csv']);
exportgraphics(f,pngFile,'Resolution',180);
savefig(f,figFile);
writetable(metrics,csvFile);

artifacts = struct('pngFile',pngFile,'figFile',figFile,'csvFile',csvFile, ...
    'metrics',metrics,'summary',summary,'envelopes',envelopes, ...
    'figureHandle',f);
if closeFigure
    close(f);
    artifacts.figureHandle = [];
end
end

function result = loadResult(resultOrFile)
if isstruct(resultOrFile)
    result = resultOrFile;
elseif ischar(resultOrFile) || ...
        (isstring(resultOrFile) && isscalar(resultOrFile))
    file = char(resultOrFile);
    if ~isfile(file)
        error('lc0NDPlotModelComparisons:missingFile', ...
            'Result file does not exist: %s',file);
    end
    loaded = load(file);
    if isfield(loaded,'result') && isstruct(loaded.result)
        result = loaded.result;
    else
        result = loaded;
    end
else
    error('lc0NDPlotModelComparisons:badResult', ...
        'RESULTORFILE must be a result struct or a MAT-file path.');
end
end

function validateResult(result)
requiredResult = {'model','comparison'};
if ~isstruct(result) || ~all(isfield(result,requiredResult)) || ...
        ~istable(result.comparison)
    error('lc0NDPlotModelComparisons:badResult', ...
        'The result must contain a model and comparison table.');
end
requiredColumns = {'case_type','alpha_deg','slip_ratio','Fz_N', ...
    'experimental_Fx_N','experimental_Fy_N','legacy_Fx_N','legacy_Fy_N', ...
    'experimental_supported','experimental_utilization'};
if ~all(ismember(requiredColumns,result.comparison.Properties.VariableNames))
    error('lc0NDPlotModelComparisons:badResult', ...
        'The comparison table is missing required columns.');
end
requiredModel = {'target_lateral','donor_longitudinal', ...
    'target_lateral_capacity','longitudinal_mu_scale'};
if ~all(isfield(result.model,requiredModel))
    error('lc0NDPlotModelComparisons:badResult', ...
        'The result model is missing capacity information.');
end
end

function comparison = normalizeComparison(comparison)
if ~ismember('experimental_extrapolated',comparison.Properties.VariableNames)
    comparison.experimental_extrapolated = false(height(comparison),1);
end
end

function value = inferOutputDirectory(result)
value = fullfile(pwd,'results');
if isfield(result,'config') && isfield(result.config,'outputDirectory') && ...
        ~isempty(result.config.outputDirectory)
    value = result.config.outputDirectory;
end
end

function value = getOption(options,name,defaultValue)
if isfield(options,name) && ~isempty(options.(name))
    value = options.(name);
else
    value = defaultValue;
end
end

function metrics = buildMetrics(comparison,pureLateral,pureLongitudinal,combined)
fxUse = (pureLongitudinal | combined) & ...
    isfinite(comparison.experimental_Fx_N) & ...
    isfinite(comparison.legacy_Fx_N);
fyUse = (pureLateral | combined) & ...
    isfinite(comparison.experimental_Fy_N) & ...
    isfinite(comparison.legacy_Fy_N);
metrics = [oneMetric("Fx",comparison.legacy_Fx_N(fxUse), ...
    comparison.experimental_Fx_N(fxUse),comparison.Fz_N(fxUse)); ...
    oneMetric("Fy",comparison.legacy_Fy_N(fyUse), ...
    comparison.experimental_Fy_N(fyUse),comparison.Fz_N(fyUse))];
end

function row = oneMetric(quantity,legacy,experimental,fz)
errorN = legacy - experimental;
errorMu = errorN./fz;
row = table(quantity,numel(errorN),mean(abs(errorN)), ...
    sqrt(mean(errorN.^2)),max(abs(errorN)),mean(abs(errorMu)), ...
    sqrt(mean(errorMu.^2)),max(abs(errorMu)), ...
    'VariableNames',{'quantity','sample_count','mae_N','rmse_N', ...
    'max_abs_N','mae_mu','rmse_mu','max_abs_mu'});
end

function summary = buildSummary(comparison,caseType,pureLateral,pureLongitudinal,combined)
summary = struct('pure_lateral_count',sum(pureLateral), ...
    'pure_longitudinal_count',sum(pureLongitudinal), ...
    'combined_count',sum(caseType == "combined"), ...
    'combined_supported_count',sum(combined & ...
    comparison.experimental_supported), ...
    'combined_extrapolated_count',sum(caseType == "combined" & ...
    comparison.experimental_extrapolated));
end

function plotPureLateral(layout,comparison,use)
nexttile(layout); hold on; grid on; box on;
if ~any(use), noData('Pure lateral'); return; end
indices = find(use);
[x,order] = sort(comparison.alpha_deg(indices));
indices = indices(order);
plot(x,comparison.experimental_Fy_N(indices)./comparison.Fz_N(indices), ...
    'LineWidth',1.6,'DisplayName','Nondimensional model');
plot(x,comparison.legacy_Fy_N(indices)./comparison.Fz_N(indices),'--', ...
    'LineWidth',1.4,'DisplayName','Current Tire2');
xlabel('\alpha [deg]'); ylabel('F_y/F_z'); title('Pure lateral');
legend('Location','best'); yline(0,':k','HandleVisibility','off');
end

function plotPureLongitudinal(layout,comparison,use)
nexttile(layout); hold on; grid on; box on;
if ~any(use), noData('Pure longitudinal'); return; end
indices = find(use);
[x,order] = sort(comparison.slip_ratio(indices));
indices = indices(order);
plot(x,comparison.experimental_Fx_N(indices)./comparison.Fz_N(indices), ...
    'LineWidth',1.6,'DisplayName','Nondimensional model');
plot(x,comparison.legacy_Fx_N(indices)./comparison.Fz_N(indices),'--', ...
    'LineWidth',1.4,'DisplayName','Current Tire2');
xlabel('\kappa'); ylabel('F_x/F_z'); title('Pure longitudinal');
legend('Location','best'); yline(0,':k','HandleVisibility','off');
end

function envelopes = plotCombinedEnvelope(layout,result,comparison,use)
nexttile(layout); hold on; grid on; box on;
if any(use)
    scatter(comparison.legacy_Fx_N(use)./comparison.Fz_N(use), ...
        comparison.legacy_Fy_N(use)./comparison.Fz_N(use),35,'o', ...
        'DisplayName','Current Tire2');
    scatter(comparison.experimental_Fx_N(use)./comparison.Fz_N(use), ...
        comparison.experimental_Fy_N(use)./comparison.Fz_N(use),45,'filled', ...
        'DisplayName','Nondimensional model');
end
reference = referenceState(result,comparison,use);
scenarios = scenarioList(result.model);
colors = lines(max(numel(scenarios),1));
scenarioEnvelopes = repmat(struct('name','','mu_x',NaN,'mu_y',NaN, ...
    'p',NaN,'x',[],'y',[]),1,numel(scenarios));
for i = 1:numel(scenarios)
    scenarioName = string(scenarios(i).name);
    [muX,muY,p] = modelCapacities(result.model,reference,scenarioName);
    [x,y] = envelope(muX,muY,p);
    scenarioEnvelopes(i) = struct('name',char(scenarioName),'mu_x',muX, ...
        'mu_y',muY,'p',p,'x',x,'y',y);
    if scenarioName == "nominal"
        lineStyle = '-'; lineWidth = 1.8;
    else
        lineStyle = '--'; lineWidth = 1.2;
    end
    plot(x,y,lineStyle,'Color',colors(i,:),'LineWidth',lineWidth, ...
        'DisplayName',sprintf('%s p-norm envelope',char(scenarioName)));
end
[muX,muY] = modelCapacities(result.model,reference,"nominal");
[x,y] = envelope(muX,muY,2);
plot(x,y,':k','LineWidth',1.4,'DisplayName','Classical p = 2 ellipse');
xlabel('F_x/F_z'); ylabel('F_y/F_z'); title('Combined-slip envelope');
axis equal; legend('Location','best');
envelopes = struct('scenarios',scenarioEnvelopes,'classical_p',2, ...
    'classical_x',x,'classical_y',y);
end

function plotUtilization(layout,comparison,use)
nexttile(layout); hold on; grid on; box on;
allCombined = string(comparison.case_type) == "combined";
if ~any(allCombined)
    noData('Combined utilization'); return;
end
z = comparison.experimental_utilization(allCombined);
scatter(comparison.slip_ratio(allCombined),comparison.alpha_deg(allCombined), ...
    55,z,'filled','DisplayName','Utilization');
finiteZ = z(isfinite(z));
if isempty(finiteZ), upper = 1; else, upper = max(1,max(finiteZ)); end
caxis([0 upper]); colorbar;
unsupported = allCombined & (~comparison.experimental_supported | ...
    comparison.experimental_extrapolated);
if any(unsupported)
    plot(comparison.slip_ratio(unsupported),comparison.alpha_deg(unsupported), ...
        'xk','LineWidth',1.5,'MarkerSize',8,'DisplayName','Unsupported/extrapolated');
end
xlabel('\kappa'); ylabel('\alpha [deg]'); title('Combined support/utilization');
legend('Location','best');
end

function plotLateralError(layout,comparison,use)
nexttile(layout); hold on; grid on; box on;
if ~any(use), noData('Lateral error'); return; end
indices = find(use);
[x,order] = sort(comparison.alpha_deg(indices));
indices = indices(order);
y = (comparison.legacy_Fy_N(indices) - comparison.experimental_Fy_N(indices))./ ...
    comparison.Fz_N(indices);
plot(x,y,'o-','LineWidth',1.2);
yline(0,':k','HandleVisibility','off');
xlabel('\alpha [deg]'); ylabel('(Tire2 - model)/F_z');
title('Pure lateral error');
end

function plotLongitudinalError(layout,comparison,use)
nexttile(layout); hold on; grid on; box on;
if ~any(use), noData('Longitudinal error'); return; end
indices = find(use);
[x,order] = sort(comparison.slip_ratio(indices));
indices = indices(order);
y = (comparison.legacy_Fx_N(indices) - comparison.experimental_Fx_N(indices))./ ...
    comparison.Fz_N(indices);
plot(x,y,'s-','LineWidth',1.2);
yline(0,':k','HandleVisibility','off');
xlabel('\kappa'); ylabel('(Tire2 - model)/F_z');
title('Pure longitudinal error');
end

function plotCombinedError(layout,comparison,use)
nexttile(layout); hold on; grid on; box on;
if ~any(use), noData('Combined error'); return; end
dx = (comparison.legacy_Fx_N(use) - comparison.experimental_Fx_N(use))./ ...
    comparison.Fz_N(use);
dy = (comparison.legacy_Fy_N(use) - comparison.experimental_Fy_N(use))./ ...
    comparison.Fz_N(use);
scatter(dx,dy,50,'filled'); xline(0,':k'); yline(0,':k');
xlabel('\Delta F_x/F_z'); ylabel('\Delta F_y/F_z');
title('Combined force error: Tire2 - model');
end

function plotMetrics(layout,metrics,summary)
nexttile(layout); hold on; grid on; box on;
bar(categorical(metrics.quantity),metrics.rmse_mu);
ylabel('RMSE of normalized force'); title('Error summary');
text(0.03,0.96,sprintf('combined: %d (%d supported)', ...
    summary.combined_count,summary.combined_supported_count), ...
    'Units','normalized','VerticalAlignment','top');
text(0.03,0.86,sprintf('extrapolated: %d',summary.combined_extrapolated_count), ...
    'Units','normalized','VerticalAlignment','top');
end

function noData(label)
axis off;
text(0.5,0.5,['No qualified data: ' label], ...
    'HorizontalAlignment','center');
end

function reference = referenceState(result,comparison,use)
if isfield(result,'config') && isfield(result.config,'reference')
    reference = result.config.reference;
else
    reference = struct();
end
if ~isfield(reference,'load_N') || isempty(reference.load_N)
    reference.load_N = median(comparison.Fz_N(use));
end
if ~isfield(reference,'pressure_psi') || isempty(reference.pressure_psi)
    reference.pressure_psi = 12;
end
if ~isfield(reference,'camber_deg') || isempty(reference.camber_deg)
    reference.camber_deg = 0;
end
end

function scenarios = scenarioList(model)
if isfield(model,'uncertainty_scenarios') && ~isempty(model.uncertainty_scenarios)
    scenarios = model.uncertainty_scenarios;
else
    scenarios = struct('name',"nominal",'rhoMu',1,'rhoStiff',1, ...
        'couplingExponent',model.coupling_exponent);
end
end

function [muX,muY,p] = modelCapacities(model,reference,scenarioName)
options = struct('outOfRange',"clamp",'pressurePsi',reference.pressure_psi, ...
    'camberDeg',reference.camber_deg,'scenario',scenarioName);
[~,~,info] = lc0NDEvaluate(model,0,0,reference.load_N,options);
qualified = model.donor_longitudinal.curve.is_qualified & ...
    isfinite(model.donor_longitudinal.curve.mu_x);
donorPeak = max(abs(model.donor_longitudinal.curve.mu_x(qualified)));
muX = donorPeak.*info.longitudinal_mu_scale;
muY = model.target_lateral_capacity.*info.lateral_mu_scale;
p = info.coupling_exponent;
end

function [x,y] = envelope(muX,muY,p)
theta = linspace(0,2*pi,361)';
x = muX.*sign(cos(theta)).*abs(cos(theta)).^(2/p);
y = muY.*sign(sin(theta)).*abs(sin(theta)).^(2/p);
end
