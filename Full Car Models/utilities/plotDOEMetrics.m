function figs = plotDOEMetrics(X,M,opts)
%PLOTDOEMETRICS Descriptive DOE plots that do not require a fitted model.
if nargin < 3, opts = struct(); end
catalog = doePlotCatalog();
plots = string(getOr(opts,'plots',catalog.id(strcmp(catalog.producer,'plotDOEMetrics'))));
visible = getOr(opts,'visible','on');
supported = ["quality","pareto_events","pareto_energy", ...
    "correlation","speed_grip","balance","sobol"];
bad = setdiff(plots,supported);
if ~isempty(bad)
    warning('plotDOEMetrics:modelPlots', ...
        'Plot(s) %s require fitted models and are produced by DOE_Fitting.', ...
        strjoin(bad,', '));
    plots = intersect(plots,supported,'stable');
end

figs = gobjects(0);
for id = plots(:)'
    f = figure('Name',"DOE - "+id,'Color','w','Visible',visible);
    figs(end+1) = f; %#ok<AGROW>
    switch id
        case "quality",       plotQuality(M)
        case "pareto_events", plotEventPareto(M)
        case "pareto_energy", plotEnergyPareto(M)
        case "correlation",   plotCorrelation(X,M)
        case "speed_grip",    plotSpeedGrip(M)
        case "balance",       plotBalance(M)
        case "sobol",         plotSobol(getOr(opts,'sobol',struct()),opts)
    end
    setPlotFont(f);
end

function plotSobol(sobol,opts)
responses = string(fieldnames(sobol));
keep = false(size(responses));
for i = 1:numel(responses)
    name = lower(responses(i));
    keep(i) = contains(name,"point") || contains(name,"score") || ...
        startsWith(name,"t_") || contains(name,"work") || ...
        contains(name,"energy") || contains(name,"understeer") || ...
        contains(name,"rebalance");
end
responses = responses(keep);
if isempty(responses)
    text(.5,.5,'No supported Sobol responses are available.', ...
        'HorizontalAlignment','center');
    axis off
    return
end

tl = tiledlayout('flow','TileSpacing','compact');
for i = 1:numel(responses)
    response = responses(i);
    S = sobol.(char(response));
    nexttile
    bar(categorical(S.parameter),[S.firstOrder S.totalOrder],'grouped');
    values = [S.firstOrder;S.totalOrder];
    ylim([min(0,1.05*min(values,[],'omitnan')) ...
        max(1,1.05*max(values,[],'omitnan'))]);
    ylabel('Sobol index');
    legend('First order','Total order','Location','best');
    xtickangle(35); grid on
    title(sobolPanelTitle(response,opts),'Interpreter','none');
end
title(tl,'Surrogate Sobol sensitivity');
end

function label = sobolPanelTitle(response,opts)
family = "quadratic";
nrmse = NaN;
if isfield(opts,'preferredModel') && isfield(opts.preferredModel,char(response))
    family = string(opts.preferredModel.(char(response)));
end
if isfield(opts,'validation') && isfield(opts.validation,char(response))
    V = opts.validation.(char(response));
    if isfield(V,char(family)) && isfield(V.(char(family)),'normalizedRMSE')
        nrmse = V.(char(family)).normalizedRMSE;
    end
end
if isfinite(nrmse)
    label = sprintf('%s (%s, CV NRMSE %.3g)',response,family,nrmse);
else
    label = sprintf('%s (%s)',response,family);
end
end
end

function plotQuality(M)
tiledlayout(2,2,'TileSpacing','compact');
nexttile; bar([sum(M.valid),sum(~M.valid)]); set(gca,'XTickLabel',{'Valid','Invalid'});
ylabel('Cases'); title('Case validity'); grid on
nexttile; histogram(M.gg_coverage(M.valid),10); xlabel('G-G retained coverage'); title('Envelope coverage'); grid on
nexttile; scatter(M.case_index,M.min_Fz_N,28,M.wheel_lift_fraction,'filled'); yline(0,'--');
xlabel('Case'); ylabel('Minimum F_z (N)'); title('Wheel-load margin'); colorbar; grid on
nexttile
eventTimes = [M.t_skid M.t_accel M.t_autox];
eventNames = categorical(repmat(["Skidpad","Accel","Autocross"],height(M),1));
boxchart(eventNames(:),eventTimes(:));
ylabel('Time (s)'); title('Event distributions'); grid on
end

function plotEventPareto(M)
ok=M.valid;
scatter(M.t_autox(ok),M.t_accel(ok),55,M.t_skid(ok),'filled');
xlabel('Autocross (s)'); ylabel('Acceleration (s)'); cb=colorbar; cb.Label.String='Skidpad (s)';
title('Event-time Pareto view'); grid on
end

function plotEnergyPareto(M)
ok=M.valid;
scatter(M.total_work_kJ(ok),M.t_autox(ok),55,M.v_mean_mps(ok),'filled');
xlabel('Autocross tractive work (kJ)'); ylabel('Autocross time (s)');
cb=colorbar; cb.Label.String='Mean speed (m/s)'; title('Performance versus energy'); grid on
end

function plotCorrelation(X,M)
metricNames = {'t_autox','t_accel','t_skid','total_work_kJ','gLat_peak_g', ...
    'gg_lat_20_g','gg_accel_20_g','gg_brake_20_g','min_Fz_N'};
vary = var(X{:,:},0,1,'omitnan') > 1e-12;
Xv=X{:,vary}; xn=string(X.Properties.VariableNames(vary));
Y=M{:,metricNames}; names=[xn string(metricNames)];
C=corrcoef([Xv Y],'Rows','pairwise');
imagesc(C,[-1 1]); axis image; colorbar; colormap(redblue());
set(gca,'XTick',1:numel(names),'YTick',1:numel(names), ...
    'XTickLabel',names,'YTickLabel',names,'TickLabelInterpreter','none');
xtickangle(55); title('Input/output correlation');
end

function plotSpeedGrip(M)
v=[10 20 30];
L=[M.gg_lat_10_g M.gg_lat_20_g M.gg_lat_30_g];
A=[M.gg_accel_10_g M.gg_accel_20_g M.gg_accel_30_g];
B=[M.gg_brake_10_g M.gg_brake_20_g M.gg_brake_30_g];
hold on
errorbar(v,mean(L,1,'omitnan'),std(L,0,1,'omitnan'),'-o','LineWidth',1.5);
errorbar(v,mean(A,1,'omitnan'),std(A,0,1,'omitnan'),'-o','LineWidth',1.5);
errorbar(v,mean(B,1,'omitnan'),std(B,0,1,'omitnan'),'-o','LineWidth',1.5);
xlabel('Vehicle speed (m/s)'); ylabel('Capability (g)');
legend('Lateral','Acceleration','Braking','Location','best'); title('G-G capability versus speed'); grid on
end

function plotBalance(M)
ok=M.valid;
scatter(M.mechanical_balance(ok),M.aero_balance(ok),55, ...
    M.understeer_proxy_25_deg(ok)-M.understeer_proxy_10_deg(ok),'filled');
xlabel('Mechanical balance, LLTD'); ylabel('Aero front fraction');
cb=colorbar; cb.Label.String='High-low speed understeer proxy (deg)';
title('Mechanical and aero balance'); grid on
end

function C = redblue()
n=256; x=linspace(0,1,n/2)';
C=[x x ones(n/2,1); ones(n/2,1) flipud(x) flipud(x)];
end

function value=getOr(s,name,default)
if isfield(s,name), value=s.(name); else, value=default; end
end
