function A = doeAnalyze(resultPath,opts)
%DOEANALYZE Extract DOE metrics, fit selected responses, and create plots.
if nargin < 2, opts = struct(); end
S = load(resultPath,'carCell','designTable');
if ~isfield(S,'carCell') || ~isfield(S,'designTable')
    error('doeAnalyze:badResults','Result file must contain carCell and designTable.');
end

metricOpts = struct('runRampMetrics',getOr(opts,'runRampMetrics',false));
if isfield(opts,'rampOpts'), metricOpts.rampOpts=opts.rampOpts; end
[metricTable,metricDetails] = doeMetrics(S.carCell,metricOpts);

valid = metricTable.valid & metricTable.t_skid <= 8 & ...
    metricTable.t_autox >= 35 & metricTable.t_autox <= 75 & ...
    metricTable.t_accel >= 3;
quality = struct('total',height(metricTable),'valid',sum(valid), ...
    'invalid',sum(~valid),'wheelLiftCases',sum(metricTable.min_Fz_N < 0), ...
    'meanGGCoverage',mean(metricTable.gg_coverage(valid),'omitnan'));

metricVars = setdiff(metricTable.Properties.VariableNames, ...
    {'case_index','valid','error_message'},'stable');
fullTable = [S.designTable,metricTable(:,metricVars)];
cleanTable = fullTable(valid,:);

allPredictors = S.designTable.Properties.VariableNames;
varying = false(size(allPredictors));
for i=1:numel(allPredictors)
    x=cleanTable.(allPredictors{i});
    varying(i)=(max(x)-min(x)) > 1e-12*max(1,abs(mean(x)));
end
predictors = allPredictors(varying);
requestedPredictors = getOr(opts,'predictors',{});
if ~isempty(requestedPredictors)
    missing = setdiff(requestedPredictors,allPredictors);
    if ~isempty(missing)
        error('doeAnalyze:unknownPredictor','Unknown predictor(s): %s',strjoin(missing,', '));
    end
    predictors = intersect(requestedPredictors,predictors,'stable');
end

defaultResponses = {'t_autox','t_accel','t_skid','total_work_kJ', ...
    'gLat_peak_g','gg_lat_10_g','gg_lat_20_g','gg_lat_30_g', ...
    'gg_accel_20_g','gg_brake_20_g','min_Fz_N', ...
    'understeer_proxy_10_deg','understeer_proxy_25_deg'};
responses = getOr(opts,'responsesWanted',defaultResponses);
responses = intersect(responses,metricVars,'stable');
if isempty(predictors), error('doeAnalyze:noPredictors','No varying predictors remain.'); end
if isempty(responses), error('doeAnalyze:noResponses','No requested response exists in the metric table.'); end

models=struct(); modelOrder=struct();
for i=1:numel(responses)
    response=responses{i};
    vars=[predictors {response}];
    T=cleanTable(:,vars);
    finiteRows=all(isfinite(T{:,:}),2);
    T=T(finiteRows,:);
    p=numel(predictors); n=height(T);
    if n <= p+2
        error('doeAnalyze:tooFewRuns', ...
            '%s has %d valid rows for %d predictors; need at least %d.', ...
            response,n,p,p+3);
    end
    nQuadratic=1+2*p+p*(p-1)/2;
    if n >= nQuadratic+5, upper='quadratic'; else, upper='linear'; end
    models.(response)=stepwiselm(T,'ResponseVar',response, ...
        'Lower','constant','Upper',upper,'Criterion','bic','Verbose',0);
    modelOrder.(response)=upper;
end

catalog=doePlotCatalog();
plots=string(getOr(opts,'plots',catalog.id(catalog.default_enabled)));
visible=getOr(opts,'visible','on');
descriptive=["quality","pareto_events","pareto_energy", ...
    "correlation","speed_grip","balance"];
figures=plotDOEMetrics(S.designTable,metricTable, ...
    struct('plots',intersect(plots,descriptive,'stable'),'visible',visible));

if any(plots=="sensitivity")
    figures(end+1)=plotSensitivity(models,cleanTable,responses,visible); %#ok<AGROW>
end
if any(plots=="main_effects")
    figures(end+1)=plotMainEffects(models,cleanTable,responses,visible); %#ok<AGROW>
end
if any(plots=="validation")
    figures(end+1)=plotValidation(models,cleanTable,responses,visible); %#ok<AGROW>
end
if any(plots=="interactions")
    figures=[figures plotInteractions(models,cleanTable,responses,visible)]; %#ok<AGROW>
end

A=struct('metricTable',metricTable,'metricDetails',metricDetails, ...
    'cleanTable',cleanTable,'models',models,'modelOrder',modelOrder, ...
    'graphCatalog',catalog,'figures',figures,'quality',quality, ...
    'settings',struct('predictors',{predictors},'responses',{responses}, ...
    'plots',plots,'resultPath',string(resultPath)));

savePath=string(getOr(opts,'savePath',""));
if strlength(savePath)>0
    analysis=A; analysis.figures=gobjects(0); %#ok<NASGU>
    save(savePath,'analysis','-v7.3');
end
end

function f=plotSensitivity(models,T,responses,visible)
f=figure('Name','DOE - Ranked sensitivities','Color','w','Visible',visible);
tl=tiledlayout('flow','TileSpacing','compact');
for i=1:numel(responses)
    response=responses{i}; mdl=models.(response);
    names=mdl.CoefficientNames(2:end); beta=mdl.Coefficients.Estimate(2:end);
    impact=zeros(numel(names),1);
    sy=std(T.(response),'omitnan');
    for j=1:numel(names)
        z=termValue(T,names{j});
        impact(j)=beta(j)*std(z,'omitnan')/max(sy,eps);
    end
    [~,order]=sort(abs(impact),'descend'); order=order(1:min(8,numel(order)));
    nexttile; barh(impact(order)); set(gca,'YTick',1:numel(order), ...
        'YTickLabel',names(order),'TickLabelInterpreter','none','YDir','reverse');
    xlabel('Standardized effect'); title(response,'Interpreter','none'); grid on
end
title(tl,'Ranked normalized sensitivity'); setPlotFont(f);
end

function f=plotMainEffects(models,T,responses,visible)
f=figure('Name','DOE - Main effects','Color','w','Visible',visible);
tl=tiledlayout('flow','TileSpacing','compact');
for i=1:numel(responses)
    response=responses{i}; mdl=models.(response);
    p=mdl.PredictorNames;
    nexttile
    if isempty(p), text(.5,.5,'Constant model','HorizontalAlignment','center'); axis off; continue, end
    xName=p{1}; x=linspace(min(T.(xName)),max(T.(xName)),60)';
    Q=meanPredictorTable(T,p,numel(x)); Q.(xName)=x;
    plot(x,predict(mdl,Q),'LineWidth',1.6); hold on; scatter(T.(xName),T.(response),15,'filled');
    xlabel(xName,'Interpreter','none'); ylabel(response,'Interpreter','none'); grid on
end
title(tl,'Main effect of first retained predictor'); setPlotFont(f);
end

function f=plotValidation(models,T,responses,visible)
f=figure('Name','DOE - Surrogate validation','Color','w','Visible',visible);
tl=tiledlayout('flow','TileSpacing','compact');
for i=1:numel(responses)
    response=responses{i}; mdl=models.(response); y=T.(response); yp=predict(mdl,T);
    nexttile; scatter(y,yp,24,'filled'); hold on
    lim=[min([y;yp]) max([y;yp])]; plot(lim,lim,'k--'); axis equal; xlim(lim); ylim(lim);
    xlabel('Simulated'); ylabel('Predicted'); title(response,'Interpreter','none'); grid on
end
title(tl,'Predicted versus simulated (training data)'); setPlotFont(f);
end

function figs=plotInteractions(models,T,responses,visible)
figs=gobjects(0);
for i=1:numel(responses)
    response=responses{i}; mdl=models.(response); p=mdl.PredictorNames;
    if numel(p)<2, continue, end
    x1=linspace(min(T.(p{1})),max(T.(p{1})),35);
    x2=linspace(min(T.(p{2})),max(T.(p{2})),35);
    [X1,X2]=meshgrid(x1,x2); Q=meanPredictorTable(T,p,numel(X1));
    Q.(p{1})=X1(:); Q.(p{2})=X2(:); Z=reshape(predict(mdl,Q),size(X1));
    f=figure('Name',['DOE interaction - ' response],'Color','w','Visible',visible);
    surf(X1,X2,Z,'EdgeColor','none'); xlabel(p{1},'Interpreter','none');
    ylabel(p{2},'Interpreter','none'); zlabel(response,'Interpreter','none'); colorbar; grid on
    setPlotFont(f); figs(end+1)=f; %#ok<AGROW>
end
end

function Q=meanPredictorTable(T,names,n)
Q=table();
for i=1:numel(names), Q.(names{i})=repmat(mean(T.(names{i}),'omitnan'),n,1); end
end

function z=termValue(T,name)
parts=strsplit(name,':'); z=ones(height(T),1);
for i=1:numel(parts)
    base=regexprep(parts{i},'\^2$',''); q=T.(base);
    if endsWith(parts{i},'^2'), q=q.^2; end
    z=z.*q;
end
end

function value=getOr(s,name,default)
if isfield(s,name), value=s.(name); else, value=default; end
end
