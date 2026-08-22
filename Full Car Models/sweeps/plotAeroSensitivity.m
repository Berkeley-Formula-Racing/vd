function plotAeroSensitivity(S,param,mode,vSlices)
% Plots aero sensitivity over the four operating axes: vCar, gLat, gLong and
% theta (the g-g polar angle).
%
% inputs  S:        output of aeroSensitivity
%         param:    'ClA' | 'CdA' | 'CoP'  (default: first one available)
%         mode:     'both' (default) | 'continuous' | 'discrete'
%                     continuous - interpolated surfaces and smooth lines
%                     discrete   - markers/stems at the solved grid points
%                                  only, nothing drawn between them
%         vSlices:  velocities (m/s) for the line plots (default: 3 spread
%                   across the solved range)
%
% Figures:
%   1  d(envelope radius)/dp over (vCar, theta)      <- theta axis
%   2  d(gLong)/dp over (vCar, gLat), accel + brake  <- gLat axis
%   3  d(gLat)/dp  over (vCar, gLong)                <- gLong axis
%   4  peak capability sensitivity vs vCar           <- vCar axis
%   5  lap time sensitivity bar chart                (if times were supplied)

if nargin < 2 || isempty(param)
    f = fieldnames(S.sens); param = f{1};
end
if nargin < 3 || isempty(mode), mode = 'both'; end
if ~isfield(S.sens,param)
    error('plotAeroSensitivity:noParam','no sensitivity computed for %s',param);
end

P  = S.sens.(param);
g  = S.grids;
v  = g.vCar;
if nargin < 4 || isempty(vSlices)
    finite = v(any(isfinite(P.r_theta_dp),2));
    if isempty(finite), finite = v; end
    vSlices = finite(round(linspace(1,numel(finite),min(3,numel(finite)))));
end

wantC = any(strcmp(mode,{'both','continuous'}));
wantD = any(strcmp(mode,{'both','discrete'}));
unit  = unitLabel(param);

%% ---- Figure 1: theta axis ----
figure('Name',sprintf('%s sensitivity: envelope vs theta',param), ...
       'Position',[80 80 1100 460]);
setPlotFont(gcf);   % must precede the subplots: these are DEFAULT properties
if wantC
    subplot(1,1+wantD,1);
    % blank out cells where the perturbation levels disagree on the sign of
    % the slope -- those are different local optima, not a gradient
    Z = P.r_theta_dp; Z(~P.r_theta_reliable) = NaN;
    surfaceOrMsg(g.theta,v,Z, ...
        '\theta (deg)  [0 accel, 90 lateral, 180 braking]','vCar (m/s)', ...
        sprintf('d(g-g radius)/d%s   [%s]  continuous',param,unit));
end
if wantD
    subplot(1,1+wantC,1+wantC); hold on
    for k = 1:numel(vSlices)
        i = find(abs(v-vSlices(k))<1e-6,1);
        if isempty(i), continue, end
        y = P.r_theta_dp(i,:);
        plot(g.theta(isfinite(y)),y(isfinite(y)),'o','MarkerSize',4, ...
            'DisplayName',sprintf('%.0f m/s',v(i)));
    end
    grid on; box on; yline(0,'k:','HandleVisibility','off');
    xlabel('\theta (deg)'); ylabel(sprintf('d(radius)/d%s  [%s]',param,unit));
    title('discrete (solved points only)'); legend('Location','best');
end

%% ---- Figure 2: gLat axis ----
figure('Name',sprintf('%s sensitivity: gLong vs gLat',param), ...
       'Position',[110 110 1100 460]);
setPlotFont(gcf);
subplot(1,2,1);
plotSlices(g.gLat,v,P.accel_at_gLat_dp,P.accel_at_gLat_reliable,vSlices,wantC,wantD, ...
    'gLat (g)',sprintf('d(gLong accel)/d%s  [%s]',param,unit), ...
    'Acceleration capability at fixed cornering');
subplot(1,2,2);
plotSlices(g.gLat,v,P.brake_at_gLat_dp,P.brake_at_gLat_reliable,vSlices,wantC,wantD, ...
    'gLat (g)',sprintf('d(gLong brake)/d%s  [%s]',param,unit), ...
    'Braking capability at fixed cornering');

%% ---- Figure 3: gLong axis ----
figure('Name',sprintf('%s sensitivity: gLat vs gLong',param), ...
       'Position',[140 140 700 460]);
setPlotFont(gcf);
plotSlices(g.gLong,v,P.gLat_at_gLong_dp,P.gLat_at_gLong_reliable,vSlices,wantC,wantD, ...
    'gLong (g)  [negative = braking]', ...
    sprintf('d(gLat)/d%s  [%s]',param,unit), ...
    'Cornering capability at fixed longitudinal g');

%% ---- Figure 4: vCar axis ----
figure('Name',sprintf('%s sensitivity: peaks vs vCar',param), ...
       'Position',[170 170 760 460]);
setPlotFont(gcf); hold on
series = {'maxGLat_dp','max lateral g'; ...
          'maxGLong_dp','max acceleration g'; ...
          'maxGBrake_dp','max braking g'};
co = lines(size(series,1));
anyBad = false;
for s = 1:size(series,1)
    y   = P.(series{s,1});
    rel = P.([erase(series{s,1},'_dp') '_reliable']);
    ok  = isfinite(y) & rel;
    bad = isfinite(y) & ~rel;
    yPlot = y; yPlot(~ok) = NaN;   % break the line across bad points
    if wantC
        plot(v,yPlot,'-','Color',co(s,:),'LineWidth',2, ...
            'DisplayName',series{s,2});
    end
    if wantD
        plot(v(ok),y(ok),'o','Color',co(s,:),'MarkerFaceColor',co(s,:), ...
            'MarkerSize',4,'HandleVisibility',onlyIf(~wantC),...
            'DisplayName',series{s,2});
    end
    if any(bad)
        anyBad = true;
        plot(v(bad),y(bad),'x','Color',[0.6 0.6 0.6],'MarkerSize',7, ...
            'LineWidth',1.2,'HandleVisibility','off');
    end
end
grid on; box on; yline(0,'k:','HandleVisibility','off');
xlabel('vCar (m/s)'); ylabel(sprintf('d(peak g)/d%s   [%s]',param,unit));
t = sprintf('Peak capability sensitivity to %s',param);
if anyBad
    t = [t,'   (grey x = perturbation levels disagree in sign)'];
end
title(t); legend('Location','best');

%% ---- Figure 5: lap time ----
tf = fieldnames(P);
tf = tf(startsWith(tf,'d_') & endsWith(tf,'_dp'));
if ~isempty(tf)
    figure('Name',sprintf('%s sensitivity: lap time',param), ...
           'Position',[200 200 620 420]);
    setPlotFont(gcf);
    % anchored strip; erase(...,'d_') would also eat the 'd_' in skidpad_dp
    names = regexprep(tf,'^d_(.*)_dp$','$1');
    yv = cellfun(@(f) P.(f),tf);
    bar(categorical(names,names),yv);
    grid on; box on;
    ylabel(sprintf('d(time)/d%s   [s / %s]',param,unit));
    title(sprintf('Lap time sensitivity to %s (negative = faster)',param));
end
end

%% helpers
function plotSlices(x,v,M,R,vSlices,wantC,wantD,xl,yl,ttl)
% R is the reliability mask matching M; unreliable points are dropped from
% the continuous line and drawn as grey x so they stay visible but unread
hold on
co = lines(numel(vSlices));
for k = 1:numel(vSlices)
    i = find(abs(v-vSlices(k))<1e-6,1);
    if isempty(i), continue, end
    y   = M(i,:);
    ok  = isfinite(y) & R(i,:);
    bad = isfinite(y) & ~R(i,:);
    if ~any(ok) && ~any(bad), continue, end
    if wantC
        yPlot = y; yPlot(~ok) = NaN;
        plot(x,yPlot,'-','Color',co(k,:),'LineWidth',1.8, ...
            'DisplayName',sprintf('%.0f m/s',v(i)));
    end
    if wantD
        plot(x(ok),y(ok),'o','Color',co(k,:),'MarkerFaceColor',co(k,:), ...
            'MarkerSize',4,'HandleVisibility',onlyIf(~wantC), ...
            'DisplayName',sprintf('%.0f m/s',v(i)));
    end
    if any(bad)
        plot(x(bad),y(bad),'x','Color',[0.6 0.6 0.6],'MarkerSize',6, ...
            'LineWidth',1.1,'HandleVisibility','off');
    end
end
grid on; box on; yline(0,'k:','HandleVisibility','off');
xlabel(xl); ylabel(yl); title(ttl); legend('Location','best');
end

function surfaceOrMsg(x,y,Z,xl,yl,ttl)
if ~any(isfinite(Z(:)))
    text(0.5,0.5,'no overlapping solved points','HorizontalAlignment','center');
    axis off; title(ttl); return
end
% pcolor tolerates the NaN holes that appear where a case could not reach an
% operating point; contourf would silently close over them
h = pcolor(x,y,Z); set(h,'EdgeColor','none');
shading interp; colorbar; colormap(parula);
xlabel(xl); ylabel(yl); title(ttl);
end

function s = onlyIf(tf)
if tf, s = 'on'; else, s = 'off'; end
end

function u = unitLabel(param)
switch param
    case 'ClA', u = 'g per m^2 ClA';
    case 'CdA', u = 'g per m^2 CdA';
    case 'CoP', u = 'g per unit front balance';
    otherwise,  u = 'per unit';
end
end
