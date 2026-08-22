function plotRampSweep(R,which)
% Plots a rampSweep result.
%
% inputs  R:     output of rampSweep
%         which: 'all' (default) | 'balance' | 'understeer' | 'ramps'
%
% Figures
%   balance     the metrics that matter over vCar: understeer gradient,
%               mechanical / aero / grip balance, axle load split, downforce
%   understeer  steer excess vs gLat for every speed, which is the raw
%               material K is fitted from -- look here before trusting K
%   ramps       per-corner and per-axle behaviour along each ramp
%
% Where a ramp was truncated (ramp_complete < 1, the car cannot hold that
% speed at the limit) the affected markers are drawn hollow. Those points are
% not at the car's limit and reading them as limit behaviour is wrong.

if nargin < 2 || isempty(which), which = 'all'; end
S = R.perSpeed;
P = R.points;
v = S.vCar;

trunc = S.ramp_complete < 0.98;
mode  = R.settings.mode;

wantB = any(strcmpi(which,{'all','balance'}));
wantU = any(strcmpi(which,{'all','understeer'}));
wantR = any(strcmpi(which,{'all','ramps'}));

%% ---- Figure 1: balance over vCar ----
if wantB
figure('Name',sprintf('ramp sweep: balance over vCar (%s)',mode), ...
       'Position',[60 60 1280 800]);
setPlotFont(gcf);   % must precede the subplots: these are DEFAULT properties

subplot(2,3,1); hold on
plot(v,S.K_linear,'-o','LineWidth',1.8,'DisplayName','linear range');
plot(v,S.K_at_limit,'-s','LineWidth',1.4,'DisplayName','top of ramp');
markTrunc(v,S.K_at_limit,trunc);
yline(0,'k:','HandleVisibility','off');
grid on; box on; xlabel('vCar (m/s)'); ylabel('K (deg/g)');
title('understeer gradient'); legend('Location','best','FontSize',7);
subtitle('above 0 = understeer','FontSize',7);

subplot(2,3,2); hold on
plot(v,S.mech_balance,'-o','LineWidth',1.8,'DisplayName','LLTD (mechanical)');
plot(v,S.aero_balance,'-s','LineWidth',1.8,'DisplayName','CoP (aero)');
plot(v,S.front_Fz_frac_limit,'-^','LineWidth',1.4,'DisplayName','front Fz share');
markTrunc(v,S.front_Fz_frac_limit,trunc);
grid on; box on; xlabel('vCar (m/s)'); ylabel('front share (-)');
title('balance terms'); legend('Location','best','FontSize',7);
subtitle('LLTD and CoP are constants by construction','FontSize',7);

subplot(2,3,3); hold on
plot(v,S.grip_balance_mid,'-o','LineWidth',1.8,'DisplayName','mid ramp');
plot(v,S.grip_balance_limit,'-s','LineWidth',1.4,'DisplayName','top of ramp');
markTrunc(v,S.grip_balance_limit,trunc);
yline(0,'k:','HandleVisibility','off');
grid on; box on; xlabel('vCar (m/s)'); ylabel('\mu_f used - \mu_r used');
title('grip balance'); legend('Location','best','FontSize',7);
subtitle('above 0 = front saturating first','FontSize',7);

subplot(2,3,4); hold on
plot(v,S.LLT_norm_balance_mid,'-o','LineWidth',1.8,'DisplayName','mid ramp');
plot(v,S.LLT_norm_balance_limit,'-s','LineWidth',1.4,'DisplayName','top of ramp');
markTrunc(v,S.LLT_norm_balance_limit,trunc);
yline(0,'k:','HandleVisibility','off');
grid on; box on; xlabel('vCar (m/s)'); ylabel('LLT/Fz front - rear');
title('load transfer, normalised by axle load');
legend('Location','best','FontSize',7);
subtitle('this is the term aero actually moves','FontSize',7);

subplot(2,3,5); hold on
yyaxis left
plot(v,S.gLat_max,'-o','LineWidth',1.8,'DisplayName','free limit');
plot(v,S.gLat_top,'-s','LineWidth',1.4,'DisplayName','ramp top');
ylabel('gLat (g)');
yyaxis right
plot(v,S.downforce_frac,'-^','LineWidth',1.4,'DisplayName','downforce share of load');
plot(v,S.drag_decel,'-v','LineWidth',1.4,'DisplayName','drag (g)');
ylabel('fraction / g');
grid on; box on; xlabel('vCar (m/s)');
title('capability and aero'); legend('Location','best','FontSize',7);

subplot(2,3,6); hold on
plot(v,S.alpha_balance_mid,'-o','LineWidth',1.8,'DisplayName','mid ramp');
plot(v,S.alpha_balance_limit,'-s','LineWidth',1.4,'DisplayName','top of ramp');
markTrunc(v,S.alpha_balance_limit,trunc);
yline(0,'k:','HandleVisibility','off');
grid on; box on; xlabel('vCar (m/s)'); ylabel('|\alpha_f| - |\alpha_r| (deg)');
title('axle slip balance'); legend('Location','best','FontSize',7);
subtitle('above 0 = understeer','FontSize',7);

sgtitle(sprintf(['ramp sweep, %s mode  |  %d speeds, %d pts/ramp  |  ' ...
    'hollow markers = ramp truncated, not at the limit'], ...
    mode, height(S), R.settings.nRamp));
end

%% ---- Figure 2: the raw understeer curves ----
if wantU
figure('Name','ramp sweep: steer excess vs gLat','Position',[90 90 1150 480]);
setPlotFont(gcf);
cmap = parula(max(height(S),2));

subplot(1,2,1); hold on
for i = 1:height(S)
    m = P.vCar == v(i);
    plot(P.gLat(m),P.steer_excess(m),'-o','Color',cmap(i,:),'LineWidth',1.5, ...
        'MarkerSize',3.5,'DisplayName',sprintf('%.0f m/s',v(i)));
end
grid on; box on; xlabel('gLat (g)'); ylabel('\delta - \delta_{Ackermann} (deg)');
title('steer excess: K is the slope of this');
legend('Location','northwest','FontSize',7);

subplot(1,2,2); hold on
for i = 1:height(S)
    m = P.vCar == v(i);
    plot(P.gLat(m),P.alpha_balance(m),'-o','Color',cmap(i,:),'LineWidth',1.5, ...
        'MarkerSize',3.5,'DisplayName',sprintf('%.0f m/s',v(i)));
end
yline(0,'k:','HandleVisibility','off');
grid on; box on; xlabel('gLat (g)'); ylabel('|\alpha_f| - |\alpha_r| (deg)');
title('axle slip balance along the ramp');

sgtitle(sprintf(['raw ramp curves (%s mode) -- check these are smooth before ' ...
    'trusting any fitted gradient'],mode));
end

%% ---- Figure 3: per-ramp detail ----
if wantR
figure('Name','ramp sweep: axle detail','Position',[120 120 1150 800]);
setPlotFont(gcf);
cmap = parula(max(height(S),2));

subplot(2,2,1); hold on
for i = 1:height(S)
    m = P.vCar == v(i);
    plot(P.gLat(m),P.mu_f_used(m),'-','Color',cmap(i,:),'LineWidth',1.6, ...
        'DisplayName',sprintf('%.0f front',v(i)));
    plot(P.gLat(m),P.mu_r_used(m),'--','Color',cmap(i,:),'LineWidth',1.2, ...
        'DisplayName',sprintf('%.0f rear',v(i)));
end
grid on; box on; xlabel('gLat (g)'); ylabel('Fy_{axle} / Fz_{axle}');
title('grip used per axle (solid front, dashed rear)');

subplot(2,2,2); hold on
for i = 1:height(S)
    m = P.vCar == v(i);
    plot(P.gLat(m),P.front_Fz_frac(m),'-','Color',cmap(i,:),'LineWidth',1.6, ...
        'DisplayName',sprintf('%.0f m/s',v(i)));
end
grid on; box on; xlabel('gLat (g)'); ylabel('front Fz / total Fz');
title('axle load split along the ramp');
legend('Location','best','FontSize',7);

subplot(2,2,3); hold on
for i = 1:height(S)
    m = P.vCar == v(i);
    plot(P.gLat(m),P.beta(m),'-','Color',cmap(i,:),'LineWidth',1.6);
end
grid on; box on; xlabel('gLat (g)'); ylabel('\beta (deg)');
title('body slip angle');

subplot(2,2,4); hold on
for i = 1:height(S)
    m = P.vCar == v(i);
    plot(P.gLat(m),P.min_Fz(m),'-','Color',cmap(i,:),'LineWidth',1.6);
end
yline(0,'r:','LineWidth',1.2,'HandleVisibility','off');
grid on; box on; xlabel('gLat (g)'); ylabel('min corner Fz (N)');
title('lightest wheel (below the red line = wheel lift)');

sgtitle(sprintf('per-ramp detail (%s mode), colour = speed',mode));
end
end

function markTrunc(x,y,trunc)
% hollow overlay on speeds whose ramp never reached the limit
if ~any(trunc), return, end
plot(x(trunc),y(trunc),'o','MarkerSize',11,'LineWidth',1.4, ...
    'MarkerFaceColor','none','Color',[0.85 0.33 0.1], ...
    'DisplayName','ramp truncated');
end
