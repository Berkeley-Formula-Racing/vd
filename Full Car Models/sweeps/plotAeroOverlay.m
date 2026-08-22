function plotAeroOverlay(S,param,event,nBins,mode)
% Where on the lap a change in an aero parameter actually buys time.
%
% One figure, four panels (vCar, gLat, gLong, theta). The y axis is LAP TIME
% DELTA per unit of the parameter -- seconds per m^2 of ClA, seconds per unit
% of balance for CoP. Negative is faster.
%
%   BARS  - TIME SAVED per bin, in seconds, by actually making the baseline
%           -> fastest change (left axis). Taller is better; a bar below the
%           zero line marks a region that change made SLOWER. Not per unit of
%           the parameter: per-unit inverts whenever the improving direction
%           is negative -- the fastest CdA case is a drag REDUCTION -- which
%           would push every bar below the axis for a change that helped. The
%           per-unit sensitivity is on the third title line instead, since
%           that is the form that compares across parameters.
%   LINE  - the running cumulative share of that same quantity, 0 to 100%
%           (right axis)
%   SHADE - the lap time-at-condition kernel, background context only
%
% Note the plotted sign is the opposite of the internal one:
% lapTimeAttribution works in dt, where negative is faster, because that
% matches the lap sim's own d(t)/dp. Everything drawn or printed by this
% function is negated into time saved, so the whole figure reads one way.
%
% Bars and line are the SAME data: the line is the bars' running total,
% normalised by the total lap time delta, so it ends at 100% by construction.
% Read it as "by this speed / this much lateral g, the change has delivered
% X% of everything it is going to deliver". The line is computed from the
% per-sample contributions rather than by summing the bars, so it does not
% inherit the bar resolution.
%
% In 'fraction' mode the BARS are divided by the total too, so they sum to 1
% -- the line is a percentage either way.
%
% The underlying envelope quantity is still d(g-g radius)/d(param), but it is
% converted to seconds by lapTimeAttribution, which walks the baseline lap
% trace and books the time saved sample by sample. That conversion is first
% order, so the subtitle prints the attributed total next to the lap sim's own
% measured delta: if those two disagree badly, the shape is still informative
% but the magnitude is not.
%
% inputs  S:      output of aeroSensitivity (needs .traces AND .times)
%         param:  'ClA' | 'CdA' | 'CoP'
%         event:  'autocross' (default) | 'endurance'
%         nBins:  bar count (default 12)
%         mode:   'seconds' (default) | 'fraction'

if nargin < 2 || isempty(param)
    f = fieldnames(S.sens); param = f{1};
end
if nargin < 3 || isempty(event), event = 'autocross'; end
if nargin < 4 || isempty(nBins), nBins = 12; end
if nargin < 5 || isempty(mode),  mode  = 'seconds'; end
if ~any(strcmpi(mode,{'seconds','fraction'}))
    error('plotAeroOverlay:badMode','mode must be ''seconds'' or ''fraction''');
end
frac = strcmpi(mode,'fraction');
if ~isfield(S.sens,param)
    error('plotAeroOverlay:noParam','no sensitivity computed for %s',param);
end

P = S.sens.(param);
if ~isfield(P,'compareCase')
    error('plotAeroOverlay:noTimes', ...
        ['need lap times to identify the fastest setting -- run ' ...
         'AeroSensitivityStudy with runLapTimes = true']);
end

base = S.baseIdx;
cmp  = P.compareCase;
dpc  = P.compareDp;
v    = S.grids.vCar(:);
th   = S.grids.theta(:).';

%% the baseline lap trace -- everything below is attributed onto it
tr = [];
if isfield(S,'traces') && ~isempty(S.traces) && numel(S.traces) >= base ...
        && ~isempty(S.traces{base}) && isfield(S.traces{base},event)
    tr = S.traces{base}.(event);
end
if isempty(tr) || ~isfield(tr,'t') || numel(tr.t) < 3
    error('plotAeroOverlay:noTrace', ...
        ['no %s trace for the baseline case, so there is no lap to attribute ' ...
         'time onto. Run AeroSensitivityStudy with runLapTimes = true and ' ...
         'pass opts.traces into aeroSensitivity.'],event);
end

%% envelope slopes -> seconds, via the trace
rBase = S.env{base}.r_theta;             % [nV x nTheta]
dDisc = (S.env{cmp}.r_theta - rBase)/dpc;  % two-case discrete slope

% Only the discrete comparison is attributed. Attributing the multi-level
% field as well meant a second pass over ~100k trace samples whose result is
% not drawn -- the scalar multi-level lap-time slope below is enough to tell
% whether the two-case comparison is representative.
aD = lapTimeAttribution(tr,v,th,rBase,dDisc);

if aD.n == 0
    error('plotAeroOverlay:noOverlap', ...
        ['the %s lap never enters the solved envelope, so no time can be ' ...
         'attributed. Check that the g-g and the event ran on the same car.'],event);
end

% what the lap sim itself measured, for the honesty check in the subtitle
mDisc = NaN; mCont = NaN; Tev = NaN;
if isfield(S,'times') && ismember(event,S.times.Properties.VariableNames)
    mDisc = (S.times.(event)(cmp) - S.times.(event)(base))/dpc;
    Tev   = S.times.(event)(base);
end
if isfield(P,['d_' event '_dp'])
    mCont = P.(['d_' event '_dp']);
end

% The trace and the event result can cover different distances: the endurance
% trace is a SINGLE lap while the endurance time is the whole ~11-lap event,
% so an unscaled attribution comes out ~11x short and looks like a broken
% model rather than a units mismatch. Scale onto the event using the ratio of
% durations -- the raw ratio, not a rounded lap count, since the trace need
% not start and stop exactly on the timing line.
nRep = 1;
if isfinite(Tev) && aD.lapTime > 0
    r = Tev/aD.lapTime;
    if r > 0.5 && r < 50, nRep = r; end
end
aD = scaleAttr(aD,nRep);

% Per-unit is the right basis for comparing PARAMETERS against each other,
% and it is what the subtitle and S.sens report. It is the wrong basis for
% the bars: it flips sign whenever the improving direction is negative. The
% fastest CdA case is a drag REDUCTION, so per-unit puts every bar below the
% axis for a change that made the car quicker, and "negative bar = slower
% section" stops meaning anything. Scaling by dp puts the bars in seconds for
% the actual baseline -> fastest change, which reads the same way for all
% three parameters: up is the change helping, down is it hurting.
aDp = scaleAttr(aD,dpc);

%% panels
axes4 = { 'vCar',  'vCar (m/s)'; ...
          'gLat',  '|gLat| (g)'; ...
          'gLong', 'gLong (g)'; ...
          'theta', '\theta (deg)  [0 accel, 90 lateral, 180 brake]'};

if frac
    yUnit = 'fraction of total time saved, per bin';
else
    yUnit = 'time saved per bin  [s]';   % seconds for the actual change
end

fig = figure('Name',sprintf('%s lap time attribution (%s)',param,event), ...
       'Position',[70 70 1250 780]);
setPlotFont(fig);

for k = 1:4
    subplot(2,2,k);
    nm = axes4{k,1};

    xD = getAx(aDp,nm);
    % Plotted sign is flipped from the internal one: lapTimeAttribution
    % returns dt where negative means faster (matching the lap sim's own
    % d(t)/dp), but a bar chart reads far better when taller = better. So the
    % bars show TIME SAVED, and a bar below the zero line is a region the
    % change made slower.
    yD = -aDp.dt;

    if isempty(xD) || ~any(isfinite(xD))
        text(0.5,0.5,'nothing attributed on this axis', ...
            'HorizontalAlignment','center'); axis off; title(nm); continue
    end
    lo = min(xD); hi = max(xD);
    if hi <= lo, hi = lo + eps(lo) + 1e-9; end

    edges = linspace(lo,hi,nBins+1);
    ctr   = (edges(1:end-1)+edges(2:end))/2;
    fine  = linspace(lo,hi,300);

    barY = binSum(xD,yD,edges);
    sD   = sum(yD);
    if frac && abs(sD) > 0, barY = barY/sD; end

    yyaxis left
    bar(ctr,barY,1.0,'FaceColor',[0.35 0.55 0.85],'FaceAlpha',0.55, ...
        'EdgeColor',[0.2 0.3 0.5],'DisplayName','time saved per bin (baseline vs fastest)');
    hold on
    yline(0,'k:','HandleVisibility','off');
    ylabel(yUnit);
    ax = gca; ax.YColor = 'k';

    yyaxis right
    % Lap time density of this panel's own variable, straight from the trace.
    % Background context only -- it says where the car spends time, which is
    % not where the time is won. Scaled to fill the panel, not to be read off
    % the axis.
    kd = lapKernel(tr,nm,fine);
    if ~isempty(kd) && any(isfinite(kd)) && max(kd) > 0
        kd(~isfinite(kd)) = 0;
        kd = 100*kd/max(kd);
        fill([fine fliplr(fine)],[kd zeros(size(kd))],[0.4 0.4 0.4], ...
            'FaceAlpha',0.12,'EdgeColor','none','DisplayName','lap time density');
    end
    % The cumulative curve: running share of the SAME per-bin quantity the
    % bars show, so it is literally their running total and lands on 100%.
    % Computed from the per-sample contributions rather than from the bars, so
    % it does not inherit the bin resolution.
    cum = 100*cumFrac(xD,yD,fine);
    plot(fine,cum,'-','Color',[0.85 0.25 0.1],'LineWidth',2.2, ...
        'DisplayName','cumulative share of total');
    yline(100,'-','Color',[0.85 0.25 0.1 0.35],'LineWidth',0.8, ...
        'HandleVisibility','off');
    ylim([0 108]); ylabel('cumulative share of \Deltat (%)');
    ax = gca; ax.YColor = [0.85 0.25 0.1];

    grid on; box on; xlabel(axes4{k,2});
    if frac
        title(sprintf('vs %s   (\\Sigmabars %.3f, curve ends %.1f%%)', ...
            nm, sum(barY,'omitnan'), cum(end)));
    else
        title(sprintf('vs %s   (\\Sigmabars %+.4g, curve ends %.1f%%)', ...
            nm, sum(barY,'omitnan'), cum(end)));
    end
    if k == 1, legend('Location','best','FontSize',7); end
end

note = '';
if P.baselineIsFastest
    note = '  [baseline already fastest; bars use the best perturbed case]';
end
rep = '';
if abs(nRep-1) > 0.01
    rep = sprintf('  |  trace \\times%.2f to cover the event',nRep);
end

% If the envelope says this change makes the car slower while the lap sim
% says it made it faster (or vice versa), the two disagree about something
% more basic than magnitude. In practice it means the parameter's effect is
% under the lap sim's own noise, so whichever case came out quickest was
% selected from scatter rather than from a real optimum -- and then the
% comparison this whole figure is built on is arbitrary. Say so on the plot.
signWarn = '';
if isfinite(mDisc) && aDp.total*(mDisc*dpc) < 0
    signWarn = ['\newline[!] attributed and measured disagree in SIGN: the effect is ' ...
                'under the lap sim noise floor and the compare case is scatter-selected'];
end

% The multi-level slope uses every perturbation level, so one noisy case
% moves it far less than it moves the two-point comparison the bars are built
% from. Shown only when the two materially disagree, which is the tell that
% the compare case is not representative of the parameter's real effect.
if isempty(signWarn) && isfinite(mCont) && isfinite(mDisc) && abs(mDisc) > 0
    r2 = mCont/mDisc;
    if r2 < 0.5 || r2 > 2
        signWarn = sprintf(['\\newline[!] multi-level slope %+.4g disagrees with the ' ...
            'two-case %+.4g [%s]: the compare case is not representative'], ...
            -mCont, -mDisc, unitLabel(param));
    end
end
% three short lines rather than two long ones: the two-line form ran off both
% edges of the figure and clipped the numbers it exists to show
% Line 2 is the bars' own quantity -- seconds for this specific change.
% Line 3 is the per-unit sensitivity, which is what compares across
% parameters and what S.sens holds. Both are signed as time SAVED.
sgtitle(sprintf(['%s on %s, case %d \\rightarrow %d (\\Delta%s = %+.4g)%s\n' ...
    'this change saves: attributed %+.4g s vs lap sim %+.4g s\n' ...
    'per unit %s: %+.4g vs %+.4g [%s]  |  covers %.0f%% of the lap%s'], ...
    param, event, base, cmp, param, dpc, note, ...
    -aDp.total, -mDisc*dpc, ...
    param, -aD.total, -mDisc, unitLabel(param), ...
    100*aD.coverage, [rep signWarn]), 'FontSize',10);
end

%% helpers
function A = scaleAttr(A,k)
% stretch a single-lap attribution onto a multi-lap event
A.dt    = A.dt*k;
A.total = A.total*k;
end

function x = getAx(A,nm)
if isempty(A.dt), x = []; else, x = A.(nm); end
end

function y = binSum(x,val,edges)
% total attributed time falling in each bin -- a sum, not a mean, so the bars
% add up to the total
n = numel(edges)-1; y = zeros(1,n);
if isempty(x), y(:) = NaN; return, end
for i = 1:n
    if i < n
        m = x >= edges(i) & x < edges(i+1);
    else
        m = x >= edges(i) & x <= edges(i+1);
    end
    if any(m), y(i) = sum(val(m),'omitnan'); end
end
end

function yq = cumFrac(x,val,xq)
% Running share of the total, evaluated at xq. Exact -- it is the empirical
% cumulative sum of the per-sample contributions, not a smoothed or
% integrated density, so it starts at 0 and lands on exactly 1 by
% construction rather than to within a quadrature error.
%
% It is NOT forced monotone. If some region of the axis costs time while
% another gains it -- extra drag hurting on the straight while downforce pays
% in the corners -- the curve rises above 100% and comes back, and that is
% worth seeing rather than hiding.
yq = nan(size(xq));
if isempty(x), return, end
x = x(:); val = val(:);
ok = isfinite(x) & isfinite(val);
x = x(ok); val = val(ok);
if numel(x) < 2, return, end

[xs,idx] = sort(x);
c = cumsum(val(idx));
tot = c(end);
if abs(tot) <= 0, return, end

% collapse duplicate x to their last cumulative value; interp1 'previous'
% rejects repeated sample points
[xu,ia] = unique(xs,'last');
if numel(xu) < 2, return, end

yq = interp1(xu,c(ia)/tot,xq,'previous',NaN);
yq = reshape(yq,size(xq));
yq(xq <  xu(1))   = 0;
yq(xq >= xu(end)) = 1;
end

function u = unitLabel(param)
switch param
    case 'ClA', u = 's per m^2';
    case 'CdA', u = 's per m^2';
    case 'CoP', u = 's per unit balance';
    otherwise,  u = 's per unit';
end
end
