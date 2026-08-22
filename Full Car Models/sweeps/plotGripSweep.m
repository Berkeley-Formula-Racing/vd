function plotGripSweep(G)
% Contour maps of a gripSweep over the (front, rear) grip grid.
%
% One panel per event, one for balance if it was evaluated, and one for the
% combined score if anything was targeted. Where a target exists its contour
% is drawn heavy: every (front,rear) pair on that line reproduces the
% measured value, which is the point -- a single event does not identify a
% pair, it identifies a CURVE of pairs. Two events crossing, or one event
% crossing the balance line, is what pins the answer down.
%
% Axes are fitted to the swept grid, so a narrow sweep fills the panel. The
% black circle marking the unmodified tyre at (1,1) is therefore only drawn
% when the grid reaches it; the subtitle notes when it does not.
%
% input  G: output of gripSweep

T = G.times;
F = G.grid.F; R = G.grid.R;
sz = size(F);

% contourf needs a 2x2 grid at minimum, and says so with an error that names
% Z rather than the sweep. A one-value axis is a plausible thing to ask for --
% sweeping only the rear, or checking a single pair -- so name it here.
if numel(unique(F(:))) < 2 || numel(unique(R(:))) < 2
    error('plotGripSweep:degenerateGrid', ...
        ['contour maps need at least 2 front and 2 rear values; this sweep ' ...
         'has %d x %d. Widen opts.front / opts.rear, or read G.times ' ...
         'directly for a single pair.'], numel(unique(F(:))),numel(unique(R(:))));
end
ev = G.settings.events;
if ischar(ev), ev = {ev}; end

panels = ev;
if ismember('balance',T.Properties.VariableNames), panels{end+1} = 'balance'; end
if ismember('score',T.Properties.VariableNames),   panels{end+1} = 'score';   end

nP = numel(panels);
nc = min(3,nP); nr = ceil(nP/nc);

fig = figure('Name','grip sweep','Position',[60 60 420*nc 380*nr]);
setPlotFont(fig);

for k = 1:nP
    subplot(nr,nc,k);
    name = panels{k};
    Z = reshape(T.(name),sz);

    if all(~isfinite(Z(:)))
        text(0.5,0.5,sprintf('%s: nothing solved',name), ...
            'HorizontalAlignment','center'); axis off; continue
    end

    % LineColor takes an RGB triplet only -- a 4-element RGBA errors out
    contourf(F,R,Z,14,'LineColor',[0.85 0.85 0.85]);
    colormap(gca,parula); c = colorbar; c.Label.String = unitFor(name);
    hold on

    % the target contour, if this panel has one
    % No legend: 'best' placement lands on the title in most panels, and the
    % red line is already named in the title and the figure subtitle.
    tv = targetValue(G,name);
    if ~isempty(tv) && tv > min(Z(:)) && tv < max(Z(:))
        contour(F,R,Z,[tv tv],'LineColor','r','LineWidth',2.5);
    end

    % the best pair, on every panel, so it can be read against each map
    if ismember('score',T.Properties.VariableNames) && isfield(G,'best')
        plot(G.best.gripF,G.best.gripR,'p','MarkerSize',15, ...
            'MarkerFaceColor',[1 0.85 0.1],'MarkerEdgeColor','k', ...
            'HandleVisibility','off');
    end

    % Baseline, the unmodified tyre. Only drawn when the sweep actually covers
    % it: it used to be drawn unconditionally, and on a grid that does not
    % reach 1 it dragged both axes out to the baseline and squashed the swept
    % region into a corner -- a 0.55-0.62 grid ended up occupying about a
    % sixth of the panel. The subtitle says so when it is off-grid.
    if baselineInGrid(F,R)
        plot(1,1,'ko','MarkerSize',7,'LineWidth',1.4,'HandleVisibility','off');
    end

    % fit the axes to what was actually swept
    xlim(fitRange(F)); ylim(fitRange(R));

    axis square; grid on; box on;
    xlabel('front grip scale'); ylabel('rear grip scale');
    title(titleFor(name,tv));
end

s = G.settings;
sub = sprintf(['front %.2f-%.2f, rear %.2f-%.2f  |  red = measured value ' ...
    '(every pair on that line reproduces it)'], ...
    min(s.front),max(s.front),min(s.rear),max(s.rear));
if ~isempty(s.balanceSpeed)
    sub = sprintf('%s  |  balance: %s at %.1f m/s (%s mode)', ...
        sub,s.balanceMetric,s.balanceSpeed,s.balanceMode);
end
if isfield(G,'best') && isfinite(G.best.score)
    sub = sprintf('%s\nbest gripF %.3f, gripR %.3f (score %.2f -- star)', ...
        sub,G.best.gripF,G.best.gripR,G.best.score);
end
% Say it rather than let the reader wonder where the reference marker went.
% Axes are fitted to the grid, so on a sweep that does not reach 1 there is no
% black circle to find.
if ~baselineInGrid(F,R)
    sub = sprintf('%s  |  unmodified tyre (1,1) is outside this grid',sub);
end
sgtitle(sprintf('Grip sweep\n%s',sub),'FontSize',10);
end

function r = fitRange(V)
% Axis limits for one swept axis. A single-point sweep has no range to fit, so
% it gets a small symmetric window instead -- xlim rejects a zero-width range.
lo = min(V(:)); hi = max(V(:));
if hi > lo
    r = [lo hi];
else
    pad = max(0.01,abs(lo)*0.05);
    r = [lo-pad lo+pad];
end
end


function tf = baselineInGrid(F,R)
tf = 1 >= min(F(:)) && 1 <= max(F(:)) && 1 >= min(R(:)) && 1 <= max(R(:));
end


function tv = targetValue(G,name)
tv = [];
if strcmp(name,'balance')
    if ~isempty(G.settings.balanceSpeed), tv = G.settings.balanceTarget; end
elseif isfield(G.target,name)
    tv = G.target.(name);
end
end

function s = titleFor(name,tv)
switch name
    case 'score',   s = 'combined match score (lower better)';
    case 'balance', s = 'balance at the target speed';
    otherwise,      s = sprintf('%s time',name);
end
if ~isempty(tv), s = sprintf('%s   [target %.4g]',s,tv); end
end

function u = unitFor(name)
switch name
    case 'score',   u = 'RMS tolerances';
    case 'balance', u = '-';
    otherwise,      u = 's';
end
end
