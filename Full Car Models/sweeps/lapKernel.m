function [dens,grid,info] = lapKernel(trace,which,grid,bw)
% Time-weighted density ("kernel") of a lap trace over vCar, gLat, gLong or
% theta -- i.e. how much time the car actually spends at each condition.
%
% This is the weighting that turns a condition-dependent sensitivity into a
% lap-relevant number: sensitivity where the car never operates is worth
% nothing.
%
% inputs  trace: struct with t, v, gLat, gLong (from aeroLapTimes)
%         which: 'vCar' | 'gLat' | 'gLong' | 'theta'
%         grid:  evaluation points (optional)
%         bw:    bandwidth (optional; Silverman's rule by default)
% outputs dens:  density on grid, integrates to 1 over the grid
%         grid:  the grid used
%         info:  struct with .bw, .totalTime, .n, .mean (time-weighted)
%
% Implemented as a weighted Gaussian KDE in base MATLAB -- ksdensity lives in
% the Statistics toolbox, which this install does not have.

dens = []; info = struct('bw',NaN,'totalTime',0,'n',0,'mean',NaN);
if nargin < 3, grid = []; end
if nargin < 4, bw = []; end

if isempty(trace) || isempty(trace.t) || numel(trace.t) < 3
    if isempty(grid), grid = linspace(0,1,10); end
    dens = nan(size(grid));
    return
end

% lims are the physically meaningful range of each variable; the default grid
% is clamped to them so the kernel is never drawn at, say, theta = -50 deg
switch lower(which)
    case 'vcar',  x = trace.v;                             lims = [0 Inf];
    case 'glat',  x = abs(trace.gLat);                     lims = [0 Inf];
    case 'glong', x = trace.gLong;                         lims = [-Inf Inf];
    case 'theta', x = atan2d(abs(trace.gLat),trace.gLong); lims = [0 180];
    otherwise, error('lapKernel:badVar','unknown variable %s',which);
end
x = x(:);

% time spent at each sample: midpoint differences of the time vector, so the
% first and last samples are not silently dropped
t = trace.t(:);
w = zeros(size(t));
if numel(t) > 1
    dt = diff(t);
    w(1)       = dt(1)/2;
    w(end)     = dt(end)/2;
    w(2:end-1) = (dt(1:end-1)+dt(2:end))/2;
end

ok = isfinite(x) & isfinite(w) & w > 0;
x = x(ok); w = w(ok);
if numel(x) < 3
    if isempty(grid), grid = linspace(0,1,10); end
    dens = nan(size(grid));
    return
end
w = w/sum(w);

if isempty(grid)
    pad = 0.05*(max(x)-min(x)+eps);
    lo = max(min(x)-pad, lims(1));
    hi = min(max(x)+pad, lims(2));
    grid = linspace(lo,hi,200);
end

% Silverman on the time-weighted sample, using Kish effective sample size
mu  = sum(w.*x);
sd  = sqrt(max(sum(w.*(x-mu).^2),eps));
nEff = 1/sum(w.^2);
if isempty(bw)
    bw = 1.06*sd*nEff^(-1/5);
    bw = max(bw,(max(x)-min(x))/200 + eps);   % never collapse to a spike
end

g = grid(:);
% (nGrid x nSample) kernel matrix; traces are a few thousand samples so this
% is cheap and avoids a loop
Z = (g - x.')/bw;
dens = (exp(-0.5*Z.^2)*(w))/(bw*sqrt(2*pi));
dens = reshape(dens,size(grid));

% normalise over the grid so different cases are comparable
A = trapz(grid,dens);
if A > 0, dens = dens/A; end

info.bw = bw;
info.totalTime = sum(trace.t(end)-trace.t(1));
info.n = numel(x);
info.mean = mu;
end
