function A = lapTimeAttribution(trace,vGrid,thGrid,r0,drdp)
% First-order attribution of a lap time change to WHERE ON THE LAP it happens.
%
% A g-g envelope that grows by dr at operating condition (v,theta) makes the
% car quicker only in the parts of the lap that (a) actually pass through that
% condition and (b) are actually limited by the envelope while they do.
%
% Model, evaluated per trace sample s:
%     phi_s  = dr(v_s,th_s) / r0(v_s,th_s)      fractional envelope gain
%     util_s = min(r_s / r0(v_s,th_s), 1)       how close to the limit it is
%     dt_s   = -0.5 * w_s * phi_s * util_s      time delta carried by sample s
%
% The -1/2 exponent is not a fudge -- t ~ a^(-1/2) holds for all three segment
% types in a QSS lap, so the same factor applies wherever on the envelope the
% limiting happens:
%     corner at fixed radius   v = sqrt(a_y R),  t = s/v      -> t ~ a^(-1/2)
%     acceleration over s      s = a t^2 / 2                  -> t ~ a^(-1/2)
%     braking over s           same                           -> t ~ a^(-1/2)
%
% util_s matters because a lap is not on the limit everywhere. Time spent
% coasting mid-straight does not shrink when the envelope grows, and weighting
% it as if it did would inflate the total and misplace it on the vCar axis.
%
% Samples falling outside the solved envelope (NaN in r0 or drdp -- unreached
% operating points, or cells masked as unreliable) carry no attribution at
% all. A.coverage reports how much of the lap that silently cost, and callers
% should surface it: a decomposition covering 40% of the lap is not a
% decomposition.
%
% inputs  trace:  struct with t, v, gLat, gLong (from aeroLapTimes)
%         vGrid:  velocity grid of the envelope   [nV x 1]
%         thGrid: theta grid of the envelope      [1 x nTheta]
%         r0:     baseline envelope radius        [nV x nTheta]
%         drdp:   d(radius)/d(param)              [nV x nTheta]
% output  A:      struct with per-sample dt and the four operating axes
%                 (vCar, gLat, gLong, theta), plus .total, .coverage,
%                 .lapTime and .n
%
% Sign convention: negative dt is faster. A parameter that helps produces a
% negative total, matching the sign of the lap sim's own measured delta.

A = struct('dt',[],'vCar',[],'gLat',[],'gLong',[],'theta',[], ...
           'total',0,'coverage',0,'lapTime',NaN,'n',0);

if isempty(trace) || ~isstruct(trace) || ~isfield(trace,'t') ...
        || numel(trace.t) < 3 || ~isfield(trace,'gLat') || ~isfield(trace,'gLong') ...
        || isempty(trace.gLat) || isempty(trace.gLong)
    return
end

t  = trace.t(:);
v  = trace.v(:);
gy = abs(trace.gLat(:));
gx = trace.gLong(:);

% the event solver can return profiles one sample longer or shorter than the
% time vector; trim rather than guess which end is stale
m  = min([numel(t) numel(v) numel(gy) numel(gx)]);
if m < 3, return, end
t = t(1:m); v = v(1:m); gy = gy(1:m); gx = gx(1:m);

% time carried by each sample: midpoint differences, so the first and last
% samples keep their share instead of being dropped
w = zeros(m,1);
dt = diff(t);
w(1)       = dt(1)/2;
w(end)     = dt(end)/2;
w(2:end-1) = (dt(1:end-1)+dt(2:end))/2;

th = atan2d(gy,gx);       % 0 accel, 90 pure lateral, 180 braking
r  = hypot(gx,gy);

% lapTime is the WHOLE lap, measured from the time vector alone. Deriving it
% from the usable samples instead would let a trace full of NaN accelerations
% report 100% coverage of a lap it barely touched.
hasT = isfinite(w) & w > 0;
A.lapTime = sum(w(hasT));

ok = hasT & isfinite(v) & isfinite(th) & isfinite(r);
if ~any(ok), return, end

% interp2 takes X across columns (theta) and Y down rows (velocity), which is
% the layout aeroEnvelope produces. No extrapolation: outside the solved
% envelope the answer is NaN, not an invented gradient.
r0q = interp2(thGrid(:).',vGrid(:),r0,   th,v,'linear',NaN);
drq = interp2(thGrid(:).',vGrid(:),drdp, th,v,'linear',NaN);

use = ok & isfinite(r0q) & isfinite(drq) & r0q > 0;
if ~any(use), return, end

phi  = drq(use)./r0q(use);
util = min(r(use)./r0q(use),1);

A.dt       = -0.5*w(use).*phi.*util;
A.vCar     = v(use);
A.gLat     = gy(use);
A.gLong    = gx(use);
A.theta    = th(use);
A.total    = sum(A.dt);
A.coverage = sum(w(use))/max(A.lapTime,eps);
A.n        = sum(use);
end
