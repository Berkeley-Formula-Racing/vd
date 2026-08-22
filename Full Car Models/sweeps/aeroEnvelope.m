function E = aeroEnvelope(car,grids)
% Extracts the performance envelope of one solved car onto common grids, so
% that cases with different peak capability can be differenced against each
% other. Called by aeroSensitivity.
%
% inputs  car:   car that has been through gg2 + makeGG
%         grids: struct with fields vCar, theta, gLat, gLong
% output  E:     struct of [nV x nX] matrices, NaN where the car cannot
%                reach that operating point (never extrapolated)
%
% Column layouts of the *_info matrices (see max_lat_accel.m etc.):
%   col 1 exitflag | col 2 long_accel | col 3 lat_accel | col 4:12 state
%   state(3) = longitudinal velocity -> col 6, state(5) = yaw rate -> col 8
%
% Lateral acceleration is derived here as v*r (col6 .* col8) rather than read
% from col 3. Both agree, but deriving it keeps this function independent of
% the column-3 convention, which was previously storing the solver residual.

g = car.g;
v  = grids.vCar(:);
nV = numel(v);

ss = filterSolved(car.ss_info);
ac = filterSolved(car.accel_info);
dc = filterSolved(car.decel_info);

E.maxGLat   = nan(nV,1);
E.maxGLong  = nan(nV,1);   % pure-ish acceleration, at the lowest gLat solved
E.maxGBrake = nan(nV,1);
E.r_theta   = nan(nV,numel(grids.theta));
E.accel_at_gLat = nan(nV,numel(grids.gLat));
E.brake_at_gLat = nan(nV,numel(grids.gLat));
E.gLat_at_gLong = nan(nV,numel(grids.gLong));

for i = 1:nV
    % rows belonging to this velocity (the gg2 grid is exact, so match tightly)
    sSel = velRows(ss,v(i));
    aSel = velRows(ac,v(i));
    dSel = velRows(dc,v(i));

    if any(sSel)
        E.maxGLat(i) = max(latG(ss(sSel,:),g));
    end
    if any(aSel)
        aLat  = latG(ac(aSel,:),g);   aLong = ac(aSel,2)/g;
        E.maxGLong(i) = max(aLong);
        E.accel_at_gLat(i,:) = interpMonotone(aLat,aLong,grids.gLat);
    end
    if any(dSel)
        dLat  = latG(dc(dSel,:),g);   dLong = dc(dSel,2)/g;
        E.maxGBrake(i) = min(dLong);
        E.brake_at_gLat(i,:) = interpMonotone(dLat,dLong,grids.gLat);
    end

    % polar form of the whole envelope: r(theta) at this velocity
    pts = [ac(aSel,2)/g latG(ac(aSel,:),g); ...
           dc(dSel,2)/g latG(dc(dSel,:),g); ...
           ss(sSel,2)/g latG(ss(sSel,:),g)];   % [gLong gLat]
    if size(pts,1) >= 2
        th = atan2d(abs(pts(:,2)),pts(:,1));
        r  = hypot(pts(:,1),pts(:,2));
        E.r_theta(i,:) = interpMonotone(th,r,grids.theta);
        % achievable lateral g at a prescribed longitudinal g
        E.gLat_at_gLong(i,:) = interpMonotone(pts(:,1),abs(pts(:,2)),grids.gLong);
    end
end
end

function a = latG(rows,g)
% lateral acceleration in g, derived as v*r from the state vector
a = rows(:,6).*rows(:,8)/g;
end

function out = filterSolved(info)
% keep only converged solves; exitflag 1 = converged, 2 = step tolerance met
if isempty(info)
    out = zeros(0,12);
else
    out = info(info(:,1) == 1 | info(:,1) == 2, :);
end
end

function sel = velRows(info,v)
if isempty(info)
    sel = false(0,1);
else
    sel = abs(info(:,6) - v) < 1e-6;
end
end

function yq = interpMonotone(x,y,xq)
% average duplicate x, sort, then interpolate with NO extrapolation -- an
% envelope must never be invented outside the range the solver actually
% reached, or sensitivity at the edges becomes fiction
x = x(:); y = y(:);
ok = isfinite(x) & isfinite(y);
x = x(ok); y = y(ok);
if numel(x) < 2
    yq = nan(size(xq));
    return
end
[xu,~,ic] = unique(round(x,9));
yu = accumarray(ic,y,[],@mean);
if numel(xu) < 2
    yq = nan(size(xq));
    return
end
yq = interp1(xu,yu,xq,'linear',NaN);
yq = reshape(yq,size(xq));
end
