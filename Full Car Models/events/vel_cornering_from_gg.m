function [x_table_corner_vel,radius,max_vel_corner_vector] = vel_cornering_from_gg(car)
% Max cornering velocity against radius, READ OFF the g-g instead of re-solved.
%
% WHY THIS REPLACES 556 fmincon SOLVES
%
% max_vel_cornering and max_lat_accel are the same optimisation written two
% ways: both maximise v*yaw_rate with throttle and all four kappa pinned at 0,
% and car.constraint5(P,r) is car.constraint1(P) plus one equality, P(3)/P(5)=r.
% At the matched velocity constraint1's feasible set therefore CONTAINS
% constraint5's, so wherever the two disagree one of them is off its optimum.
% Measured on the calibration car over all 556 radii they disagree by up to
% 0.844 m/s^2 with both solves feasible to better than 1.4e-7 -- never a
% tolerance difference and never a constraint difference, just two local optima
% of one problem. Which one you land on depends on the seed: max_lat_accel
% takes a low-steer branch from 14 m/s up while the warm-started cornering
% chain rides a high-steer one, and the g-g's own lateral envelope is
% non-monotone in velocity at 3 to 11 of its 29 rows on every car on file.
%
% Underneath that, the tyre has no lateral peak -- Fy rises monotonically in
% |alpha| out past 40 deg -- so "maximum lateral acceleration" here is set by
% the steer and lat_vel box bounds rather than by the tyres, and no amount of
% optimiser work will make two solvers agree on which corner of that box to
% sit in. Solving it once is the only way they can agree.
%
% Solving v^2/r = latmax(v) against the g-g's own envelope cannot disagree with
% the g-g, because it IS the g-g. That is the point: Track_Solver clips lap
% velocity to this table at every station while the g-g supplies the
% longitudinal capability, so any gap between the two parks the interpolant
% query outside the hull where create_scattered_interpolants2 clamps it
% silently.
%
% WHAT IT COSTS, AND IT IS NOT A LEVEL SHIFT. Across 9 grip_sweep cars the
% change runs -0.118 to +0.548 s on autocross and -1.21 to +8.36 s on
% endurance -- FASTER on 4 of the 9, because on those the solved table was the
% pessimistic one and max_vel_cornering had fallen into the worse basin. The
% largest moves are on cars whose g-g lost velocity rows, where this table
% bridges the gap linearly. So re-run the grip calibration after this lands
% rather than assuming an offset; it re-shapes the grip-to-laptime surface.
% Cost is 0.2 s against 8-15 s, which takes the Events2 constructor from about
% 22 s to about 8.
%
% The envelope comes from car.ss_info, which is what max_lat_accel actually
% returned, NOT from car.longAccelLookup. The lookups sit 0.1 m/s^2 below it at
% 25 of 29 velocities because gg2 samples linspace(0.1,maxLat-0.1,20), and up
% to 0.98 m/s^2 below at the rest where makeGG dropped a node. Note the
% consequence: capping the lap at exactly v^2/r = latmax(v) leaves every
% grip-limited apex sitting 0.1 m/s^2 ABOVE the lookup's lateral extent, so the
% envelope clamp in create_scattered_interpolants2 still fires there. Closing
% that needs gg2's 0.1 inset removed, which is a separate decision.

if isempty(car.ss_info)
    error('vel_cornering_from_gg:noGG', ...
        ['car has no g-g results -- run gg2 and makeGG first, or call ' ...
         'vel_cornering_sweep to solve the table on its own.']);
end

radius = 4.5:0.1:60;

% ss_info carries one identical row per lateral grid column, so 20 copies of
% each velocity; collapse before interpolating. Column 6 is long_vel and
% column 3 is x(3)*x(5), the lateral acceleration at the limit.
env = unique(car.ss_info(:,[6 3]),'rows');
[vRows,~,rowIndex] = unique(env(:,1));
latRow = accumarray(rowIndex,env(:,2),[],@max);

% A missing g-g row is a domain gap, not a request to draw a line between the
% nearest solved rows.  gg2/makeGG preserve that information in ggMask.  Keep
% the old clamped behaviour for legacy cars that have no mask, but when a mask
% is present only search within contiguous solved velocity segments.
[expectedVelocity,validExpected,hasMask] = envelopeMask(car,vRows);
if hasMask && any(~validExpected)
    segmentIds = contiguousSegments(validExpected);
else
    expectedVelocity = vRows;
    validExpected = true(size(vRows));
    segmentIds = {1:numel(vRows)};
end
segments = cell(size(segmentIds));
for s = 1:numel(segmentIds)
    idx = segmentIds{s};
    segments{s} = struct('velocity',expectedVelocity(idx), ...
        'lateral',latAtVelocity(expectedVelocity(idx),vRows,latRow), ...
        'startsAtMinimum',idx(1) == 1);
end

v_top = cappedEnvelopeTopVelocity(car,expectedVelocity,validExpected,hasMask);

n = numel(radius);
max_vel_corner_vector = nan(1,n);
for i = 1:n
    r = radius(i);
    for s = 1:numel(segments)
        segment = segments{s};
        v = segment.velocity;
        lat = segment.lateral;
        if isempty(v), continue, end

        % g(v) = v^2/r - latmax(v) is convex on every breakpoint interval,
        % since v^2/r is convex and latmax is linear there.  For the first
        % segment, preserve the historical flat clamp below its first solved
        % velocity.  Later segments start at a real solved row; crossing a
        % missing interval is deliberately never allowed.
        if segment.startsAtMinimum
            lo = 0;
            if gap(v(1),r,v,lat) >= 0
                hi = v(1);
            else
                lo = v(1);
                hi = [];
                for k = 2:numel(v)
                    if gap(v(k),r,v,lat) >= 0, hi = v(k); break, end
                    lo = v(k);
                end
            end
        else
            lo = v(1);
            hi = [];
            if gap(lo,r,v,lat) == 0
                hi = lo;
            elseif gap(lo,r,v,lat) < 0
                for k = 2:numel(v)
                    if gap(v(k),r,v,lat) >= 0, hi = v(k); break, end
                    lo = v(k);
                end
            end
        end
        if isempty(hi), continue, end

        for k = 1:60
            mid = 0.5*(lo+hi);
            if gap(mid,r,v,lat) > 0, hi = mid; else, lo = mid; end
        end
        max_vel_corner_vector(i) = 0.5*(lo+hi);
        break
    end

    % A speed-limited result is valid only when the g-g has a solved top row.
    % Without a mask this retains the legacy max_vel clamp.  With a mask, a
    % missing top row means the correct result is unknown, not max_vel.
    if isnan(max_vel_corner_vector(i)) && (~hasMask || all(validExpected))
        max_vel_corner_vector(i) = v_top;
    end
end

% Nothing reads the state-vector table -- Events2 stores it and event_plotter
% only touches radius_vector and max_vel_corner_vector -- and there is no
% per-radius solve to report a state from any more, so it carries the three
% columns that are actually defined here.
x_table_corner_vel = array2table([radius(:) max_vel_corner_vector(:) ...
    max_vel_corner_vector(:).^2./radius(:)], ...
    'VariableNames',{'radius','max_vel_corner','lat_accel'});

end


function v_top = cappedEnvelopeTopVelocity(car,expectedVelocity,validExpected,hasMask)
v_top = car.max_vel;
if ~hasMask
    return
end

solvedVelocity = expectedVelocity(validExpected);
if ~isempty(solvedVelocity)
    v_top = min(v_top,max(solvedVelocity));
end
end


function g = gap(v,r,vRows,latRow)
% >0 means the radius demands more lateral acceleration than the g-g has at
% that speed
g = v^2/r - lininterp1(vRows,latRow,v);
end

function [expectedVelocity,validExpected,hasMask] = envelopeMask(car,vRows)
hasMask = (isobject(car) && isprop(car,'ggMask')) || ...
    (isstruct(car) && isfield(car,'ggMask'));
if ~hasMask || ~isstruct(car.ggMask) || ...
        ~isfield(car.ggMask,'velocity') || ...
        ~isfield(car.ggMask,'lateral') || isempty(car.ggMask.velocity)
    expectedVelocity = vRows;
    validExpected = true(size(vRows));
    hasMask = false;
    return
end
expectedVelocity = unique(round(double(car.ggMask.velocity(:)),6));
rawVelocity = round(double(car.ggMask.velocity(:)),6);
rawMask = logical(car.ggMask.lateral(:));
validExpected = false(size(expectedVelocity));
for k = 1:numel(expectedVelocity)
    validExpected(k) = any(rawMask(abs(rawVelocity-expectedVelocity(k)) <= 1e-9));
end
% A mask can contain a row identity that is not represented in ss_info after
% older saved g-g results are loaded.  It is still invalid for this lookup.
for k = 1:numel(expectedVelocity)
    validExpected(k) = validExpected(k) && any(abs(vRows-expectedVelocity(k)) <= 1e-9);
end
end

function segments = contiguousSegments(validRows)
indices = find(validRows);
segments = {};
if isempty(indices), return, end
start = 1;
for k = 2:numel(indices)
    if indices(k) ~= indices(k-1)+1
        segments{end+1} = indices(start:k-1); %#ok<AGROW>
        start = k;
    end
end
segments{end+1} = indices(start:end);
end

function values = latAtVelocity(queryVelocity,vRows,latRow)
values = nan(size(queryVelocity));
for k = 1:numel(queryVelocity)
    hit = find(abs(vRows-queryVelocity(k)) <= 1e-9,1);
    if ~isempty(hit), values(k) = latRow(hit); end
end
end
