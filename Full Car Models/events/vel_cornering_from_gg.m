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
env    = unique(car.ss_info(:,[6 3]),'rows');
vRows  = env(:,1);
latRow = env(:,2);

% Breakpoints of the piecewise-linear envelope, plus the ends of the search.
% lininterp1 holds its end values flat, which is the treatment clamped_query
% already gives a velocity outside the solved range.
v_top = car.max_vel;
brk = unique([0; vRows; v_top]);
brk = brk(brk <= v_top);

n = numel(radius);
max_vel_corner_vector = zeros(1,n);
for i = 1:n
    r = radius(i);

    % g(v) = v^2/r - latmax(v) is convex on every breakpoint interval, since
    % v^2/r is convex and latmax is linear there. A convex g cannot go positive
    % between two negative endpoints, so scanning the breakpoints for the FIRST
    % sign change and bisecting inside it finds the lowest root exactly. That
    % matters because the envelope is not monotone in velocity, so a plain
    % bisection over the whole range could settle on a higher crossing and
    % report a cornering speed the car cannot hold.
    lo = brk(1);
    hi = [];
    for k = 2:numel(brk)
        if gap(brk(k),r,vRows,latRow) >= 0, hi = brk(k); break, end
        lo = brk(k);
    end
    if isempty(hi)
        max_vel_corner_vector(i) = v_top;   % speed-limited, not grip-limited
        continue
    end

    for k = 1:60
        mid = 0.5*(lo+hi);
        if gap(mid,r,vRows,latRow) > 0, hi = mid; else, lo = mid; end
    end
    max_vel_corner_vector(i) = 0.5*(lo+hi);
end

% Nothing reads the state-vector table -- Events2 stores it and event_plotter
% only touches radius_vector and max_vel_corner_vector -- and there is no
% per-radius solve to report a state from any more, so it carries the three
% columns that are actually defined here.
x_table_corner_vel = array2table([radius(:) max_vel_corner_vector(:) ...
    max_vel_corner_vector(:).^2./radius(:)], ...
    'VariableNames',{'radius','max_vel_corner','lat_accel'});

end


function g = gap(v,r,vRows,latRow)
% >0 means the radius demands more lateral acceleration than the g-g has at
% that speed
g = v^2/r - lininterp1(vRows,latRow,v);
end
