function [F_accel,F_braking] = create_scattered_interpolants2(vel_matrix_accel,vel_matrix_braking)
% input: vel_matrix_accel and vel_matrix_braking obtained from g-g diagram,
%   contains velocity, lateral acceleration, and longitudinal acceleration
% output: F_accel and F_braking, each a handle called as F(lat_accel,long_vel)
%   returning the max possible accel/braking there. These used to be bare
%   scatteredInterpolants; they are handles now because the query has to be
%   projected onto the solved envelope before it is interpolated -- see
%   envelope_lookup below for why. Callers are unaffected: the call syntax is
%   the same, and nothing reads the interpolant's own properties.

long_g_accel = vel_matrix_accel(:,1);
lat_g_accel = vel_matrix_accel(:,2);
vel_accel = vel_matrix_accel(:,3);

long_g_braking = vel_matrix_braking(:,1);
lat_g_braking = vel_matrix_braking(:,2);
vel_braking = vel_matrix_braking(:,3);

F_accel   = envelope_lookup(lat_g_accel,vel_accel,long_g_accel);
F_braking = envelope_lookup(lat_g_braking,vel_braking,long_g_braking);

end


function F = envelope_lookup(lat,vel,long)
% F(lat_accel,long_vel) -> longitudinal capability, with the query projected
% onto the region the g-g actually solved before it is interpolated.
%
% WHY. scatteredInterpolant extrapolates linearly outside the convex hull of
% the cloud and gives no indication that it has, and the lap solver asks for
% states a long way outside it: on the 2024 autocross 20% of the acceleration
% queries and 47% of the braking queries fell outside, because the backward
% pass ran the velocity up unbounded and lat = v^2*kappa followed it. What came
% back is a plane continued off the boundary triangle, not a capability.
% Measured over one endurance lap, F_braking returned a POSITIVE number 22112
% times -- an acceleration, out of the braking table, peaking at +150 m/s^2 --
% and returned decelerations beyond the car's own solved peak 12505 times,
% worst -44 m/s^2 against a peak of -23. Those fed straight into the profile.
%
% Projection is the honest reading of an infeasible query. If the solver asks
% for more lateral than the car can make at that speed then the car is already
% at the limit, so hand back the longitudinal capability at the most lateral it
% CAN make -- the end of that velocity row. Velocity is clamped to the solved
% range for the same reason; note gg2 only solves up to max_vel-0.5, so the top
% 0.5 m/s of the speed range has no samples at all and always clamps.
%
% 'nearest' extrapolation is the second line of defence, for what the clamp
% cannot reach. The cloud is not convex -- the per-velocity max lateral
% reverses five times in each table, worst 0.58 m/s^2 between 21 and 22 m/s,
% and 26 of 1160 grid nodes are missing to non-convergence -- so a clamped
% query can still land in a notch or hole the triangulation bridges. Nearest at
% least bounds the answer by a value the g-g actually solved for.
%
% Clamping to the SAMPLED envelope rather than the convex hull is deliberate:
% the hull bulges above the samples at 14 of 29 velocities, by up to 0.07 g.
%
% THIS CHANGES LAP TIMES. Strictly inside the sampled envelope it is
% bit-for-bit the old answer -- verified zero difference over 8000 interior
% samples -- so only states that were already ill-posed move.
F_grid = scatteredInterpolant([lat vel],long,'linear','nearest');

vel_rows = unique(vel);
lat_lo = zeros(size(vel_rows));
lat_hi = zeros(size(vel_rows));
for i = 1:numel(vel_rows)
    row = (vel == vel_rows(i));
    lat_lo(i) = min(lat(row));
    lat_hi(i) = max(lat(row));
end

F = @(lat_q,vel_q) clamped_query(F_grid,vel_rows,lat_lo,lat_hi,lat_q,vel_q);

end


function z = clamped_query(F_grid,vel_rows,lat_lo,lat_hi,lat_q,vel_q)
% abs() because the lookups hold only the positive-lateral half of the g-g
% (makeGG keeps p1 and p3, not the mirrored p2/p4). The solver already passes
% v^2*abs(kappa), but a signed query would otherwise be extrapolated into empty
% space rather than reflected.
%
% Scalar queries only: lininterp1 uses find(...,1,'last'), so a vectorised
% vel_q would silently take one bracket for the whole vector. Every caller is
% scalar -- Events2 and Events3 both call this once per track sample.
vel_q = min(max(vel_q,vel_rows(1)),vel_rows(end));
lat_q = min(max(abs(lat_q),lininterp1(vel_rows,lat_lo,vel_q)), ...
                           lininterp1(vel_rows,lat_hi,vel_q));
z = F_grid(lat_q,vel_q);
end

