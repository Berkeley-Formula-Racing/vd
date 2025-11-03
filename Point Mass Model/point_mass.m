clear;

%% Car params
m = 168 + 80; %car + driver
cla = 3.97;
cda = 1.48;
cop = 0.418;
torque_fn = KTM450(); 
max_braking_torque = 850; %N
rolling_radius = 0.1956;
wheel_radius = rolling_radius;
weight_r = 0.512;


%% constants
rolling_resistance = 0.025;
rho = 1.2;
g = 9.81;

%% tire forces

track_friction_scale = 0.55;
mu = @(Fz) track_friction_scale*2.727 .* (Fz./100) .^ (-0.096);

wf_static = 1 - weight_r;      % front static weight fraction
bf_aero   = cop;               % front aero DF fraction (you set cop as % front)
br_aero   = 1 - bf_aero;       % rear aero DF fraction

DF_total  = @(v) 0.5*rho*cla .* (v.^2);   % total downforce [N]

Fz_front_ax = @(v) (wf_static*m*g) + bf_aero*DF_total(v);   % front axle normal [N]
Fz_rear_ax  = @(v) (weight_r*m*g) + br_aero*DF_total(v);    % rear axle normal [N]

% per-tire normals
Fz_front_t  = @(v) 0.5 * Fz_front_ax(v);
Fz_rear_t   = @(v) 0.5 * Fz_rear_ax(v);

% per-axle longitudinal capacities (μ depends on per-tire load in lbf)
mu_front    = @(v) mu( Fz_front_t(v) / 4.4482216153 );
mu_rear     = @(v) mu( Fz_rear_t(v)  / 4.4482216153 );

Fcap_front_ax = @(v) 2 .* mu_front(v) .* Fz_front_t(v);     % front axle Fx cap [N]
Fcap_rear_ax  = @(v) 2 .* mu_rear(v)  .* Fz_rear_t(v);      % rear axle Fx cap [N]
Fcap_total_ax = @(v) Fcap_front_ax(v) + Fcap_rear_ax(v);    % four-tire Fx cap [N]

% forces
Fz_fun        = @(v) m*g + 0.5*rho*cla .* (v.^2);
Fz_tire_lbf   = @(v) (Fz_fun(v)./4) / 4.4482216153;
Ft_fun        = @(v) mu(Fz_tire_lbf(v)) .* Fz_fun(v);
Fx_fun        = @(v) 0.5*cda*rho .* (v.^2) + rolling_resistance .* Fz_fun(v);
Fy_fun        = @(v,r) m .* (v.^2) ./ r;
gfun          = @(v,r) (Fy_fun(v,r).^2) - (Ft_fun(v).^2) + (Fx_fun(v).^2);

%% powertrain
gears        = [32/16 30/18 28/20 26/22 24/24];  % gear ratios
final_drive  = 35/11; % final drive ratio
eta_driveline= 0.97; % overall driveline efficiency
rpm_redline  = 11500; % engine redline
rpm_idle     = 2000; % below this we say 0 torque
twopi_over_60 = 2*pi/60;
RPM = torque_fn(1,:).';
TQ = torque_fn(2,:).';
[RPM, iu] = unique(RPM,'stable');  TQ = TQ(iu);
[RPM, is] = sort(RPM);             TQ = TQ(is);

% lb-ft -> N·m :
TQ = TQ * 1.3558179483314004;

F = griddedInterpolant(RPM, TQ, 'pchip', 'linear');  % linear extrapolation
torque_fn = @(rpm) max(0, F(rpm));  


%wheel force at gear and vel
engine_force_for_gear = @(v,gr) ...
    max(0, torque_fn( min(rpm_redline, max(rpm_idle, (v(:)./wheel_radius).*gr.*final_drive*60/(2*pi))) ) ...
        .* gr .* final_drive .* eta_driveline ./ wheel_radius );
%best wheel force for vel
engine_force = @(v) ...
    max( cell2mat( arrayfun(@(gr) engine_force_for_gear(v,gr), gears, 'UniformOutput', false) ), [], 2 );

%% Straight-line acceleration / braking capacities 
% Available tire longitudinal capacity on a straight (Fy=0) is Ft_fun(v).
% Net forward force = min(engine_force, Ft_fun) - aero/rolling drag.
% Net braking force  = min(brake_cap, Ft_fun) + aero/rolling drag (drag helps).

% If you want to cap brake torque, set a finite value here
brake_is_tire_limited = false;
if brake_is_tire_limited
    brake_cap_force = Inf;   
else
    brake_cap_force = max_braking_torque / wheel_radius;
end

engine_force_RWD = @(v) min(engine_force(v), Fcap_rear_ax(v) );

amax_fun = @(v) max(0, engine_force_RWD(v) - Fx_fun(v)) ./ m;

abrk_fun = @(v) max(0, min(brake_cap_force,  Ft_fun(v)) + Fx_fun(v) ) ./ m;  %magnitude

%% Build lookup tables (velocity vs accel)
v_grid = linspace(0, 33, 901)';
accel_lookup = amax_fun(v_grid);
decel_lookup = abrk_fun(v_grid); %magnitude - will need to use a >= 0 in backwards pass

accel_F = griddedInterpolant(v_grid, accel_lookup, 'linear', 'nearest');
decel_F = griddedInterpolant(v_grid, decel_lookup, 'linear', 'nearest');

%% track solving
load('michigantrack2024.mat');
autocross_track = [arclength; curvature];

% Solve gfun(v, r_i) = 0 for each track point to get lateral limit w/ aero & mu(v)
s   = arclength(:);
k   = curvature(:);
r_i = 1 ./ max(1e-9, abs(k));  %avoid div/0

v_latlim = zeros(size(r_i));
v_guess  = 25;

for i = 1:numel(r_i)
    r = r_i(i);
    a = 0; b = v_guess;
    gi0 = gfun(a, r); gib = gfun(b, r);
    tries = 0;
    while gi0 * gib > 0 && tries < 15
        b = b * 2.0;
        gib = gfun(b, r);
        tries = tries + 1;
    end
    if gi0 * gib > 0
        % fallback: no real root found -> use large speed (straight), limited by accel/brake later
        v_apex = b;
    else
        v_apex = fzero(@(v) gfun(v, r), [a, b]);
    end
    if ~isfinite(v_apex) || v_apex <= 0, v_apex = b; end
    v_latlim(i) = v_apex;
    v_guess = v_apex;     % carry forward a reasonable next guess
end

% If curvature==0 (very large r), v_latlim might be huge; cap to grid max for sanity
v_latlim = min(v_latlim, v_grid(end));

%% Forward/Backward pass along arclength
% Use v_{k+1} = sqrt( v_k^2 + 2*a(v_k)*Δs ) for forward (accel)
% and      v_k = min( v_k, sqrt( v_{k+1}^2 + 2*a_brake(v_{k+1})*Δs ) ) for backward (braking)
N  = numel(s);
ds = diff(s); %segment lengths
if size(ds,1) < size(ds,2), ds = ds'; end

v_profile = zeros(N,1);
v_profile(1) = min(1.0, v_latlim(1)); % small nonzero start

% Forward accel
for i = 1:N-1
    a_here = accel_F( max(0, v_profile(i)) );
    v_next = sqrt( max(0, v_profile(i)^2 + 2*a_here*ds(i)) );
    v_profile(i+1) = min(v_next, v_latlim(i+1));
end

% Backward braking
for i = N-1:-1:1
    ab_here = decel_F( max(0, v_profile(i+1)) );  % magnitude
    v_cap   = sqrt( max(0, v_profile(i+1)^2 + 2*ab_here*ds(i)) );
    v_profile(i) = min( v_profile(i), v_cap );
    v_profile(i) = min( v_profile(i), v_latlim(i) );
end

%% compute lap time and quick plot
% Integrate dt = ds / v over the track
epsv = 0.2;
dt_segments = [ds ./ max(v_profile(1:end-1), epsv); 0];  
lap_time = sum(dt_segments);

fprintf('Estimated lap time: %.2f s\n', lap_time);

figure; 
subplot(3,1,1); plot(s, v_profile, 'LineWidth', 1.2); grid on; ylabel('v [m/s]'); title('Speed Profile');
subplot(3,1,2); plot(v_grid, accel_lookup, 'LineWidth', 1.2); grid on; ylabel('a_{accel} [m/s^2]'); xlabel('v [m/s]'); title('Accel lookup');
subplot(3,1,3); plot(v_grid, decel_lookup, 'LineWidth', 1.2); grid on; ylabel('a_{brake} [m/s^2]'); xlabel('v [m/s]'); title('Decel lookup (magnitude)');