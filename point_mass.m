clear;

%Car params
m = 168 + 80; %car + driver
cla = 3.97;
cda = 1.48;
torque_fn = KTM450();
max_braking_torque = 1;
mu = @(Fz)2.727 * (Fz/100)^-0.096;
rolling_radius = 0.1956;
%track params
track_mu = 1.4;
rolling_resistance = 0.025;
rho = 1.2;
g = 9.81;

load('michigantrack2024.mat');
autocross_track = [arclength; curvature];

[extrema,extrema_indices] = curvature_apexes(arclength,curvature);

arclength = [0 arclength];

extrema_radii = 1./extrema;

apex_velocity = zeros(size(extrema_radii));

Fz_fun  = @(v) m*g + 0.5*rho*cla*v^2;   
Fz_tire_lbf = @(v) (Fz_fun(v)/4) / 4.4482216153;
Ft_fun  = @(v) mu(Fz_tire_lbf(v)).*Fz_fun(v);       
Fx_fun  = @(v) 0.5*cda*rho*v^2 + rolling_resistance*Fz_fun(v); 
Fy_fun = @(v, r) m*v^2 / r;
gfun = @(v, r) Fy_fun(v, r)^2 - Ft_fun(v)^2 + Fx_fun(v)^2;

vmax_guess = 25;

for i = 1:numel(extrema_radii)
    r = extrema_radii(i);
    a = 0;
    b = vmax_guess;
    gi0 = gfun(a, r);
    gib = gfun(b, r);
    tries = 0;
    while gi0 * gib > 0 && tries < 15
        b = b * 2.5;
        gib = gfun(b, r);
        tries = tries + 1;
    end
    
    % Solve
    v_apex = fzero(@(v) gfun(v, r), [a,b]);
    
    % Numerical safety
    if ~isfinite(v_apex) || v_apex < 0
        v_apex = 0;
    end
    apex_velocity(i) = v_apex;
end



