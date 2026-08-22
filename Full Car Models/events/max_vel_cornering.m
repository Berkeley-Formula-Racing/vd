function [x_corner_vel,max_vel_corner,vel_corner_guess,ceqMax] = max_vel_cornering(radius,max_vel,car,x0)
% uses fmincon to minimize the objective function subject to constraints
% optimizes lateral acceleration with no velocity constraint and no
%   longitudinal acceleration constraint (in contrast to skidpad)

% isempty as well as nargin: the sweep hands back [] when a solve did not
% converge, meaning "I have no usable seed, use your own"
if nargin < 4 || isempty(x0) % no initial guess supplied
    % initial guesses
    steer_angle_guess = 5; % degrees
    throttle_guess = 0;
    long_vel_guess = sqrt(9.81*0.3*radius); % approximation method to help convergence
    lat_vel_guess = 0.2; 
    yaw_rate_guess = 0.3;

    kappa_1_guess = 0;
    kappa_2_guess = 0;
    kappa_3_guess = 0;
    kappa_4_guess = 0;
    
    x0 = [steer_angle_guess,throttle_guess,long_vel_guess,lat_vel_guess,yaw_rate_guess,kappa_1_guess,...
        kappa_2_guess,kappa_3_guess,kappa_4_guess];
end

% bounds
steer_angle_bounds = [0,22];
throttle_bounds = [0,0]; 
long_vel_bounds = [0,max_vel];
lat_vel_bounds = [-3,3];
yaw_rate_bounds = [0,2];
kappa_1_bounds = [0,0];
kappa_2_bounds = [0,0];
kappa_3_bounds = [0,0];
kappa_4_bounds = [0,0];

A = [];
b = [];
Aeq = [0 1 0 0 0 0 0 0 0
       0 0 0 0 0 1 0 0 0
       0 0 0 0 0 0 1 0 0
       0 0 0 0 0 0 0 1 0
       0 0 0 0 0 0 0 0 1];
beq = [0 0 0 0 0];
lb = [steer_angle_bounds(1),throttle_bounds(1),long_vel_bounds(1),lat_vel_bounds(1),...
    yaw_rate_bounds(1),kappa_1_bounds(1),kappa_2_bounds(1),kappa_3_bounds(1),kappa_4_bounds(1)];
ub = [steer_angle_bounds(2),throttle_bounds(2),long_vel_bounds(2),lat_vel_bounds(2),...
    yaw_rate_bounds(2),kappa_1_bounds(2),kappa_2_bounds(2),kappa_3_bounds(2),kappa_4_bounds(2)];

% default algorithm is interior-point
% Left at the 1e-2 default. Tightening this was the obvious suspect for the
% sweep's sensitivity, but measured over radii 4.5-60 the results are
% BIT-IDENTICAL from 1e-2 down to 1e-6 -- the constraint tolerance is not
% what is binding. The instability comes from genuinely distinct local optima
% (at r = 20 the warm-started solve finds 18.628 m/s where the analytic seed
% finds 18.302), so which one you land on depends on where you start, not on
% how tightly you converge. The guards in vel_cornering_sweep address that;
% this does not.
opts = setOptimoptions(5000);
% objective function: longitudinal velocity times yaw rate (v*v/r = v^2/r)
f = @(P) -P(3)*P(5);                                     
% no longitidunal acceleration constraint
constraint = @(P) car.constraint5(P,radius);
% fval: objective function value (v^2/r)
[x,fval,exitflag] = fmincon(f,x0,A,b,Aeq,beq,lb,ub,constraint,opts);
% Only hand back a seed the caller can trust. vel_cornering_sweep uses this
% as the next radius's initial guess, so returning the raw x from a
% non-converged solve propagates a bad state down the rest of the sweep.
% Empty tells the caller to fall back to its own deterministic guess.
if exitflag == 1 || exitflag == 2
    vel_corner_guess = x;
else
    vel_corner_guess = [];
end

% Steady-state residual at the returned point. The caller selects between
% competing solves on this rather than on velocity: a warm-started chain that
% has drifted reports a HIGHER cornering speed than a fresh solve while
% sitting further outside the constraints -- measured, the drifting chain's
% worst point had max|ceq| = 2e-2, which does not even satisfy the 1e-2
% tolerance it was solved under, and was 0.36 m/s "faster" for it.
[~,ceq] = constraint(x);
ceqMax = max(abs(ceq));

[engine_rpm,beta,lat_accel,long_accel,yaw_accel,wheel_accel,omega,current_gear,...
Fzvirtual,Fz,alpha,T,gamma] = car.equations(x);  
gamma = transpose(gamma);
%disp(gamma)
max_vel_corner = x(3);

x_corner_vel = [exitflag long_accel x(3)*x(5) x omega(1:4) engine_rpm current_gear beta...
    Fz(1:4) alpha(1:4) T(1:4)];

