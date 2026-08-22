function [x_accel,long_accel,long_accel_guess] = max_long_accel_cornering(long_vel_guess,lat_accel_value,car,x0)
% uses fmincon to maximize longitudinal acceleration at a fixed longitudinal
% velocity and a fixed lateral acceleration
%
% inputs  long_vel_guess:  the velocity to solve at (it is fixed, not a guess
%                          -- the name is historical)
%         lat_accel_value: the lateral acceleration to hold, m/s^2
%         car:             Car object
%         x0:              optional starting point. Supplying one REPLACES the
%                          default start; the analytic start below is always
%                          tried as well.
%
% SOLVED FROM TWO STARTS, NOT ONE -- the same fix as max_lat_accel, for the
% same reason.
%
% The default start scaled the yaw rate with the requested lateral
% acceleration and left everything else at a straight-line value: steer 0,
% lateral velocity 0. That is precisely the mixture max_lat_accel warns
% against. At 7 m/s and 5.9 m/s^2 the start asks the car to yaw at 0.84 rad/s
% with the wheels pointed straight ahead, which puts the front tyres at +5.5
% deg of slip and the rears at -5.3: the two axles make lateral force in
% OPPOSITE directions, so the start is nowhere near the steady-state manifold
% the constraints define and fmincon spends the whole budget looking for it.
%
% Measured on the calibration car (grip 0.60) over the 29 x 20 grid gg2 uses:
%
%   start                       failures   median evals   p99 evals   max
%   fixed (steer 0, lat vel 0)     17/580            90        1006  1009
%   analytic (below)                0/580            82         326   867
%
% All 17 failures burned the 1000-evaluation cap or gave up (16 exitflag 0,
% 1 exitflag -2) and all sit between 6 and 13 m/s, where the yaw rate the
% constraint demands is largest. None of them is infeasible: the analytic
% start converges at every one, worst residual 3.2e-6.
%
% Taking the better of the two also changes the answer at 109 of the 562
% points that already converged, but only 1 of those moves by more than
% 0.1 m/s^2 (by 1.19). The recovered cells are the point here, not the
% accuracy of the ones that already worked.
%
% Across the 38 cars in the grip and aero sweep caches (1872 nodes) the same
% start takes the failures from 7 to 2 and improves 308 points by up to 1.83
% m/s^2, with none worse.
%
% Selection is by FEASIBILITY first and objective value only as a tie-break,
% the same rule max_lat_accel and vel_cornering_sweep use -- a solve that has
% drifted through the constraint boundary reports MORE acceleration while
% being less feasible, so ranking on value alone rewards the wrong answer.

FEAS_TOL = 1e-4;

% The lateral constraint pins the yaw rate outright (v*r = a_y, and the bounds
% below fix it), so it is the one part of the old start that was already
% right. Everything else has to be built to match it.
yaw_rate = lat_accel_value/long_vel_guess;

seeds = {};
if nargin >= 4 && ~isempty(x0)
    x0(3) = long_vel_guess;
    x0(5) = yaw_rate;
    seeds{end+1} = x0;
else
    seeds{end+1} = [0, 0, long_vel_guess, 0, yaw_rate, 0, 0, 0.01, 0.01];
end

% Analytic start: the steady turn that actually produces this yaw rate.
% Assume equal slip angle at both axles, which makes the steer angle the
% Ackermann value W_b*r/v and leaves the body slip to follow from the rear
% axle, lat_vel = l_r*r (from alpha_r = (lat_vel - l_r*r)/v). Both come from
% the same yaw rate, so the three agree with each other -- that is the whole
% point. Clamped to the bounds below because at 5 m/s the Ackermann angle
% runs past the 25 deg steering limit.
steer_seed   = min(max(rad2deg(car.W_b*yaw_rate/long_vel_guess),0),25);
lat_vel_seed = min(max(car.l_r*yaw_rate,-3),3);
seeds{end+1} = [steer_seed, 1, long_vel_guess, lat_vel_seed, yaw_rate, 0, 0, 0.02, 0.02];

% bounds
steer_angle_bounds = [0,25];
throttle_bounds = [0,1];
long_vel_bounds = [long_vel_guess,long_vel_guess];
lat_vel_bounds = [-3,3];
yaw_rate_bounds = [yaw_rate,yaw_rate];
kappa_1_bounds = [0,0];
kappa_2_bounds = [0,0];
kappa_3_bounds = [0,0.2];
kappa_4_bounds = [0,0.2];

A = [];
b = [];
Aeq = [0 0 1 0 0 0 0 0 0
       0 0 0 0 0 1 0 0 0
       0 0 0 0 0 0 1 0 0];
beq = [long_vel_guess 0 0];
lb = [steer_angle_bounds(1),throttle_bounds(1),long_vel_bounds(1),lat_vel_bounds(1),...
    yaw_rate_bounds(1),kappa_1_bounds(1),kappa_2_bounds(1),kappa_3_bounds(1),kappa_4_bounds(1)];
ub = [steer_angle_bounds(2),throttle_bounds(2),long_vel_bounds(2),lat_vel_bounds(2),...
    yaw_rate_bounds(2),kappa_1_bounds(2),kappa_2_bounds(2),kappa_3_bounds(2),kappa_4_bounds(2)];

% objective function: longitudinal acceleration (forwards)
f = @(P) -car.long_accel(P);

% constrained to lateral acceleration value
constraint = @(P) car.constraint4(P,lat_accel_value);

% default algorithm is interior-point
%
% Left at 1000. The old note here said raising it bought nothing because these
% solves genuinely diverge; that is no longer the whole story. The median solve
% takes 90 evaluations and the failures all burned the cap, so the cap was
% never the constraint -- the start was. From the analytic start the worst
% solve on the same grid takes 867 and none hits the cap.
options = setOptimoptions(1000);

x = []; exitflag = 0; bestVal = -inf; bestRes = inf;
for k = 1:numel(seeds)
    [xk,~,flagk] = fmincon(f,min(max(seeds{k},lb),ub),A,b,Aeq,beq,lb,ub,constraint,options);

    % one extra constraint evaluation per candidate, against a hundred or more
    % inside the solve -- the residual is what the selection turns on and
    % fmincon does not hand it back
    [ck,ceqk] = constraint(xk);
    resk = max([max(abs(ceqk)), max(ck), 0]);
    valk = car.long_accel(xk);

    if isempty(x) || preferNew(flagk,resk,valk,exitflag,bestRes,bestVal,FEAS_TOL)
        x = xk; exitflag = flagk; bestRes = resk; bestVal = valk;
    end
end

[engine_rpm,beta,lat_accel,long_accel,yaw_accel,wheel_accel,omega,current_gear,...
Fzvirtual,Fz,alpha,T] = car.equations(x);

long_accel_guess = x;

% generate vector of control variable values
% column 3 is the achieved lateral acceleration v*r, NOT the lat_accel output
% of equations() -- that one is the lateral residual the constraint drives to
% zero, and storing it here made every logged point read as gLat = 0
x_accel = [exitflag long_accel x(3)*x(5) x omega(1:4) engine_rpm current_gear beta...
    Fz(1:4) alpha(1:4) T(1:4)];

end


function take = preferNew(flag,res,val,bFlag,bRes,bVal,tol)
% Feasibility first, objective value only between two solves that both cleared
% it. Ranking on value alone is what lets a solve that has crept through the
% constraint boundary win, since those report MORE acceleration, not less;
% ranking on residual alone would accept a trivially feasible but badly
% suboptimal point. Same rule as max_lat_accel.
okNew = (flag == 1) || (flag == 2);
okOld = (bFlag == 1) || (bFlag == 2);
if okNew ~= okOld, take = okNew; return, end     % a converged solve always wins
if ~okNew,         take = res < bRes; return, end % neither: less infeasible

goodNew = res <= tol;
goodOld = bRes <= tol;
if goodNew && goodOld
    take = val > bVal;                            % both legal -> more accel wins
elseif goodNew ~= goodOld
    take = goodNew;                               % only one is legal
else
    take = res < bRes;                            % neither legal -> less illegal
end
end
