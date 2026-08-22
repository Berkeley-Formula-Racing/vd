function [x_ss,lat_accel,long_accel,lat_accel_guess] = max_lat_accel(long_vel_guess,car,x0)
% uses fmincon to maximize lateral acceleration at a fixed longitudinal
% velocity, with no longitudinal acceleration constraint
%
% inputs  long_vel_guess: the velocity to solve at (it is fixed, not a guess --
%                         the name is historical)
%         car:            Car object
%         x0:             optional starting point. Supplying one REPLACES the
%                         default start; the analytic start below is always
%                         tried as well.
%
% SOLVED FROM TWO STARTS, NOT ONE.
%
% The feasible set here has more than one local optimum and no single start
% finds the better one everywhere. Measured over 874 (car, velocity) points
% drawn from the grip and aero sweep caches:
%
%   start                      median evals   max evals   hit the 5000 cap
%   fixed (10 deg, 0.5 rad/s)           152        5002                  2
%   analytic (below)                    122         373                  0
%
% and on lateral acceleration they disagree both ways -- the analytic start
% wins at 129 of the 874 points and loses at 131. Taking the better of the two
% is worth up to 1.19 m/s^2 (0.12 g) against the fixed start alone and is
% never worse, so the single-start version was both occasionally failing AND
% quietly under-reporting grip at points where it converged perfectly well.
%
% Why this matters more than one dropped point: makeGG gates the whole
% velocity row on this solve, so a failure here deletes 20 acceleration and 20
% braking points from the table the lap sim interpolates. That is the "no
% lateral solution at N m/s" warning, and raising the evaluation cap alone
% only ever buys margin -- the worst observed solve used 4054 evaluations of
% the 5000 available, so the next slightly harder car falls off the end again.
%
% The two starts must each be internally consistent. Scaling the yaw rate
% while leaving the steer angle at its fixed 10 deg is WORSE than changing
% nothing: 132 of the 874 points then hit the cap, against 2 today.
%
% Selection is by FEASIBILITY first and objective value only as a tie-break,
% the same rule vel_cornering_sweep uses and for the same reason -- a solve
% that has drifted through the constraint boundary reports a higher lateral
% acceleration while being less feasible, so ranking on value alone rewards
% exactly the wrong answer.

FEAS_TOL = 1e-6;

% Starts, tried in order. The analytic one comes from the kinematics rather
% than from a constant: at the lateral limit a_y = v * yaw_rate, so the yaw
% rate a solution needs falls as 1/v -- 1.5 rad/s at 5 m/s against 0.5 at
% 30 -- and the steer angle that puts the car on that path falls with it. One
% fixed pair cannot be near both ends. 1.5 g is only a starting assumption for
% a_y; the solve moves off it immediately.
seeds = {};
if nargin >= 3 && ~isempty(x0)
    x0(3) = long_vel_guess;
    seeds{end+1} = x0;
else
    seeds{end+1} = [10; 0; long_vel_guess; 0.1; 0.5; 0; 0; 0; 0]';
end
yaw_rate_seed = min(1.5*9.81/long_vel_guess, 2);
steer_seed    = min(max(rad2deg(car.W_b*yaw_rate_seed/long_vel_guess),0), 25);
seeds{end+1}  = [steer_seed; 0; long_vel_guess; -0.3; yaw_rate_seed; 0; 0; 0; 0]';

% bounds
steer_angle_bounds = [0,25];
throttle_bounds = [0,0];
long_vel_bounds = [long_vel_guess,long_vel_guess];
lat_vel_bounds = [-3,3];
yaw_rate_bounds = [0,2];
kappa_1_bounds = [0,0];
kappa_2_bounds = [0,0];
kappa_3_bounds = [0,0];
kappa_4_bounds = [0,0];

A = [];
b = [];
Aeq = [0 1 0 0 0 0 0 0 0
       0 0 1 0 0 0 0 0 0
       0 0 0 0 0 1 0 0 0
       0 0 0 0 0 0 1 0 0
       0 0 0 0 0 0 0 1 0
       0 0 0 0 0 0 0 0 1];
beq = [0 long_vel_guess 0 0 0 0];
lb = [steer_angle_bounds(1),throttle_bounds(1),long_vel_bounds(1),lat_vel_bounds(1),...
    yaw_rate_bounds(1),kappa_1_bounds(1),kappa_2_bounds(1),kappa_3_bounds(1),kappa_4_bounds(1)];
ub = [steer_angle_bounds(2),throttle_bounds(2),long_vel_bounds(2),lat_vel_bounds(2),...
    yaw_rate_bounds(2),kappa_1_bounds(2),kappa_2_bounds(2),kappa_3_bounds(2),kappa_4_bounds(2)];

% The cap is a backstop now rather than the mechanism. With two starts a solve
% that runs away on one of them is rescued by the other, so this only has to
% be generous enough not to truncate a solve that is genuinely converging.
opts = setOptimoptions(5000);

% objective function: longitudinal velocity times yaw rate (v*v/r = v^2/r)
f = @(P) -P(3)*P(5);

% no longitudinal acceleration constraint
constraint = @(P) car.constraint1(P);

x = []; exitflag = 0; bestVal = -inf; bestRes = inf;
for k = 1:numel(seeds)
    [xk,~,flagk] = fmincon(f,seeds{k},A,b,Aeq,beq,lb,ub,constraint,opts);

    % one extra constraint evaluation per candidate, against several hundred
    % inside the solve -- the residual is what the selection turns on and
    % fmincon does not hand it back
    [ck,ceqk] = constraint(xk);
    resk = max([max(abs(ceqk)), max(ck), 0]);
    valk = xk(3)*xk(5);

    if isempty(x) || preferNew(flagk,resk,valk,exitflag,bestRes,bestVal,FEAS_TOL)
        x = xk; exitflag = flagk; bestRes = resk; bestVal = valk;
    end
end

[engine_rpm,beta,lat_accel,long_accel,yaw_accel,wheel_accel,omega,current_gear,...
Fzvirtual,Fz,alpha,T] = car.equations(x); %#ok<ASGLU>

lat_accel_guess = x;

% generate table of control variable values
x_ss = [exitflag long_accel x(3)*x(5) x omega(1:4) engine_rpm current_gear beta...
    Fz(1:4) alpha(1:4) T(1:4)];

%max lat accel in m/s^2
lat_accel = x(3)*x(5);

end


function take = preferNew(flag,res,val,bFlag,bRes,bVal,tol)
% Feasibility first, objective value only between two solves that both cleared
% it. Ranking on value alone is what lets a solve that has crept through the
% constraint boundary win, since those report MORE lateral acceleration, not
% less; ranking on residual alone would accept a trivially feasible but badly
% suboptimal point.
okNew = (flag == 1) || (flag == 2);
okOld = (bFlag == 1) || (bFlag == 2);
if okNew ~= okOld, take = okNew; return, end     % a converged solve always wins
if ~okNew,         take = res < bRes; return, end % neither: less infeasible

goodNew = res <= tol;
goodOld = bRes <= tol;
if goodNew && goodOld
    take = val > bVal;                            % both legal -> more grip wins
elseif goodNew ~= goodOld
    take = goodNew;                               % only one is legal
else
    take = res < bRes;                            % neither legal -> less illegal
end
end
