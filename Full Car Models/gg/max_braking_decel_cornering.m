function [x_braking,long_decel,braking_decel_guess] = max_braking_decel_cornering(long_vel_guess,lat_accel_value,car,x0)
% uses fmincon to maximize braking deceleration at a fixed longitudinal
% velocity and a fixed lateral acceleration
%
% inputs  long_vel_guess:  the velocity to solve at (fixed, not a guess)
%         lat_accel_value: the lateral acceleration to hold, m/s^2
%         car:             Car object
%         x0:              optional starting point. Supplying one REPLACES the
%                          default start; the analytic starts below are always
%                          tried as well.
%
% SOLVED FROM THREE STARTS -- see max_lat_accel for the pattern and
% max_long_accel_cornering for the accel half of the same fix.
%
% The default start describes a car going nearly straight: yaw rate 0.01 rad/s
% against the ay/v the lateral constraint demands, which at 5 m/s is 1.5. It is
% at least self-consistent, which is why it failed less often than the accel
% start did, but it leaves the whole turn to be discovered by the optimiser and
% it lands on a poor local optimum over most of the top of the speed range.
%
% Measured on the calibration car (grip 0.60) over the 29 x 20 grid gg2 uses:
%
%   starts                                failures   better   worst under-report
%   default only                             9/580        -                    -
%   default + brakes-off analytic            7/580      175        7.59 m/s^2
%   default + light-brake analytic           5/580       80        0.62 m/s^2
%   all three                                5/580      211        7.59 m/s^2
%
% "better" counts points that already converged and where a different start
% found MORE braking at the same lateral acceleration and a smaller constraint
% residual. The large numbers are not noise: at 30 m/s and 18.55 m/s^2 lateral
% the default start returns 6.75 m/s^2 of deceleration with the wheel straight,
% and the brakes-off start returns 14.33 with 11.7 deg of steer -- both are
% legal steady states (residual 3.5e-9 and 5.7e-14), but the second uses the
% rearward component of the front tyre lateral force and the first does not.
% The g-g has been under-reporting braking capability at high speed, which
% matters more than the dropped cells.
%
% The two analytic starts differ ONLY in whether the brakes are on, and they
% find different optima: brakes-off finds the steered branch, light-brake is
% the robust one (8 failures on its own against 75). Neither alone reaches the
% floor of 5. Dropping the original start instead of keeping it costs 56 points
% that get WORSE, which is why all three are kept.
%
% Across the 38 cars in the grip and aero sweep caches (1872 nodes) the
% brakes-off start alone takes the failures from 127 to 73 and improves 476
% points by up to 10.47 m/s^2, with none worse. The default start fails 6.8%
% of the time on that spread, four times its rate on the calibration car.
%
% 5 of the 9 original failures are not solver failures at all and no start
% fixes them -- see the note on the steering bound below.
%
% Selection is by FEASIBILITY first and objective value only as a tie-break;
% exitflag alone is not enough, since fmincon's stopping test is relative and
% does return exitflag 1 at points whose worst residual is order 1.

FEAS_TOL = 1e-4;

% The lateral constraint requires v*r = a_y, so this is the yaw rate every
% feasible point has. Unlike the accel solver the bounds do not pin it, which
% is why the old start could sit at 0.01 and still be accepted as a start.
yaw_rate = lat_accel_value/long_vel_guess;

seeds = {};
if nargin >= 4 && ~isempty(x0)
    x0(3) = long_vel_guess;
    seeds{end+1} = x0;
else
    seeds{end+1} = [0.01,-0.1,long_vel_guess,-0.1,0.01,-0.01,-0.01,-0.01,-0.01];
end

% Analytic starts: the steady turn that actually produces this yaw rate.
% Assume the same slip angle at both axles, which makes the steer angle the
% Ackermann value W_b*r/v and leaves the body slip to follow from the rear
% axle, lat_vel = v*alpha + l_r*r. The slip angle itself is the one number
% that has to be estimated: 0.35 deg per m/s^2 of lateral acceleration, read
% off the converged solves (1.1 deg at 3.4 m/s^2, 3.0 at 8.4, 5.2 at 14.7),
% capped at 8 deg where the tyre is past its peak and a larger angle is a
% worse guess rather than a better one. Left as a constant deliberately:
% dividing it by the grip scaling, which is what the linear tyre model
% suggests, is WORSE across the 38 grip and aero sweep cars -- 94 failures of
% 1872 against 73, and half as many improved points.
alpha_seed   = -min(0.35*lat_accel_value, 8)*pi/180;
steer_seed   = min(max(rad2deg(car.W_b*yaw_rate/long_vel_guess),0),22);
lat_vel_seed = min(max(long_vel_guess*alpha_seed + car.l_r*yaw_rate,-3),3);
% brakes off: satisfies the four wheel torque balances exactly at kappa = 0,
% and is the only one of the two that finds the steered braking branch
seeds{end+1} = [steer_seed, 0, long_vel_guess, lat_vel_seed, min(yaw_rate,2), 0, 0, 0, 0];
% light brake: 0.05 of the pedal against a slip ratio that roughly matches it,
% so the wheel equations start near balance rather than 100 N-m out, which is
% where the old start put them
seeds{end+1} = [steer_seed, -0.05, long_vel_guess, lat_vel_seed, min(yaw_rate,2), ...
    -0.002, -0.002, -0.002, -0.002];

% bounds
%
% NOTE the steering limit is 22 here and 25 in max_lat_accel and
% max_long_accel_cornering. gg2 sizes its lateral grid from max_lat_accel, so
% at low speed it asks this solver for lateral accelerations it cannot reach:
% capped at 22 deg the car makes 6.66 m/s^2 at 5 m/s and 9.51 at 6, against
% the 7.59 and 10.65 gg2 asks for. Those cells are infeasible, not unsolved,
% and they are 5 of the 9 failures. Raising this to 25 or lowering gg2's
% request is a separate decision from the starts.
steer_angle_bounds = [0,22];
throttle_bounds = [-1,0];
long_vel_bounds = [long_vel_guess-0.1,long_vel_guess+0.1];
lat_vel_bounds = [-3,3];
yaw_rate_bounds = [0,2];
kappa_1_bounds = [-0.2,0];
kappa_2_bounds = [-0.2,0];
kappa_3_bounds = [-0.2,0];
kappa_4_bounds = [-0.2,0];

A = [];
b = [];
Aeq = [0 0 1 0 0 0 0 0 0];
beq = [long_vel_guess];
lb = [steer_angle_bounds(1),throttle_bounds(1),long_vel_bounds(1),lat_vel_bounds(1),...
    yaw_rate_bounds(1),kappa_1_bounds(1),kappa_2_bounds(1),kappa_3_bounds(1),kappa_4_bounds(1)];
ub = [steer_angle_bounds(2),throttle_bounds(2),long_vel_bounds(2),lat_vel_bounds(2),...
    yaw_rate_bounds(2),kappa_1_bounds(2),kappa_2_bounds(2),kappa_3_bounds(2),kappa_4_bounds(2)];

% objective function: longitudinal acceleration (backwards)
f = @(P) car.long_accel(P);
% constrained to lateral acceleration value
constraint = @(P) car.constraint4(P,lat_accel_value);
% left at 1000 -- see the note in max_long_accel_cornering
opts = setOptimoptions(1000);

x = []; exitflag = 0; bestVal = -inf; bestRes = inf;
for k = 1:numel(seeds)
    [xk,~,flagk] = fmincon(f,min(max(seeds{k},lb),ub),A,b,Aeq,beq,lb,ub,constraint,opts);

    % one extra constraint evaluation per candidate, against a hundred or more
    % inside the solve -- the residual is what the selection turns on and
    % fmincon does not hand it back
    [ck,ceqk] = constraint(xk);
    resk = max([max(abs(ceqk)), max(ck), 0]);
    valk = -car.long_accel(xk);    % more braking is a bigger number

    if isempty(x) || preferNew(flagk,resk,valk,exitflag,bestRes,bestVal,FEAS_TOL)
        x = xk; exitflag = flagk; bestRes = resk; bestVal = valk;
    end
end

[engine_rpm,beta,lat_accel,long_accel,yaw_accel,wheel_accel,omega,current_gear,...
Fzvirtual,Fz,alpha,T] = car.equations(x);

braking_decel_guess = x;

% generate vector of control variable values
% column 3 is the achieved lateral acceleration v*r -- see the note in
% max_long_accel_cornering.m
x_braking = [exitflag long_accel x(3)*x(5) x omega(1:4) engine_rpm current_gear beta...
    Fz(1:4) alpha(1:4) T(1:4)];

long_decel = long_accel; %make this line less confusing

end


function take = preferNew(flag,res,val,bFlag,bRes,bVal,tol)
% Feasibility first, objective value only between two solves that both cleared
% it. Ranking on value alone is what lets a solve that has crept through the
% constraint boundary win, since those report MORE braking, not less; ranking
% on residual alone would accept a trivially feasible but badly suboptimal
% point. Same rule as max_lat_accel.
okNew = (flag == 1) || (flag == 2);
okOld = (bFlag == 1) || (bFlag == 2);
if okNew ~= okOld, take = okNew; return, end     % a converged solve always wins
if ~okNew,         take = res < bRes; return, end % neither: less infeasible

goodNew = res <= tol;
goodOld = bRes <= tol;
if goodNew && goodOld
    take = val > bVal;                            % both legal -> more braking wins
elseif goodNew ~= goodOld
    take = goodNew;                               % only one is legal
else
    take = res < bRes;                            % neither legal -> less illegal
end
end
