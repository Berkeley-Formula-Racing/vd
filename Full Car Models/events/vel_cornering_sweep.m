function [x_table_corner_vel,radius,max_vel_corner_vector] = vel_cornering_sweep(car)
% Max cornering velocity against radius, used by Autocross and Endurance to
% cap speed through each corner.
%
% WHY THIS IS MORE THAN A FOR LOOP
%
% Warm-starting each radius from the previous solution helps convergence, but
% it chains all 556 solves together: a difference at any link becomes the next
% initial guess and accumulates. Measured on the plain chained version, a 1e-9
% change in front tire force moved this vector by up to 0.59 m/s (6%) over one
% contiguous band, r = 13.6-29.9 m, and moved autocross by 0.11 s. That noise
% floor was larger than the CdA and CoP effects an aero sweep was trying to
% resolve, which is why those came out non-monotonic in the parameter.
%
% The drift is not the chain finding better optima. It creeps along and
% through the constraint boundary: the drifted chain reports HIGHER cornering
% speeds while being far less feasible, and its worst point had
% max|ceq| = 2e-2 -- outside the 1e-2 tolerance it was solved under. Selecting
% on speed therefore rewards exactly the wrong solve. This selects on
% FEASIBILITY first and uses speed only to break ties between solves that are
% both properly converged.
%
% Three mechanisms, in order of importance:
%
%   1. FEASIBILITY GATE. A candidate must converge and satisfy
%      max|ceq| <= FEAS_TOL to be preferred on speed. Calibrated from data:
%      exactly 1 of 556 chained points exceeds 1e-6, and it is the bad one.
%
%   2. PERIODIC RESTART. Every RESET_EVERY radii the solve is also attempted
%      from the deterministic analytic seed, which bounds how far the chain
%      can drift before being checked against a reproducible reference.
%      RESET_EVERY is a robustness knob, NOT a tuned constant -- the gate
%      decides what is accepted, so changing it trades a little cost for a
%      little drift and cannot change which answers are legal.
%
%   3. MONOTONICITY. A larger radius can never allow less speed, so a
%      violation means that solve fell off the branch and a fresh attempt is
%      made regardless of where the periodic restarts fall.
%
% A restart is only ever COMMITTED to if it converges. The analytic seed fails
% outright around r = 6 and r = 10, and an earlier version that restarted
% blindly injected exactly the garbage the restart existed to remove -- which
% is why a 25-radius interval behaved far worse than a 50-radius one.

RESET_EVERY = 50;     % see note 2: robustness knob, not load-bearing
FEAS_TOL    = 1e-6;   % see note 1

radius = 4.5:0.1:60;
n = numel(radius);
max_vel_corner_vector = zeros(1,n);
x_matrix = [];
guess = [];           % first radius has no predecessor -> analytic seed
nRestart = 0; nSwapped = 0;

for i = 1:n
    [xA,vA,gA,cA] = max_vel_cornering(radius(i),car.max_vel,car,guess);
    okA = (xA(1) == 1 || xA(1) == 2);

    % try a fresh solve when the chain is due a check, has broken physics, or
    % has produced something that is not properly feasible
    dueRestart = mod(i-1,RESET_EVERY) == 0;
    brokeMono  = i > 1 && vA < max_vel_corner_vector(i-1);
    if ~isempty(guess) && (dueRestart || brokeMono || ~okA || cA > FEAS_TOL)
        nRestart = nRestart + 1;
        [xB,vB,gB,cB] = max_vel_cornering(radius(i),car.max_vel,car);
        okB = (xB(1) == 1 || xB(1) == 2);
        if pick(okA,cA,vA,okB,cB,vB,FEAS_TOL)
            xA = xB; vA = vB; gA = gB; cA = cB; %#ok<NASGU>
            nSwapped = nSwapped + 1;
        end
    end

    max_vel_corner_vector(i) = vA;
    x_matrix = [x_matrix; xA]; %#ok<AGROW>
    guess = gA;
end

if nSwapped > 0
    fprintf(['  vel_cornering_sweep: %d restarts attempted, %d preferred over ' ...
             'the chained solve\n'],nRestart,nSwapped);
end

x_table_corner_vel = generate_table(x_matrix);

% visual check if necessary
% plot(radius,max_vel_corner_vector);

end

function takeB = pick(okA,cA,vA,okB,cB,vB,feasTol)
% Feasibility gate first, speed only as a tie-break between two solves that
% both cleared it. Ranking on speed alone is what let a drifted, constraint-
% violating solve win; ranking on residual alone would happily accept a
% trivially feasible but badly suboptimal point.
if ~okB, takeB = false; return, end          % never commit to a failed solve
if ~okA, takeB = true;  return, end
goodA = cA <= feasTol;  goodB = cB <= feasTol;
if goodA && goodB
    takeB = vB > vA;                          % both legal -> faster wins
elseif goodB
    takeB = true;                             % only B is legal
elseif goodA
    takeB = false;
else
    takeB = cB < cA;                          % neither legal -> less illegal
end
end
