function car = makeGG(paramArr,car)
ggpoints = [];
ss_info = [];
accel_info = [];
decel_info = [];
longAccelLookup = [];
longDecelLookup = [];
arr = reshape(paramArr,[numel(paramArr) 1]);
for i = 1:numel(arr)
    pSet = arr(i);
    p1 = []; p2 = []; p3 = []; p4 = []; x0 = []; x1 = []; x2 = [];

    if converged(pSet.maxLatFlag)
        x0 = pSet.maxLatx;
    end

    %x: long %y: lat %z: velo
    if converged(pSet.maxLatFlag) && converged(pSet.maxLongFlag)
        p1 = [pSet.maxLongLongAccel pSet.maxLongLatAccel pSet.longVel];
        p2 = [pSet.maxLongLongAccel -pSet.maxLongLatAccel pSet.longVel];
        x1 = pSet.maxLongAccelx;
    end

    if converged(pSet.maxLatFlag) && converged(pSet.maxBrakeFlag)
        p3 = [pSet.maxBrakeLongDecel pSet.maxBrakeLatAccel pSet.longVel];
        p4 = [pSet.maxBrakeLongDecel -pSet.maxBrakeLatAccel pSet.longVel];
        x2 = pSet.maxBrakeDecelx;
    end

    ggpoints = [ggpoints; p1; p2; p3; p4];
    ss_info = [ss_info; x0];
    accel_info = [accel_info; x1];
    decel_info = [decel_info; x2];
    longAccelLookup = [longAccelLookup; p1];
    longDecelLookup = [longDecelLookup; p3];
end

car.ggPoints = ggpoints;
car.ss_info = ss_info;
car.accel_info = accel_info;
car.decel_info = decel_info;
car.longAccelLookup = longAccelLookup;
car.longDecelLookup = longDecelLookup;

% Say what was thrown away. A dropped velocity removes its whole row -- 20
% accel and 20 brake points -- from the table the lap sim interpolates, and
% until this printed there was no way to tell a 29-velocity g-g from a
% 27-velocity one short of counting rows by hand. That is how a systematic
% non-convergence at 9 m/s survived a full aero sweep unnoticed.
allV  = unique(round([arr.longVel],6));
gotV  = [];
if ~isempty(ss_info), gotV = unique(round(ss_info(:,6),6)); end
lostV = setdiff(allV,gotV);
if ~isempty(lostV)
    % Report the exitflag alongside the velocity. This warning is intermittent
    % by nature -- it depends on the car -- so the flag is the difference
    % between a reproducible diagnosis and another round of guessing. 0 means
    % max_lat_accel ran out of function evaluations, -2 that it found no
    % feasible point at all; those want different fixes. max_lat_accel solves
    % from two starts and reports the better, so reaching here means BOTH
    % failed, and the cheapest next move is to add a third start there.
    lostFlag = zeros(size(lostV));
    for k = 1:numel(lostV)
        m = abs(round([arr.longVel],6) - lostV(k)) < 1e-9;
        f = [arr(m).maxLatFlag];
        lostFlag(k) = f(1);
    end
    warning('makeGG:droppedVelocities', ...
        ['g-g solved %d of %d velocities -- no lateral solution at %s m/s ' ...
         '(fmincon exitflag %s), so those rows are missing from the accel ' ...
         'and brake lookups too.'], ...
        numel(gotV),numel(allV),strjoin(compose('%.2f',lostV),', '), ...
        strjoin(compose('%d',lostFlag),', '));
end
end

function ok = converged(exitflag)
% fmincon exitflag 1 = first-order optimality met, 2 = step below tolerance.
% BOTH are successful terminations. This used to test == 1 only, which threw
% away every exitflag-2 point -- and because the accel and brake rows are
% gated on the lateral flag too, one discarded max_lat_accel took its whole
% velocity row of 20 accel and 20 brake points with it.
%
% Measured: the exitflag-2 points are not marginal. At v = 9 m/s the residual
% is max|ceq| = 1e-10, an order of magnitude BETTER than the worst
% exitflag-1 point on the same grid (1.4e-9).
%
% This also makes the g-g agree with the rest of the toolchain --
% aeroEnvelope.filterSolved, ggMetrics and rampSweep all accept 1 or 2.
ok = (exitflag == 1) | (exitflag == 2);
end
