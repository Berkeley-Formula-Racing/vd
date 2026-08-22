function R = rampSweep(car,opts)
% Constant-speed lateral-acceleration ramps across a range of vCar, logging
% the balance metrics that decide how the car actually behaves.
%
% At each speed the car is solved at a series of prescribed steady-state
% lateral accelerations, from near zero up to its own limit. This is a
% ramp-steer test run backwards -- prescribe Ay and solve for the steer angle
% that produces it, rather than sweeping steer and recording whatever Ay comes
% out. Two reasons: the ramp lands on exactly the accelerations asked for, so
% speeds are directly comparable, and near the limit the delta -> Ay map goes
% flat, where a steer sweep wastes most of its points.
%
% Above the limit the map inverts and there are two steer angles for one Ay
% (the stable low-steer branch and the unstable high-steer one). The solve
% MINIMISES steer angle, which selects the stable branch -- the one a driver
% is actually on during a ramp.
%
% inputs  car:  a Car. Does NOT need to have been through gg2/makeGG; this is
%               a standalone sweep that solves its own points.
%         opts: optional struct
%           .speeds    vCar values, m/s          (default 5:2.5:30)
%           .nRamp     points per ramp           (default 12)
%           .ayMaxFrac ramp to this fraction of the solved limit (default 0.97)
%           .ayMinFrac start the ramp here       (default 0.04)
%           .linearFrac  portion of the ramp treated as the linear range when
%                        fitting the understeer gradient (default 0.40)
%           .nBisect   bisection steps used to pin the top of a truncated
%                      ramp (default 6, 0 disables). Without this the
%                      sustainable limit is only known to the ramp's own
%                      resolution, and the _limit metrics land at a different
%                      fraction of the limit at every speed.
%           .ceqTol    largest steady-state residual accepted (default 1e-2,
%                      matching fmincon's ConstraintTolerance). Points above
%                      it are discarded rather than logged, because fmincon
%                      reports exitflag 2 at infeasible points too.
%           .mode      'balanced' (default) | 'coast'
%                        balanced - longitudinal acceleration held at zero by
%                                   letting the powertrain trim out drag. This
%                                   is the textbook constant-speed ramp and
%                                   the right basis for balance numbers.
%                        coast    - throttle and slip ratios pinned at zero,
%                                   exactly matching the steady-state branch
%                                   of the g-g. The car is then decelerating
%                                   at drag/M, which at 30 m/s is not small
%                                   and loads the rear axle.
%           .verbose   progress printing         (default true)
%
% output  R: struct with
%           .points   table, one row per solved ramp point: every field of
%                     Car.metrics plus the ramp-derived columns below
%           .perSpeed table, one row per speed: the scalars you plot over vCar
%           .settings what was actually run
%
% A note on what moves and what does not, with static aero:
%   LLTD  (front share of lateral load transfer) is pure geometry and roll
%         stiffness split -- both load transfers scale linearly with Ay, so
%         their ratio is constant in BOTH Ay and speed. A flat mech_balance
%         line is correct, not a bug.
%   CoP   (front downforce fraction) is a constant too, by definition of the
%         static aero model.
% So the balance shift with speed comes from neither of those directly. It
% comes from downforce growing the axle loads while the load transfer at a
% given Ay does not, and from drag transferring load rearward. That is what
% LLT_norm_front / LLT_norm_rear and front_Fz_frac track.

if nargin < 2, opts = struct(); end
if numel(opts) ~= 1
    error('rampSweep:optsNotScalar','opts must be a scalar struct, got %d-by-%d', ...
        size(opts,1),size(opts,2));
end

speeds     = getOr(opts,'speeds',5:2.5:30);
nRamp      = getOr(opts,'nRamp',12);
ayMaxFrac  = getOr(opts,'ayMaxFrac',0.97);
ayMinFrac  = getOr(opts,'ayMinFrac',0.04);
linearFrac = getOr(opts,'linearFrac',0.40);
mode       = getOr(opts,'mode','balanced');
ceqTol     = getOr(opts,'ceqTol',1e-2);
nBisect    = getOr(opts,'nBisect',6);
verbose    = getOr(opts,'verbose',true);

if ~any(strcmpi(mode,{'balanced','coast'}))
    error('rampSweep:badMode','mode must be ''balanced'' or ''coast''');
end
holdSpeed = strcmpi(mode,'balanced');

speeds = speeds(:).';
g = car.g;

rows = {};
perSpeed = struct([]);

for iv = 1:numel(speeds)
    v = speeds(iv);

    % --- limit at this speed, which sets the top of the ramp ---
    lastwarn('');
    ayLim = NaN; xLim = [];
    try
        [~,ayLim,~,xLim] = max_lat_accel(v,car);
    catch ME
        if verbose, fprintf('  v=%5.1f  limit solve failed (%s)\n',v,ME.message); end
    end
    if ~isfinite(ayLim) || ayLim <= 0
        if verbose, fprintf('  v=%5.1f  no lateral limit found, skipping\n',v); end
        continue
    end

    ayTargets = linspace(ayMinFrac*ayLim, ayMaxFrac*ayLim, nRamp);

    % --- walk the ramp upward, warm starting each point from the last ---
    x0 = seedState(v,ayTargets(1),xLim);
    ramp = struct([]);
    misses = 0;
    for k = 1:nRamp
        ay = ayTargets(k);

        % Continuation first, then two fallback seeds. Without the retries a
        % merely awkward warm start is indistinguishable from a point the car
        % genuinely cannot reach, and in 'balanced' mode that distinction is
        % the whole story at high speed.
        seeds = {x0, seedState(v,ay,xLim), [2,0,v,0.05,ay/v,0,0,0,0]};
        got = false;
        for si = 1:numel(seeds)
            xs = seeds{si};
            xs(3) = v;
            xs(5) = ay/v;         % the yaw rate the constraint will enforce
            [x,exitflag,ceqMax] = solveRampPoint(car,xs,v,ay,holdSpeed);
            % Exitflag alone is not enough. fmincon's exitflag 2 means the
            % step got small, which it also does when it is stuck at an
            % INFEASIBLE point -- and that is exactly what happens where a
            % constant-speed ramp runs out of power. Accepting those silently
            % booked residuals of 4+ N-m as converged solutions.
            if (exitflag == 1 || exitflag == 2) && ceqMax <= ceqTol
                got = true; break
            end
        end
        if ~got
            misses = misses + 1;
            % three in a row means the ramp is over, not that one point was
            % awkward; keep going and every remaining solve is wasted
            if misses >= 3, break, end
            continue
        end
        misses = 0;
        x0 = x;                   % continuation

        m = car.metrics(x);
        m = addRampFields(m,car,x,ay,exitflag,ceqMax,v,g);
        if isempty(ramp), ramp = m; else, ramp(end+1) = m; end %#ok<AGROW>
    end

    if isempty(ramp)
        if verbose, fprintf('  v=%5.1f  ramp did not converge anywhere\n',v); end
        continue
    end

    % If the ramp stopped short, its top point is wherever the coarse target
    % grid happened to land, which quantises the sustainable limit to the ramp
    % resolution -- ramp_complete comes out as 0.97 / 0.885 / 0.80 and nothing
    % between. Every metric tagged _limit is then read at a different fraction
    % of the limit at each speed, so they cannot be compared down a column.
    % Bisect for the actual edge instead.
    ayOK = max([ramp.gLat])*g;
    nxt  = ayTargets(ayTargets > ayOK + 1e-9);
    if nBisect > 0 && ~isempty(nxt)
        ayBad = nxt(1);
        best  = [];
        for b = 1:nBisect
            aym = 0.5*(ayOK+ayBad);
            hit = false;
            for xs = {x0, seedState(v,aym,xLim)}
                z = xs{1}; z(3) = v; z(5) = aym/v;
                [x,exitflag,ceqMax] = solveRampPoint(car,z,v,aym,holdSpeed);
                if (exitflag == 1 || exitflag == 2) && ceqMax <= ceqTol
                    hit = true; break
                end
            end
            if hit
                ayOK = aym;
                best = addRampFields(car.metrics(x),car,x,aym,exitflag,ceqMax,v,g);
                x0 = x;
            else
                ayBad = aym;
            end
        end
        if ~isempty(best), ramp(end+1) = best; end %#ok<AGROW>
    end

    rows{end+1} = struct2table(ramp); %#ok<AGROW>
    s = summariseRamp(ramp,v,ayLim,linearFrac,car,g);
    if isempty(perSpeed), perSpeed = s; else, perSpeed(end+1) = s; end %#ok<AGROW>

    if verbose
        flag = '';
        if s.ramp_complete < 0.98
            why = 'grip';
            if s.power_limited, why = 'POWER'; end
            flag = sprintf('  <-- stops at %.0f%% of the free limit (%s)', ...
                100*s.ramp_complete, why);
        end
        fprintf(['  v=%5.1f  %2d/%2d pts   gLat %5.3f/%5.3f   K_lin %+7.3f deg/g ' ...
                 '(R2 %.3f)   K_top %+7.3f   LLTD %.3f   frontFz %.3f%s\n'], ...
            v, numel(ramp), nRamp, s.gLat_top, s.gLat_max, s.K_linear, s.K_r2, ...
            s.K_at_limit, s.mech_balance, s.front_Fz_frac_limit, flag);
    end
end

if isempty(rows)
    error('rampSweep:noSolution', ...
        ['not one ramp point converged across %d speeds. Check that the car ' ...
         'solves at all -- max_lat_accel on its own is the place to start.'], ...
        numel(speeds));
end

R.points   = vertcat(rows{:});
R.perSpeed = struct2table(perSpeed);
R.settings = struct('speeds',speeds,'nRamp',nRamp,'ayMaxFrac',ayMaxFrac, ...
    'ayMinFrac',ayMinFrac,'linearFrac',linearFrac,'mode',lower(mode));
end

%% ------------------------------------------------------------------------
function [x,exitflag,ceqMax] = solveRampPoint(car,x0,v,ayTarget,holdSpeed)
% One steady-state point: the smallest steer angle that produces ayTarget.

% bounds mirror max_lat_accel, except that 'balanced' mode has to be allowed
% throttle and rear slip ratio to trim out drag
if holdSpeed
    throttle_bounds = [0,1];
    kappa_r_bounds  = [0,0.25];
else
    throttle_bounds = [0,0];
    kappa_r_bounds  = [0,0];
end

lb = [ 0, throttle_bounds(1), v, -3, 0, 0, 0, kappa_r_bounds(1), kappa_r_bounds(1)];
ub = [25, throttle_bounds(2), v,  3, 2, 0, 0, kappa_r_bounds(2), kappa_r_bounds(2)];

% long_vel pinned, front slip ratios zero (undriven, unbraked), rear slip
% ratios symmetric. The symmetry goes in Aeq rather than being forced inside
% the objective -- objective and constraints must evaluate the same state.
Aeq = [0 0 1 0 0 0 0 0 0
       0 0 0 0 0 1 0 0 0
       0 0 0 0 0 0 1 0 0
       0 0 0 0 0 0 0 1 -1];
beq = [v 0 0 0];

opts = setOptimoptions(1500);
f = @(P) P(1);                       % minimum steer -> the stable branch
nl = @(P) rampConstraint(car,P,ayTarget,holdSpeed);

[x,~,exitflag] = fmincon(f,x0,[],[],Aeq,beq,lb,ub,nl,opts);

[~,ceq] = rampConstraint(car,x,ayTarget,holdSpeed);
ceqMax = max(abs(ceq));
end

function [c,ceq] = rampConstraint(car,P,ayTarget,holdSpeed)
% Steady state at a prescribed lateral acceleration. Built here rather than
% added to Car.m so this sweep stays self-contained: constraint4 does the Ay
% target but cannot also hold longitudinal acceleration at zero, which is
% what makes a constant-SPEED ramp a constant-speed ramp.
[engine_rpm,beta,lat_accel,long_accel,yaw_accel,wheel_accel,~,~,Fzvirtual] = ...
    car.equations(P);
c   = [engine_rpm-13000, abs(beta)-20, -Fzvirtual(1:4)];
ceq = [P(3)*P(5)-ayTarget, lat_accel, yaw_accel, wheel_accel(1:4)];
if holdSpeed
    ceq = [ceq, long_accel];
end
end

function x0 = seedState(v,ay,xLim)
% start from the limit solution if there is one, scaled back toward the
% bottom of the ramp; otherwise a bland guess
if ~isempty(xLim) && numel(xLim) == 9 && all(isfinite(xLim))
    x0 = xLim(:).';
    frac = 0.15;
    x0(1) = max(x0(1)*frac,0.5);     % steer
    x0(4) = x0(4)*frac;              % lat vel
else
    x0 = [2, 0, v, 0.05, ay/v, 0, 0, 0, 0];
end
x0(3) = v;
x0(5) = ay/v;
end

%% ------------------------------------------------------------------------
function m = addRampFields(m,car,x,ayTarget,exitflag,ceqMax,v,g)
% ramp-specific columns layered on top of Car.metrics

m.ay_target = ayTarget/g;            % g, what was asked for
m.exitflag  = exitflag;
m.max_ceq   = ceqMax;

% --- steering ---
% the two front wheels differ through Ackermann; the bicycle-model comparison
% against the Ackermann angle wants their average
[d1,d2] = frontSteer(car,x(1));
m.steer_avg   = (d1+d2)/2;
m.steer_inner = d1;
m.steer_outer = d2;
% radius comes from the state, not from a nominal: R = v/r
m.radius = v/max(abs(x(5)),eps);
m.steer_ackermann = (180/pi)*car.W_b/m.radius;
% the classic understeer measure: steer held above what pure geometry needs
m.steer_excess = m.steer_avg - m.steer_ackermann;

% --- axle slip angles ---
% magnitudes, so the result does not depend on the sign convention: the front
% carrying MORE slip than the rear is understeer either way
af = (abs(m.alpha_1)+abs(m.alpha_2))/2;
ar = (abs(m.alpha_3)+abs(m.alpha_4))/2;
m.alpha_f = af;
m.alpha_r = ar;
m.alpha_balance = af - ar;           % >0 understeer, <0 oversteer

% --- how much of each axle's grip is spent ---
FzF = max(m.Fz_front_axle,eps);
FzR = max(m.Fz_rear_axle, eps);
m.mu_f_used = abs(m.Fy_front)/FzF;
m.mu_r_used = abs(m.Fy_rear)/FzR;
m.grip_balance = m.mu_f_used - m.mu_r_used;   % >0 front closer to saturation

% --- load transfer normalised by the axle it acts on ---
% this is the balance term that actually moves with speed: the transfer at a
% given Ay is fixed, but downforce keeps growing the axle load underneath it
m.LLT_norm_front = abs(m.LLT_front)/FzF;
m.LLT_norm_rear  = abs(m.LLT_rear)/FzR;
m.LLT_norm_balance = m.LLT_norm_front - m.LLT_norm_rear;

% --- gains ---
sa = m.steer_avg;
if abs(sa) > 1e-6
    m.yaw_rate_gain = x(5)/sa;       % (rad/s) per deg of steer
    m.Ay_gain       = m.gLat/sa;     % g per deg of steer
else
    m.yaw_rate_gain = NaN;
    m.Ay_gain       = NaN;
end

% --- aero at this speed ---
m.downforce_frac = m.downforce/max(m.downforce + car.M*g,eps);
m.drag_decel     = car.aero.drag(v)/(car.M*g);   % g, what coast mode bleeds

% roll angle uses the SAME hardcoded 0.68 deg/g gradient that
% Camber_Evaluation applies. It is not derived from carConfig spring rates,
% so treat it as a placeholder, not a result.
m.roll_angle = abs(m.gLat)*0.68;
end

function [d1,d2] = frontSteer(car,steer_angle)
% reproduces the Ackermann split inside Car.equations
a_poly = [1.000170221169974e-05,-6.101610470854148e-05,0.009416280591935,0.999321611021116];
if steer_angle >= 0
    d1 = steer_angle;
    d2 = abs(car.ackermann)*polyval(a_poly,steer_angle)*steer_angle;
else
    d2 = steer_angle;
    d1 = -(abs(car.ackermann)*polyval(a_poly,-steer_angle)*(-steer_angle));
end
if car.ackermann < 0
    [d1,d2] = deal(d2,d1);
end
end

%% ------------------------------------------------------------------------
function s = summariseRamp(ramp,v,ayLim,linearFrac,~,g)
% collapse one ramp into the scalars that get plotted over vCar

gLat  = [ramp.gLat].';
exc   = [ramp.steer_excess].';
yawr  = [ramp.yaw_rate].';
steer = [ramp.steer_avg].';

s.vCar     = v;
s.nPoints  = numel(ramp);
% gLat_max is the FREE lateral limit, solved with no constraint on
% longitudinal acceleration. gLat_top is how far up the ramp actually got.
% In 'balanced' mode those two separate at high speed and the gap is real
% physics, not a solver failure: holding the speed constant costs the rear
% tires drag/M of longitudinal force, and eventually there is not enough of
% the friction ellipse left to also make the lateral force. Every field
% below tagged _limit is taken at gLat_top, so when ramp_complete < 1 they
% describe the highest SUSTAINABLE point, not the car's limit.
%
% Treat a truncated gLat_top as a LOWER bound. Right at the constant-speed
% feasibility boundary the solve itself degrades, so bisection can fail to
% place points that are in fact reachable. K_at_limit inherits that: it is
% trustworthy where ramp_complete is near 1 and noisy where it is not,
% because it is then fitted at a different fraction of the limit per speed.
% K_linear does not have this problem -- it is fitted well below the boundary.
s.gLat_max = ayLim/g;
s.gLat_top = max(gLat);
s.ramp_complete = s.gLat_top/max(s.gLat_max,eps);

% --- understeer gradient ---
% linear range: the bottom linearFrac of the ramp. Fitting the whole ramp
% would fold the limit rolloff into the gradient and report a number the car
% never has at any single Ay.
lin = gLat <= linearFrac*s.gLat_max;
[s.K_linear,s.K_r2] = fitSlope(gLat(lin),exc(lin));
% terminal gradient, over the top quarter of the ramp's SPAN rather than a
% fixed three points: near the limit the steer excess curves hard, so a
% three-point fit rides the local noise and swings tens of deg/g between
% adjacent speeds. NaN unless there are enough points to mean anything.
span = s.gLat_top - min(gLat);
if span > 0
    top = gLat >= s.gLat_top - 0.25*span;
else
    top = false(size(gLat));
end
if sum(top) >= 3
    s.K_at_limit = fitSlope(gLat(top),exc(top));
else
    s.K_at_limit = NaN;
end
% yaw gain in the same linear range, for reference
s.yaw_gain_linear = fitSlope(steer(lin),yawr(lin));

% --- balance, at the limit and at mid-ramp ---
% mid-ramp is where a driver spends most of a lap; the limit is where the car
% decides the corner. They are frequently not the same story.
[~,iMid] = min(abs(gLat - 0.5*s.gLat_max));
iLim = numel(gLat);

s.mech_balance      = ramp(iLim).LLTD;        % constant in Ay and in speed
s.aero_balance      = ramp(iLim).CoP;         % constant, static aero
s.downforce         = ramp(iLim).downforce;
s.downforce_frac    = ramp(iLim).downforce_frac;
s.drag_decel        = ramp(iLim).drag_decel;

s.front_Fz_frac_mid   = ramp(iMid).front_Fz_frac;
s.front_Fz_frac_limit = ramp(iLim).front_Fz_frac;

s.LLT_norm_balance_mid   = ramp(iMid).LLT_norm_balance;
s.LLT_norm_balance_limit = ramp(iLim).LLT_norm_balance;

s.grip_balance_mid    = ramp(iMid).grip_balance;
s.grip_balance_limit  = ramp(iLim).grip_balance;
s.alpha_balance_mid   = ramp(iMid).alpha_balance;
s.alpha_balance_limit = ramp(iLim).alpha_balance;

s.beta_limit    = ramp(iLim).beta;
s.steer_limit   = ramp(iLim).steer_avg;
s.min_Fz_limit  = ramp(iLim).min_Fz;      % negative => a wheel has lifted
s.roll_limit    = ramp(iLim).roll_angle;

% --- why the ramp stopped ---
% throttle pinned at the stop means the car ran out of power to hold this
% speed, not out of tyre. At high speed that is what ends the ramp, and it
% is a different engineering problem from running out of grip.
s.throttle_top   = ramp(iLim).throttle;
s.power_limited  = s.throttle_top > 0.99 && s.ramp_complete < 0.98;

% --- solve health, so a bad number is visible as a bad number ---
% exitflag is NOT ordered by quality (1 best, 2 acceptable, 0 and -2 bad), so
% min() would call a converged ramp the worst one. Count instead.
s.max_ceq  = max([ramp.max_ceq]);
s.n_exit1  = sum([ramp.exitflag] == 1);
s.n_exit2  = sum([ramp.exitflag] == 2);
end

function [k,r2] = fitSlope(x,y)
% least squares slope with R^2, base MATLAB only
k = NaN; r2 = NaN;
x = x(:); y = y(:);
ok = isfinite(x) & isfinite(y);
x = x(ok); y = y(ok);
if numel(x) < 2 || max(x)-min(x) <= 0, return, end
p  = polyfit(x,y,1);
k  = p(1);
yh = polyval(p,x);
ss = sum((y-mean(y)).^2);
if ss > 0, r2 = 1 - sum((y-yh).^2)/ss; else, r2 = NaN; end
end

function v = getOr(s,f,d)
if isfield(s,f), v = s.(f); else, v = d; end
end
