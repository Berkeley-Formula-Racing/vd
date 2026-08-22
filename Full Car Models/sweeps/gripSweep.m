function G = gripSweep(baseCell,opts)
% Sweeps front and rear grip scaling and reports event times plus handling
% balance, so the sim can be correlated against measured track data.
%
% Front and rear separately because the events separate them: skidpad is
% close to pure lateral balance, accel is almost entirely rear longitudinal,
% autocross mixes both. One event cannot identify two parameters.
%
% TIMES ALONE ARE NOT ENOUGH. Many (front,rear) pairs give nearly the same
% lap time with very different handling -- losing front grip and losing rear
% grip both cost time. That is why this also targets a BALANCE at a chosen
% speed: it is the constraint that separates pairs the stopwatch cannot.
%
% inputs  baseCell: carConfig() output, a SINGLE case
%         opts:     optional struct
%           .front        front multipliers (default 0.90:0.05:1.10)
%           .rear         rear multipliers  (default same as .front)
%           .events       events to run (default skidpad/accel/autocross;
%                         endurance costs more than the other three together)
%           .numWorkers   workers (default 8)
%           .loadFromPrev reuse a previous solve (default false)
%           .cacheName    .mat to read/write (default "grip_sweep.mat")
%           .verbose      progress printing (default true)
%           .refineGrid   mesh points per axis for locating the optimum
%                         between grid nodes (default 201)
%           .eventParams  REQUIRED. Second output of carConfig -- the driver
%                         factors, course geometry and lap count.
%
%         --- matching targets, all optional ---
%           .target       measured event times, e.g.
%                           struct('skidpad',4.95,'autocross',48.2)
%           .timeTol      time miss treated as "one unit" of error, as a
%                         fraction (default 0.01, i.e. 1%)
%           .balanceSpeed vCar (m/s) at which balance is evaluated. Setting
%                         this turns the balance target on.
%           .balanceTarget wanted balance value at that speed (default 0,
%                         i.e. neutral -- the car "rebalances" there)
%           .balanceMetric which number to match (default 'grip_balance_limit')
%                           grip_balance_limit  mu used front - rear at the
%                                               top of the ramp; >0 means the
%                                               front saturates first
%                           K_linear            understeer gradient, deg/g
%                           alpha_balance_limit |alpha_f| - |alpha_r|, deg
%           .balanceTol   balance miss treated as one unit (default 0.02)
%           .balanceMode  'coast' (default) | 'balanced', passed to rampSweep
%
% output  G: struct with
%           .times    table: case, gripF, gripR, one column per event,
%                     balance if requested, error columns and .score
%           .best     the lowest-score GRID row, when anything was targeted
%           .bestInterp the optimum INTERPOLATED between grid nodes: the
%                     better fit, but predicted rather than solved. Empty when
%                     the grid is one point or any case failed to solve.
%           .cars     solved carCell, so nothing needs re-solving
%           .settings what was run
%
% Scoring: every term is (error / its tolerance), so seconds and a
% dimensionless balance can share one number. score is the RMS of those, and
% score = 1 means "off by about one tolerance on average". Raw seconds would
% let a 1500 s endurance drown out a 5 s skidpad, and the short events are
% exactly the ones that separate front from rear.
%
% A caution on reading the result. Matching times does not make the tire model
% right -- it makes two numbers absorb every other error in the sim. If the
% best fit needs an axle far from 1, that is evidence about something else
% (aero, load transfer, the driver-conservatism factors baked into
% Events2.Autocross and .Endurance), not a measurement of your grip.

if nargin < 2, opts = struct(); end
if numel(opts) ~= 1
    error('gripSweep:optsNotScalar','opts must be a scalar struct, got %d-by-%d', ...
        size(opts,1),size(opts,2));
end
if size(baseCell,1) ~= 1
    error('gripSweep:notSingle', ...
        ['baseCell must hold exactly one case but holds %d. Set every swept ' ...
         'vector in carConfig back to a scalar -- this function builds the ' ...
         'grid.'], size(baseCell,1));
end

front      = getOr(opts,'front',0.90:0.05:1.10);
rear       = getOr(opts,'rear',front);
events     = cellstr(getOr(opts,'events',{'skidpad','accel','autocross'}));
numWorkers = getOr(opts,'numWorkers',8);
cacheName  = getOr(opts,'cacheName',"grip_sweep.mat");
loadPrev   = getOr(opts,'loadFromPrev',false);
verbose    = getOr(opts,'verbose',true);
% fine mesh used to locate the interpolated optimum between grid nodes. 201
% points per axis is ~40k interpolant evaluations -- microseconds -- and
% resolves the optimum to well under a thousandth of a grip step.
refineGrid = getOr(opts,'refineGrid',201);

% Required, and deliberately not defaulted here: the driver factors and the
% course geometry belong in carConfig, and a default in this file would be a
% second place for them to disagree from.
eventParams = getOr(opts,'eventParams',[]);
if isempty(eventParams)
    error('gripSweep:noEventParams', ...
        ['opts.eventParams is required. Take it from carConfig:\n' ...
         '    [baseCell,eventParams] = carConfig();\n' ...
         '    opts.eventParams = eventParams;']);
end

target        = getOr(opts,'target',struct());
timeTol       = getOr(opts,'timeTol',0.01);
balanceSpeed  = getOr(opts,'balanceSpeed',[]);
balanceTarget = getOr(opts,'balanceTarget',0);
balanceMetric = getOr(opts,'balanceMetric','grip_balance_limit');
balanceTol    = getOr(opts,'balanceTol',0.02);
balanceMode   = getOr(opts,'balanceMode','coast');

front = front(:).'; rear = rear(:).';
if any(front <= 0) || any(rear <= 0)
    error('gripSweep:badLevels','grip multipliers must be positive');
end
if ~any(strcmp(balanceMetric,{'grip_balance_limit','K_linear','alpha_balance_limit'}))
    error('gripSweep:badMetric',['balanceMetric must be grip_balance_limit, ' ...
        'K_linear or alpha_balance_limit']);
end
doBalance = ~isempty(balanceSpeed);
if doBalance && ~isscalar(balanceSpeed)
    error('gripSweep:badSpeed','balanceSpeed must be a single velocity');
end

%% case grid -- full cross product
[F,R] = ndgrid(front,rear);
gripF = F(:); gripR = R(:);
n = numel(gripF);

carCell = cell(n,2);
for i = 1:n
    c = baseCell{1,1};  c.gripScaleF = gripF(i);  c.gripScaleR = gripR(i);
    % the accel-configuration car gets the same scaling: same tyres on the
    % same surface, only a different wing
    a = baseCell{1,2};  a.gripScaleF = gripF(i);  a.gripScaleR = gripR(i);
    carCell{i,1} = c;   carCell{i,2} = a;
end

if verbose
    fprintf('grip sweep: %d cases (%d front x %d rear), events: %s\n', ...
        n,numel(front),numel(rear),strjoin(events,', '));
    if doBalance
        fprintf('  balance target: %s = %+.4g at %.1f m/s (tol %.3g, %s mode)\n', ...
            balanceMetric,balanceTarget,balanceSpeed,balanceTol,balanceMode);
    end
end

%% solve, or reload
haveGG = false; haveTimes = false; haveBal = false;
times = []; balance = [];
if loadPrev && isfile(cacheName)
    prev = load(cacheName);
    if all(isfield(prev,{'carCell','gripF','gripR'})) && numel(prev.gripF) == n && ...
            max(abs(prev.gripF(:)-gripF)) < 1e-12 && ...
            max(abs(prev.gripR(:)-gripR)) < 1e-12 && ...
            ~isempty(prev.carCell{1,1}.ss_info)
        carCell = prev.carCell; haveGG = true;
        if isfield(prev,'times') && ~isempty(prev.times)
            times = prev.times; haveTimes = all(ismember(events, ...
                times.Properties.VariableNames)) && ...
                ~any(all(ismissing(times{:,events}),1));
        end
        if isfield(prev,'balance') && isfield(prev,'balanceKey') && doBalance && ...
                isequal(prev.balanceKey,{balanceMetric,balanceSpeed,balanceMode})
            balance = prev.balance; haveBal = true;
        end
        if verbose
            fprintf('  loaded %s: g-g=%d times=%d balance=%d\n', ...
                cacheName,haveGG,haveTimes,haveBal);
        end
    else
        warning('gripSweep:staleCache','%s does not match this grid -- re-solving.',cacheName);
    end
end

t0 = tic;
job = simLog.start('grip sweep');

%% g-g
if ~haveGG
    % Parallelise across CASES with a serial gg2 inside. gg2's own parfor
    % runs over ~29 velocities and one of them (around 9 m/s) takes an order
    % of magnitude longer than the rest, so most workers sit idle waiting on
    % the straggler. Splitting by case instead keeps every worker busy.
    gg = cell(n,1);
    parfor (i = 1:n, numWorkers)
        gg{i} = makeGG(gg2(carCell{i,1},0),carCell{i,1});
    end
    for i = 1:n, carCell{i,1} = gg{i}; end
    if verbose, fprintf('  g-g done (%.0f s)\n',toc(t0)); end
end

%% events
if ~haveTimes
    if verbose, fprintf('  events ...\n'); end
    times = aeroLapTimes(carCell,numWorkers,events,eventParams);
    if verbose, fprintf('  events done (%.0f s)\n',toc(t0)); end
end

%% balance at the chosen speed
if doBalance && ~haveBal
    if verbose, fprintf('  balance ramps at %.1f m/s ...\n',balanceSpeed); end
    balance = nan(n,1);
    ro = struct('speeds',balanceSpeed,'nRamp',10,'mode',balanceMode, ...
        'verbose',false,'nBisect',4);
    cars1 = carCell(:,1);
    parfor (i = 1:n, numWorkers)
        try
            Ri = rampSweep(cars1{i},ro);
            balance(i) = Ri.perSpeed.(balanceMetric)(1);
        catch ME
            % one unsolvable ramp must not take the sweep down; it shows up
            % as NaN and is excluded from the score
            warning('gripSweep:rampFailed','case %d ramp failed: %s',i,ME.message);
        end
    end
    if verbose, fprintf('  balance done (%.0f s)\n',toc(t0)); end
end

if ~haveGG || ~haveTimes || (doBalance && ~haveBal)
    % saved so a cached balance can be rejected when the metric, speed or
    % mode changes -- the times stay valid but the balance column does not
    balanceKey = {balanceMetric,balanceSpeed,balanceMode};
    save(cacheName,'carCell','gripF','gripR','times','balance','balanceKey','-v7.3');
end

%% assemble
T = table((1:n)',gripF,gripR,'VariableNames',{'case','gripF','gripR'});
for k = 1:numel(events)
    T.(events{k}) = times.(events{k});
end
if doBalance
    T.balance = balance(:);
end

%% score against the targets
tf = fieldnames(target);
tf = tf(ismember(tf,events));
if ~isempty(fieldnames(target)) && isempty(tf)
    warning('gripSweep:targetUnmatched', ...
        'target names none of the events that were run (%s)',strjoin(events,', '));
end

terms = zeros(n,0);
for k = 1:numel(tf)
    e = tf{k};
    err = T.(e) - target.(e);
    T.(['err_' e])    = err;
    T.(['pcterr_' e]) = 100*err/target.(e);
    % relative miss, expressed in units of timeTol. Grown rather than
    % preallocated: at most four events plus balance, so the column count is
    % single digits and preallocating would only obscure which terms exist.
    terms(:,end+1) = (err/target.(e))/timeTol; %#ok<AGROW>
end
if doBalance
    berr = T.balance - balanceTarget;
    T.err_balance = berr;
    % absolute, not relative: the target is routinely 0 (neutral), and a
    % fractional error against 0 is undefined
    terms(:,end+1) = berr/balanceTol;
end

if ~isempty(terms)
    % NaN in any term (a failed event or ramp) makes that case unscoreable
    T.score = sqrt(mean(terms.^2,2));
    T.score(any(~isfinite(terms),2)) = NaN;
    [~,b] = min(T.score);
    G.best = T(b,:);

    % The grid almost never lands on the optimum, and solving a finer grid
    % costs a g-g per node. Interpolate instead -- but interpolate the EVENT
    % TIMES and re-score on a fine mesh, NOT the score surface itself. The
    % times are smooth and close to linear in grip; score is an RMS of their
    % errors, so it has a kink wherever a term crosses its target. Scoring
    % interpolated times puts that kink exactly where it belongs, which is
    % usually within one grid cell of the answer -- interpolating score
    % directly would round the kink off and pull the minimum away from it.
    G.bestInterp = refineBest(T,tf,target,timeTol,doBalance,balanceTarget, ...
        balanceTol,front,rear,refineGrid);

    if verbose && isfinite(T.score(b))
        fprintf('\n  best GRID node: gripF %.3f, gripR %.3f   score %.2f\n', ...
            T.gripF(b),T.gripR(b),T.score(b));
        for k = 1:numel(tf)
            fprintf('        %-10s sim %8.3f  target %8.3f  (%+.2f%%)\n', ...
                tf{k},T.(tf{k})(b),target.(tf{k}),T.(['pcterr_' tf{k}])(b));
        end
        if doBalance
            fprintf('        %-10s sim %8.4f  target %8.4f  (%+.4f)\n', ...
                'balance',T.balance(b),balanceTarget,T.err_balance(b));
        end

        BI = G.bestInterp;
        if ~isempty(fieldnames(BI)) && isfinite(BI.score)
            fprintf('\n  best INTERPOLATED: gripF %.4f, gripR %.4f   score %.2f  (%s)\n', ...
                BI.gripF,BI.gripR,BI.score,BI.method);
            for k = 1:numel(tf)
                fprintf('        %-10s pred %8.3f  target %8.3f  (%+.2f%%)\n', ...
                    tf{k},BI.(tf{k}),target.(tf{k}), ...
                    100*(BI.(tf{k})-target.(tf{k}))/target.(tf{k}));
            end
            if doBalance
                fprintf('        %-10s pred %8.4f  target %8.4f  (%+.4f)\n', ...
                    'balance',BI.balance,balanceTarget,BI.balance-balanceTarget);
            end
            fprintf(['        those are INTERPOLATED, not solved. Confirm with a\n' ...
                     '        1x1 sweep at that pair before trusting them.\n']);
        end

        % edge test on whichever pair is being recommended
        if ~isempty(fieldnames(BI)) && isfinite(BI.score)
            eF = BI.gripF; eR = BI.gripR;
        else
            eF = T.gripF(b); eR = T.gripR(b);
        end
        tolF = 1e-9 + 1e-6*max(abs(front));
        tolR = 1e-9 + 1e-6*max(abs(rear));
        onEdge = eF <= min(front)+tolF || eF >= max(front)-tolF || ...
                 eR <= min(rear)+tolR  || eR >= max(rear)-tolR;
        if onEdge
            fprintf(['  NOTE: the best pair sits on the edge of the grid, so the true\n' ...
                     '        optimum is probably outside it -- widen .front / .rear\n' ...
                     '        and re-run before believing this pair.\n']);
        end
        % close the loop: these are the lines that lock the fit into the car
        fprintf(['\n  to keep this fit, put these in carConfig.m:\n' ...
                 '    tireParams.grip_scaling_front = %.4f;\n' ...
                 '    tireParams.grip_scaling_rear  = %.4f;\n'], eF, eR);
    end
end

G.times    = T;
G.cars     = carCell;
G.target   = target;
G.settings = struct('front',front,'rear',rear,'events',{events}, ...
    'numWorkers',numWorkers,'cacheName',cacheName,'timeTol',timeTol, ...
    'balanceSpeed',balanceSpeed,'balanceTarget',balanceTarget, ...
    'balanceMetric',balanceMetric,'balanceTol',balanceTol,'balanceMode',balanceMode);
G.grid     = struct('F',F,'R',R);

% one row for the whole sweep: the elapsed time reflects whatever was actually
% solved (the cache flags say how much was reused), and the best pair is the
% headline result worth keeping next to the timing
best = struct();
if isfield(G,'best'), best = G.best; end
simLog.finish(job,'events',events,'workers',numWorkers,'nCases',n, ...
    'details',sprintf('loadedGG=%d loadedTimes=%d cache=%s%s',haveGG,haveTimes, ...
        char(cacheName), gripSummary(best)));
end

function s = gripSummary(best)
if isfield(best,'gripF') && isfield(best,'gripR')
    s = sprintf(' best=%.3f/%.3f',best.gripF,best.gripR);
else
    s = '';
end
end

function v = getOr(s,f,d)
if isfield(s,f), v = s.(f); else, v = d; end
end


function BI = refineBest(T,tf,target,timeTol,doBalance,balanceTarget, ...
    balanceTol,front,rear,nMesh)
% Locate the best (gripF,gripR) BETWEEN grid nodes.
%
% Interpolates each targeted quantity -- the event times, and balance if it
% was evaluated -- over the swept grid, then rebuilds the score on a fine mesh
% using exactly the same formula the grid scoring uses. The minimum of that
% mesh is the answer.
%
% Interpolating the quantities rather than the score matters. Event times are
% smooth and nearly linear in grip, so they interpolate honestly. score is an
% RMS of their errors: it has a kink wherever a term crosses its target, and
% the optimum usually sits ON that kink. Interpolating score directly would
% smooth the kink away and slide the minimum off it.
%
% Returns an empty struct when there is nothing to interpolate (a single grid
% point, or a NaN anywhere in the quantities that would poison the mesh).

BI = struct();
nF = numel(front); nR = numel(rear);
if nF*nR < 2, return, end

% every quantity the score depends on
q = tf(:).';
if doBalance, q = [q {'balance'}]; end
if isempty(q), return, end

% A NaN anywhere means a case failed to solve. Interpolating across it would
% invent capability the sweep never demonstrated, so refuse rather than guess.
for k = 1:numel(q)
    if any(~isfinite(T.(q{k})))
        warning('gripSweep:refineNaN', ...
            ['not interpolating the optimum: %s has a non-finite value, so ' ...
             'at least one case did not solve. The grid best still stands.'], q{k});
        return
    end
end

% method by how much support each axis has. Spline needs 4 points; makima is
% well behaved at 3; 2 leaves nothing but a straight line anyway.
if nF > 1 && nR > 1, nMin = min(nF,nR); else, nMin = max(nF,nR); end
if     nMin >= 4, method = 'spline';
elseif nMin == 3, method = 'makima';
else,             method = 'linear';
end

fq = linspace(min(front),max(front),nMesh);
rq = linspace(min(rear), max(rear), nMesh);

% Build one interpolant per quantity and evaluate on the mesh. ndgrid put
% front along dim 1, so reshape matches without a transpose.
vals = struct();
if nF > 1 && nR > 1
    [FQ,RQ] = ndgrid(fq,rq);
    for k = 1:numel(q)
        M = reshape(T.(q{k}),[nF nR]);
        Fi = griddedInterpolant({front,rear},M,method);
        vals.(q{k}) = Fi(FQ,RQ);
    end
elseif nF > 1          % rear is a single value
    FQ = fq(:); RQ = repmat(rear,numel(fq),1);
    for k = 1:numel(q)
        Fi = griddedInterpolant(front,T.(q{k})(:).',method);
        vals.(q{k}) = Fi(FQ);
    end
else                   % front is a single value
    RQ = rq(:); FQ = repmat(front,numel(rq),1);
    for k = 1:numel(q)
        Fi = griddedInterpolant(rear,T.(q{k})(:).',method);
        vals.(q{k}) = Fi(RQ);
    end
end

% same score as the grid: every term is (error / its tolerance)
terms = [];
for k = 1:numel(tf)
    e = tf{k};
    terms = cat(3,terms,((vals.(e)-target.(e))/target.(e))/timeTol);
end
if doBalance
    terms = cat(3,terms,(vals.balance-balanceTarget)/balanceTol);
end
score = sqrt(mean(terms.^2,3));

[bestScore,idx] = min(score(:));
BI.gripF  = FQ(idx);
BI.gripR  = RQ(idx);
BI.score  = bestScore;
BI.method = method;
for k = 1:numel(q)
    BI.(q{k}) = vals.(q{k})(idx);
end
end
