%% run_grip_sweep
% Fits front and rear tire grip so the sim matches the real car.

clear classes
setup_paths

%% grid
opts = struct();
opts.front = 0.595:0.005:0.62;
opts.rear = 0.595:0.005:0.62;

opts.events = {'skidpad','accel','autocross'};

%% real data
MEASURED = struct();
MEASURED.skidpad   = 4.9;
MEASURED.accel     = 4.80;
MEASURED.autocross = 48;

opts.timeTol = 0.01;

%% rebalance

opts.balanceSpeed  = 13.75;         % m/s, where the car swaps understeer -> oversteer
opts.balanceTarget = 0;             % 0 = neutral at that speed
opts.balanceMetric = 'K_linear';    % understeer gradient, deg/g
opts.balanceTol    = 0.01;          % deg/g
opts.balanceMode   = 'coast';

%% run

opts.numWorkers   = 16;
opts.cacheName    = "grip_sweep.mat";
opts.loadFromPrev = false;   % reuse a previous solve of the SAME grid

opts.target = MEASURED;
if isempty(fieldnames(MEASURED)) && isempty(opts.balanceSpeed)
    warning('run_grip_sweep:noTargets', ...
        ['nothing to match: MEASURED is empty and no balanceSpeed is set, so ' ...
         'this will just tabulate the grid without scoring it.']);
end

[baseCell,eventParams] = carConfig();
opts.eventParams = eventParams;

G = gripSweep(baseCell,opts);

%% ---------------- results ----------------
disp(G.times(:,intersect(G.times.Properties.VariableNames, ...
    [{'gripF','gripR'},opts.events,{'balance','score'}],'stable')));

plotGripSweep(G);
