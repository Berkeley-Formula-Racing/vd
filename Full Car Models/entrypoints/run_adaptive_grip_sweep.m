function result = run_adaptive_grip_sweep(targetTimes,options)
%RUN_ADAPTIVE_GRIP_SWEEP Adaptively fit front/rear grip to event times.
%
%   result = run_adaptive_grip_sweep(targetTimes)
%   result = run_adaptive_grip_sweep(targetTimes,options)
%
% TARGETTIMES is a scalar struct with positive times in seconds:
%   struct('autocross',48,'skidpad',4.9,'accel',4.8)
%
% The first pass uses the configured front/rear grip grids. Each following
% pass shrinks both axes around gripSweep's best interpolated candidate. A
% final 1-by-1 sweep solves that candidate instead of reporting interpolation
% alone. The default car is the single zero-ride-height baseline selected from
% carConfig(), because carConfig() otherwise returns its ride-height grid.
%
% Useful options:
%   initialFront / initialRear  starting grip vectors
%   maxIterations               adaptive passes, default 5
%   minSpan                     stop when both grid spans are at most this
%   shrinkFactor                span multiplier after each pass, default .5
%   timeTol                     relative time error unit, default .01
%   numWorkers                  gripSweep workers, default 8
%   refineGrid                  interpolation mesh, default 101
%   plot                        plot the final 2-D adaptive grid, default true
%   verbose                     print progress, default true
%   loadFromPrev                reuse matching gripSweep caches, default false
%   cachePrefix                 cache filename prefix
%   outputPath                  result MAT file, default adaptive_grip_sweep.mat
%   saveResults                 save result, default true
%
% carConfigFcn and sweepFcn are optional dependency-injection hooks used by
% tests and by callers that already have a configured model. The normal path
% uses @carConfig and @gripSweep.

if nargin < 1 || isempty(targetTimes)
    targetTimes = defaultTargets();
end
if nargin < 2 || isempty(options)
    options = struct();
end

targetTimes = validateTargets(targetTimes);
options = normalizeOptions(options);

entrypointDir = fileparts(mfilename('fullpath'));
fullCarModelsDir = fileparts(entrypointDir);
repoRoot = fileparts(fullCarModelsDir);
if options.bootstrap
    clear classes %#ok<CLCLS>
    run(fullfile(entrypointDir,'bootstrap.m'));
end

[baseCell,eventParams] = resolveBaseline(options);
if size(baseCell,1) ~= 1
    error('run_adaptive_grip_sweep:notSingleBaseline', ...
        'The adaptive sweep requires exactly one baseline car case.');
end

cachePrefix = options.cachePrefix;
if isempty(cachePrefix)
    cachePrefix = fullfile(repoRoot,'adaptive_grip_sweep');
end
outputPath = options.outputPath;
if isempty(outputPath)
    outputPath = fullfile(repoRoot,'adaptive_grip_sweep.mat');
end

front = options.initialFront;
rear = options.initialRear;
history = emptyHistory();
lastSweep = [];
candidate = [];
stopReason = 'maxIterations';

for iteration = 1:options.maxIterations
    sweepOptions = struct( ...
        'front',front, ...
        'rear',rear, ...
        'events',{{'skidpad','accel','autocross'}}, ...
        'target',targetTimes, ...
        'timeTol',options.timeTol, ...
        'numWorkers',options.numWorkers, ...
        'refineGrid',options.refineGrid, ...
        'eventParams',eventParams, ...
        'cacheName',sprintf('%s_iter_%02d.mat',char(cachePrefix),iteration), ...
        'loadFromPrev',options.loadFromPrev, ...
        'verbose',options.verbose);

    lastSweep = options.sweepFcn(baseCell,sweepOptions);
    candidate = selectCandidate(lastSweep);
    edge = edgeStatus(candidate,front,rear);

    h = struct( ...
        'iteration',iteration, ...
        'front',front, ...
        'rear',rear, ...
        'frontSpan',max(front)-min(front), ...
        'rearSpan',max(rear)-min(rear), ...
        'candidate',candidate, ...
        'onFrontEdge',edge.front, ...
        'onRearEdge',edge.rear, ...
        'onEdge',edge.any, ...
        'cacheName',sweepOptions.cacheName);
    history(end+1,1) = h; %#ok<AGROW>

    if options.verbose
        fprintf(['adaptive grip pass %d/%d: F %.5f-%.5f, R %.5f-%.5f, ' ...
            'candidate %.5f/%.5f, score %.4g\n'], ...
            iteration,options.maxIterations,min(front),max(front), ...
            min(rear),max(rear),candidate.gripF,candidate.gripR, ...
            candidate.score);
    end

    if edge.any
        stopReason = 'gridEdge';
        warning('run_adaptive_grip_sweep:gridEdge', ...
            ['The best grip candidate remains on the adaptive grid edge ' ...
             '(front=%d, rear=%d). Widen initialFront/initialRear before ' ...
             'trusting the fit.'],edge.front,edge.rear);
        break
    end

    if h.frontSpan <= options.minSpan && h.rearSpan <= options.minSpan
        stopReason = 'spanTolerance';
        break
    end

    if iteration < options.maxIterations
        front = shrinkAxis(front,candidate.gripF,options.shrinkFactor);
        rear = shrinkAxis(rear,candidate.gripR,options.shrinkFactor);
    end
end

if isempty(candidate)
    error('run_adaptive_grip_sweep:noCandidate', ...
        'The adaptive sweep did not produce a candidate.');
end

% The interpolated point is useful for locating the answer, but this solve is
% the result that should be carried into carConfig or another study.
confirmOptions = struct( ...
    'front',candidate.gripF, ...
    'rear',candidate.gripR, ...
    'events',{{'skidpad','accel','autocross'}}, ...
    'target',targetTimes, ...
    'timeTol',options.timeTol, ...
    'numWorkers',options.numWorkers, ...
    'refineGrid',1, ...
    'eventParams',eventParams, ...
    'cacheName',sprintf('%s_confirm.mat',char(cachePrefix)), ...
    'loadFromPrev',false, ...
    'verbose',options.verbose);
confirmationSweep = options.sweepFcn(baseCell,confirmOptions);
confirmed = selectCandidate(confirmationSweep);

result = struct();
result.target = targetTimes;
result.candidate = candidate;
result.confirmed = confirmed;
result.best = confirmed;
result.status = stopReason;
result.history = history;
result.lastSweep = lastSweep;
result.confirmationSweep = confirmationSweep;
result.settings = publicSettings(options,cachePrefix,outputPath);

if options.verbose
    fprintf('\nconfirmed grip fit: front %.5f, rear %.5f, score %.4g\n', ...
        confirmed.gripF,confirmed.gripR,confirmed.score);
    fprintf('  skidpad   sim %8.4f  target %8.4f\n', ...
        confirmed.skidpad,targetTimes.skidpad);
    fprintf('  accel     sim %8.4f  target %8.4f\n', ...
        confirmed.accel,targetTimes.accel);
    fprintf('  autocross sim %8.4f  target %8.4f\n', ...
        confirmed.autocross,targetTimes.autocross);
end

if options.plot
    plotGripSweep(lastSweep);
end

if options.saveResults
    save(char(outputPath),'result','-v7.3');
    if options.verbose
        fprintf('saved adaptive grip result: %s\n',char(outputPath));
    end
end
end

function targets = defaultTargets()
targets = struct('autocross',48,'skidpad',4.9,'accel',4.8);
end

function targets = validateTargets(targets)
if ~isstruct(targets) || numel(targets) ~= 1
    error('run_adaptive_grip_sweep:badTargets', ...
        'targetTimes must be a scalar struct.');
end

required = {'autocross','skidpad','accel'};
for k = 1:numel(required)
    name = required{k};
    if ~isfield(targets,name)
        error('run_adaptive_grip_sweep:missingTarget', ...
            'targetTimes.%s is required.',name);
    end
    value = targets.(name);
    if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value) || value <= 0
        error('run_adaptive_grip_sweep:badTarget', ...
            'targetTimes.%s must be a finite positive scalar in seconds.',name);
    end
end

targets = struct('autocross',double(targets.autocross), ...
    'skidpad',double(targets.skidpad),'accel',double(targets.accel));
end

function options = normalizeOptions(options)
if ~isstruct(options) || numel(options) ~= 1
    error('run_adaptive_grip_sweep:badOptions', ...
        'options must be a scalar struct.');
end

options.bootstrap = getOption(options,'bootstrap',true);
options.carConfigFcn = getOption(options,'carConfigFcn',@carConfig);
options.sweepFcn = getOption(options,'sweepFcn',@gripSweep);
options.initialFront = axisValues(getOption(options,'initialFront',0.58:0.01:0.66),'initialFront');
options.initialRear = axisValues(getOption(options,'initialRear',0.58:0.01:0.66),'initialRear');
options.maxIterations = getOption(options,'maxIterations',5);
options.minSpan = getOption(options,'minSpan',0.001);
options.shrinkFactor = getOption(options,'shrinkFactor',0.5);
options.timeTol = getOption(options,'timeTol',0.01);
options.numWorkers = getOption(options,'numWorkers',8);
options.refineGrid = getOption(options,'refineGrid',101);
options.plot = getOption(options,'plot',true);
options.verbose = getOption(options,'verbose',true);
options.loadFromPrev = getOption(options,'loadFromPrev',false);
options.cachePrefix = getOption(options,'cachePrefix','');
options.outputPath = getOption(options,'outputPath','');
options.saveResults = getOption(options,'saveResults',true);

if ~isscalar(options.maxIterations) || options.maxIterations < 1 || ...
        options.maxIterations ~= floor(options.maxIterations)
    error('run_adaptive_grip_sweep:badIterations', ...
        'maxIterations must be a positive integer.');
end
if ~isscalar(options.minSpan) || ~isfinite(options.minSpan) || options.minSpan <= 0
    error('run_adaptive_grip_sweep:badSpan', ...
        'minSpan must be a finite positive scalar.');
end
if ~isscalar(options.shrinkFactor) || ~isfinite(options.shrinkFactor) || ...
        options.shrinkFactor <= 0 || options.shrinkFactor >= 1
    error('run_adaptive_grip_sweep:badShrink', ...
        'shrinkFactor must be between zero and one.');
end
if ~isscalar(options.timeTol) || ~isfinite(options.timeTol) || options.timeTol <= 0
    error('run_adaptive_grip_sweep:badTimeTol', ...
        'timeTol must be a finite positive scalar.');
end
if ~isscalar(options.numWorkers) || options.numWorkers < 0 || ...
        options.numWorkers ~= floor(options.numWorkers)
    error('run_adaptive_grip_sweep:badWorkers', ...
        'numWorkers must be a nonnegative integer.');
end
if ~isscalar(options.refineGrid) || options.refineGrid < 1 || ...
        options.refineGrid ~= floor(options.refineGrid)
    error('run_adaptive_grip_sweep:badRefineGrid', ...
        'refineGrid must be a positive integer.');
end
end

function values = axisValues(values,name)
if ~isnumeric(values) || isempty(values) || any(~isfinite(values(:))) || any(values(:) <= 0)
    error('run_adaptive_grip_sweep:badAxis', ...
        '%s must contain finite positive grip multipliers.',name);
end
values = sort(unique(double(values(:).')));
if numel(values) < 2
    error('run_adaptive_grip_sweep:shortAxis', ...
        '%s must contain at least two grip values.',name);
end
end

function [baseCell,eventParams] = resolveBaseline(options)
if isfield(options,'baseCell') && ~isempty(options.baseCell)
    baseCell = options.baseCell;
    eventParams = getOption(options,'eventParams',[]);
    if isempty(eventParams)
        error('run_adaptive_grip_sweep:noEventParams', ...
            'eventParams is required when options.baseCell is supplied.');
    end
    return
end

[allCars,eventParams,designTable] = options.carConfigFcn();
if ~istable(designTable) || ...
        ~all(ismember({'static_front_ride_height_in','static_rear_ride_height_in'}, ...
                      designTable.Properties.VariableNames))
    error('run_adaptive_grip_sweep:noRideHeightColumns', ...
        'carConfigFcn must return the scalar-baseline ride-height columns.');
end

mask = abs(designTable.static_front_ride_height_in) <= 1e-12 & ...
       abs(designTable.static_rear_ride_height_in) <= 1e-12;
if nnz(mask) ~= 1
    error('run_adaptive_grip_sweep:baselineNotUnique', ...
        'Expected exactly one zero-ride-height baseline, found %d.',nnz(mask));
end
baseCell = allCars(mask,:);
end

function candidate = selectCandidate(G)
candidate = struct('gripF',NaN,'gripR',NaN,'score',NaN, ...
    'skidpad',NaN,'accel',NaN,'autocross',NaN,'source','');

if isfield(G,'bestInterp') && ~isempty(fieldnames(G.bestInterp)) && ...
        isfield(G.bestInterp,'score') && isfinite(G.bestInterp.score)
    source = 'interpolated';
    values = G.bestInterp;
elseif isfield(G,'best') && istable(G.best) && height(G.best) >= 1
    source = 'grid';
    row = G.best(1,:);
    values = struct('gripF',row.gripF,'gripR',row.gripR,'score',row.score, ...
        'skidpad',row.skidpad,'accel',row.accel,'autocross',row.autocross);
else
    error('run_adaptive_grip_sweep:noCandidate', ...
        'gripSweep did not return a usable best candidate.');
end

for name = {'gripF','gripR','score','skidpad','accel','autocross'}
    if isfield(values,name{1})
        candidate.(name{1}) = double(values.(name{1}));
    end
end
candidate.source = source;
if any(~isfinite([candidate.gripF,candidate.gripR,candidate.score]))
    error('run_adaptive_grip_sweep:noCandidate', ...
        'gripSweep returned a non-finite grip candidate or score.');
end
end

function edge = edgeStatus(candidate,front,rear)
tolF = 1e-9 + 1e-6*max(abs(front));
tolR = 1e-9 + 1e-6*max(abs(rear));
edge.front = candidate.gripF <= min(front)+tolF || candidate.gripF >= max(front)-tolF;
edge.rear = candidate.gripR <= min(rear)+tolR || candidate.gripR >= max(rear)-tolR;
edge.any = edge.front || edge.rear;
end

function values = shrinkAxis(values,center,shrinkFactor)
span = max(values)-min(values);
nextSpan = span*shrinkFactor;
halfSpan = nextSpan/2;
values = linspace(center-halfSpan,center+halfSpan,numel(values));
end

function history = emptyHistory()
history = struct('iteration',{},'front',{},'rear',{},'frontSpan',{}, ...
    'rearSpan',{},'candidate',{},'onFrontEdge',{},'onRearEdge',{}, ...
    'onEdge',{},'cacheName',{});
end

function settings = publicSettings(options,cachePrefix,outputPath)
settings = struct('initialFront',options.initialFront, ...
    'initialRear',options.initialRear,'maxIterations',options.maxIterations, ...
    'minSpan',options.minSpan,'shrinkFactor',options.shrinkFactor, ...
    'timeTol',options.timeTol,'numWorkers',options.numWorkers, ...
    'refineGrid',options.refineGrid,'cachePrefix',cachePrefix, ...
    'outputPath',outputPath);
end

function value = getOption(options,name,default)
if isfield(options,name)
    value = options.(name);
else
    value = default;
end
end
