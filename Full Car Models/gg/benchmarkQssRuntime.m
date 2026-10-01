function report = benchmarkQssRuntime(model,options)
%BENCHMARKQSSRUNTIME Measure one QSS configuration through selected events.
%   REPORT = BENCHMARKQSSRUNTIME(MODEL,OPTIONS) returns repeatable runtime and
%   validity evidence without writing benchmark artifacts. MODEL is forwarded
%   explicitly to carConfig. OPTIONS may contain repeats, ggOptions, and
%   events. The default is a three-repeat, production-grid, all-event run.

if nargin < 1 || isempty(model), model = "legacy"; end
if nargin < 2 || isempty(options), options = struct(); end
model = validateModel(model);
options = validateOptions(options);
eventNames = options.events;

phaseTemplate = struct('configurationSeconds',0,'ggSolveSeconds',0, ...
    'ggAssemblySeconds',0,'eventPreparationSeconds',0);
phaseSeconds = repmat(phaseTemplate,options.repeats,1);
eventSeconds = eventSecondsTemplate(eventNames,options.repeats);
eventTimes = eventTimesTemplate(eventNames,options.repeats);
repeatReports = repmat(repeatTemplate(),options.repeats,1);
totalSeconds = zeros(options.repeats,1);
wallSeconds = zeros(options.repeats,1);

for repeatIndex = 1:options.repeats
    wallTimer = tic;

    timer = tic;
    [cars,eventParams] = carConfig('FullFactorial',[],char(model));
    lapCar = cars{1,1};
    accelCar = cars{1,2};
    phaseSeconds(repeatIndex).configurationSeconds = toc(timer);

    timer = tic;
    rawGG = gg2(lapCar,0,options.ggOptions);
    phaseSeconds(repeatIndex).ggSolveSeconds = toc(timer);

    timer = tic;
    solvedCar = makeGG(rawGG,lapCar);
    phaseSeconds(repeatIndex).ggAssemblySeconds = toc(timer);

    timer = tic;
    eventSim = Events2(solvedCar,accelCar,eventParams);
    phaseSeconds(repeatIndex).eventPreparationSeconds = toc(timer);

    for eventIndex = 1:numel(eventNames)
        eventName = char(eventNames(eventIndex));
        methodName = eventName;
        methodName(1) = upper(methodName(1));
        timer = tic;
        eventSim.(methodName)();
        eventSeconds.(eventName)(repeatIndex) = toc(timer);
        eventTimes.(eventName)(repeatIndex) = eventTime(eventSim,eventName);
    end

    repeatReports(repeatIndex).counts = solverCounts(rawGG,solvedCar.ggMask);
    repeatReports(repeatIndex).solvers = solverEvidence(rawGG);
    repeatReports(repeatIndex).masks = solvedCar.ggMask;
    repeatReports(repeatIndex).eventTimes = eventTimeRow(eventTimes,eventNames,repeatIndex);
    repeatReports(repeatIndex).stationProfiles = stationProfiles(eventSim,eventNames);
    totalSeconds(repeatIndex) = sum(struct2array(phaseSeconds(repeatIndex))) + ...
        sumEventSeconds(eventSeconds,eventNames,repeatIndex);
    wallSeconds(repeatIndex) = toc(wallTimer);
end

report = struct();
report.timestamp = datetime('now');
report.configuration = struct('model',model,'ggOptions',options.ggOptions, ...
    'events',eventNames,'repeats',options.repeats);
report.environment = environmentEvidence();
report.inputs = inputEvidence();
report.phases = collapsePhaseSeconds(phaseSeconds);
report.eventSeconds = eventSeconds;
report.eventTimes = eventTimes;
report.totalSeconds = totalSeconds;
report.wallSeconds = wallSeconds;
report.repeats = repeatReports;
report.counts = sumCounts(repeatReports);
report.solvers = sumSolverEvidence(repeatReports);
report.masks = {repeatReports.masks}.';
report.stationProfiles = {repeatReports.stationProfiles}.';
end

function model = validateModel(model)
model = string(model);
if ~isscalar(model) || strlength(model) == 0
    error('benchmarkQssRuntime:badModel','model must be one nonempty string.');
end
end

function options = validateOptions(options)
if ~isstruct(options) || ~isscalar(options)
    error('benchmarkQssRuntime:badOptions','options must be a scalar struct.');
end
options.repeats = optionOr(options,'repeats',3);
validateattributes(options.repeats,{'numeric'}, ...
    {'scalar','integer','positive','finite'},mfilename,'options.repeats');
if isfield(options,'events')
    options.events = string(options.events);
else
    options.events = ["skidpad","accel","autocross","endurance"];
end
if ~isempty(options.events) && ~isvector(options.events)
    error('benchmarkQssRuntime:badEvents','events must be a string vector.');
end
options.events = options.events(:).';
allowedEvents = ["skidpad","accel","autocross","endurance"];
if any(~ismember(options.events,allowedEvents)) || numel(unique(options.events)) ~= numel(options.events)
    error('benchmarkQssRuntime:badEvents', ...
        'events must contain each supported event at most once.');
end
productionOptions = ggProductionOptions();
options.ggOptions = optionOr(options,'ggOptions',productionOptions);
if ~isstruct(options.ggOptions) || ~isscalar(options.ggOptions)
    error('benchmarkQssRuntime:badGGOptions','ggOptions must be a scalar struct.');
end
if ~isfield(options.ggOptions,'fastScreening')
    options.ggOptions.fastScreening = false;
end
if ~isfield(options.ggOptions,'continuation')
    options.ggOptions.continuation = productionOptions.continuation;
end
if ~isfield(options.ggOptions,'maxVelocity')
    options.ggOptions.maxVelocity = productionOptions.maxVelocity;
end
end

function template = repeatTemplate()
template = struct('counts',emptyCounts(),'solvers',emptySolverEvidence(), ...
    'masks',struct(),'eventTimes',struct(),'stationProfiles',struct());
end

function counts = emptyCounts()
counts = struct('requested',0,'accepted',0, ...
    'lateral',struct('requested',0,'accepted',0), ...
    'acceleration',struct('requested',0,'accepted',0), ...
    'braking',struct('requested',0,'accepted',0));
end

function counts = solverCounts(paramArr,mask)
requested = numel(paramArr);
counts = emptyCounts();
counts.lateral = struct('requested',requested,'accepted',nnz(mask.lateral));
counts.acceleration = struct('requested',requested,'accepted',nnz(mask.acceleration));
counts.braking = struct('requested',requested,'accepted',nnz(mask.braking));
counts.requested = counts.lateral.requested + counts.acceleration.requested + ...
    counts.braking.requested;
counts.accepted = counts.lateral.accepted + counts.acceleration.accepted + ...
    counts.braking.accepted;
end

function evidence = emptySolverEvidence()
branch = struct('functionEvaluations',0,'equationCalls',0,'residuals',zeros(0,1));
evidence = struct('lateral',branch,'acceleration',branch,'braking',branch, ...
    'functionEvaluations',struct('total',0),'equationCalls',struct('total',0));
end

function evidence = solverEvidence(paramArr)
evidence = emptySolverEvidence();
arr = paramArr(:);
evidence.lateral = branchEvidence(arr,'maxLat');
evidence.acceleration = branchEvidence(arr,'maxLong');
evidence.braking = branchEvidence(arr,'maxBrake');
evidence.functionEvaluations.total = evidence.lateral.functionEvaluations + ...
    evidence.acceleration.functionEvaluations + evidence.braking.functionEvaluations;
evidence.equationCalls.total = evidence.lateral.equationCalls + ...
    evidence.acceleration.equationCalls + evidence.braking.equationCalls;
end

function evidence = branchEvidence(arr,prefix)
evaluationField = [prefix 'FunctionEvaluations'];
equationField = [prefix 'EquationCalls'];
residualField = [prefix 'Residual'];
evidence = struct('functionEvaluations',0,'equationCalls',0,'residuals',zeros(0,1));
if isempty(arr) || ~all(isprop(arr,evaluationField))
    return
end
evidence.functionEvaluations = sum([arr.(evaluationField)]);
if all(isprop(arr,equationField))
    evidence.equationCalls = sum([arr.(equationField)]);
end
if all(isprop(arr,residualField))
    residuals = [arr.(residualField)].';
    evidence.residuals = residuals(isfinite(residuals));
end
end

function totals = sumCounts(repeats)
totals = emptyCounts();
for index = 1:numel(repeats)
    totals.requested = totals.requested + repeats(index).counts.requested;
    totals.accepted = totals.accepted + repeats(index).counts.accepted;
    names = ["lateral","acceleration","braking"];
    for name = names
        totals.(name).requested = totals.(name).requested + repeats(index).counts.(name).requested;
        totals.(name).accepted = totals.(name).accepted + repeats(index).counts.(name).accepted;
    end
end
end

function totals = sumSolverEvidence(repeats)
totals = emptySolverEvidence();
for index = 1:numel(repeats)
    names = ["lateral","acceleration","braking"];
    for name = names
        totals.(name).functionEvaluations = totals.(name).functionEvaluations + ...
            repeats(index).solvers.(name).functionEvaluations;
        totals.(name).equationCalls = totals.(name).equationCalls + ...
            repeats(index).solvers.(name).equationCalls;
        totals.(name).residuals = [totals.(name).residuals; ...
            repeats(index).solvers.(name).residuals]; %#ok<AGROW>
    end
end
totals.functionEvaluations.total = totals.lateral.functionEvaluations + ...
    totals.acceleration.functionEvaluations + totals.braking.functionEvaluations;
totals.equationCalls.total = totals.lateral.equationCalls + ...
    totals.acceleration.equationCalls + totals.braking.equationCalls;
end

function seconds = collapsePhaseSeconds(repeats)
seconds = struct();
names = fieldnames(repeats);
for index = 1:numel(names)
    name = names{index};
    seconds.(name) = reshape([repeats.(name)],[],1);
end
end

function values = eventSecondsTemplate(eventNames,repeats)
values = struct();
for eventName = eventNames
    values.(eventName) = zeros(repeats,1);
end
end

function values = eventTimesTemplate(eventNames,repeats)
values = struct();
for eventName = eventNames
    values.(eventName) = nan(repeats,1);
end
end

function value = eventTime(eventSim,eventName)
if isprop(eventSim,'times') && isfield(eventSim.times,eventName)
    value = eventSim.times.(eventName);
else
    value = NaN;
end
end

function row = eventTimeRow(eventTimes,eventNames,index)
row = struct();
for eventName = eventNames
    row.(eventName) = eventTimes.(eventName)(index);
end
end

function profiles = stationProfiles(eventSim,eventNames)
profiles = struct();
for eventName = intersect(eventNames,["autocross","endurance"],'stable')
    eventData = eventSim.(eventName);
    profiles.(eventName) = struct('longVel',eventData.long_vel, ...
        'longAccel',eventData.long_accel,'latAccel',eventData.lat_accel);
end
end

function seconds = sumEventSeconds(eventSeconds,eventNames,index)
seconds = 0;
for eventName = eventNames
    seconds = seconds + eventSeconds.(eventName)(index);
end
end

function evidence = environmentEvidence()
evidence = struct('matlabVersion',string(version), ...
    'computer',string(computer),'cores',safeCoreCount(), ...
    'gitRevision',gitOutput('rev-parse HEAD'), ...
    'dirtyFiles',splitlines(gitOutput('status --porcelain')));
end

function count = safeCoreCount()
try
    count = feature('numcores');
catch
    count = NaN;
end
end

function evidence = inputEvidence()
names = {'carConfig','gg2','makeGG','Events2','michigantrack2024.mat','2024endurancetrack.mat'};
evidence = repmat(struct('name',"",'path',"",'sha256',""),numel(names),1);
for index = 1:numel(names)
    path = which(names{index});
    evidence(index).name = string(names{index});
    evidence(index).path = string(path);
    evidence(index).sha256 = fileHash(path);
end
end

function output = gitOutput(arguments)
[status,output] = system(['git -c safe.directory=C:/VD ' arguments]);
if status ~= 0
    output = '';
end
output = string(strtrim(output));
end

function output = fileHash(path)
output = "";
if isempty(path) || ~isfile(path)
    return
end
fileId = fopen(path,'r');
cleanup = onCleanup(@() fclose(fileId)); %#ok<NASGU>
bytes = fread(fileId,Inf,'*uint8');
digest = java.security.MessageDigest.getInstance('SHA-256');
digest.update(bytes);
output = string(lower(reshape(dec2hex(uint8(digest.digest()),2).',1,[])));
end

function value = optionOr(options,name,defaultValue)
if isfield(options,name) && ~isempty(options.(name))
    value = options.(name);
else
    value = defaultValue;
end
end
