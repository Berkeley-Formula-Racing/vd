function state = runAdaptiveDOE(study)
%RUNADAPTIVEDOE Run or resume a checkpointed adaptive vehicle DOE study.

requireFitrgp()
study = applyParallelGuard(study);
outputDirectory = string(study.output.directory);
if ~isfolder(outputDirectory), mkdir(char(outputDirectory)); end
checkpointPath = fullfile(outputDirectory,string(study.output.checkpoint));
resultsPath = fullfile(outputDirectory,string(study.output.results));
job = simLog.start('adaptive DOE');
started = tic;

try
    [~,eventParams,~,baselineTable] = carConfig();
    resolverStudy = study;
    % Task 1 validates worker counts for parallel execution, while Task 7
    % defines zero as the serial-run sentinel.  Bounds resolution itself is
    % independent of workers, so preserve the execution setting below.
    if resolverStudy.numWorkers == 0, resolverStudy.numWorkers = 1; end
    resolved = doeResolveStudy(resolverStudy,baselineTable);
    resolvedStudy = attachStudySettings(resolved,study);

    if study.resume && isfile(checkpointPath)
        [state,~] = loadDOECheckpoint(checkpointPath,resolvedStudy);
        validateResumeStudy(state.resolvedStudy,resolvedStudy)
        state = normalizeState(state,study,eventParams);
    else
        state = newState(resolvedStudy,study,eventParams);
        [state.pendingU,physicalDesign] = doeInitialDesign( ...
            resolvedStudy,study.initialCases,study.randomSeed);
        state.pendingSelection = initialSelection(height(physicalDesign));
        saveDOECheckpoint(checkpointPath,state)
    end

    elapsedBefore = state.elapsed;
    state = backfillRamps(state,checkpointPath);
    while height(state.designTable) < state.resolvedStudy.maxCases
        if isempty(state.pendingU)
            state = scheduleNextBatch(state,checkpointPath);
        end

        remaining = state.resolvedStudy.maxCases-height(state.designTable);
        n = min([state.resolvedStudy.batchSize,remaining,size(state.pendingU,1)]);
        batchU = state.pendingU(1:n,:);
        batchSelection = state.pendingSelection(1:n,:);
        state.pendingU(1:n,:) = [];
        state.pendingSelection(1:n,:) = [];

        batchTable = state.resolvedStudy.toPhysical(batchU);
        batchCars = carConfig("Explicit",batchTable);
        caseIndices = height(state.designTable) + (1:n)';
        results = doeRunBatch(batchCars,state.eventParams,state.resolvedStudy,caseIndices);
        state = appendBatch(state,batchU,batchTable,batchSelection,results);
        state.batchNumber = state.batchNumber + 1;
        state.elapsed = elapsedBefore + toc(started);
        saveDOECheckpoint(checkpointPath,state)
    end

    state.elapsed = elapsedBefore + toc(started);
    state.study = study;
    writeResults(resultsPath,state)
    finishLog(job,state,resultsPath,"complete")
catch ME
    if exist('state','var')
        try
            state.elapsed = elapsedBefore + toc(started);
            finishLog(job,state,resultsPath,"failed")
        catch
        end
    else
        try
            simLog.finish(job,'details',"status=failed; " + string(ME.identifier));
        catch
        end
    end
    rethrow(ME)
end
end

function state = newState(resolvedStudy,study,eventParams)
nParameters = height(resolvedStudy.parameters);
state = struct();
state.resolvedStudy = resolvedStudy;
state.study = study;
state.U = zeros(0,nParameters);
state.designTable = resolvedStudy.toPhysical(zeros(0,nParameters));
state.carCell = cell(0,2);
state.metricTable = emptyMetricTable();
state.rampData = cell(0,1);
state.pointData = cell(0,1);
state.caseStatus = strings(0,1);
state.rampBackfillDiagnostics = emptyRampBackfillDiagnostics();
state.selectionHistory = emptySelectionTable();
state.pendingU = zeros(0,nParameters);
state.pendingSelection = emptySelectionTable();
state.batchNumber = 0;
state.elapsed = 0;
state.randomState = seededState(study.randomSeed);
state.eventParams = eventParams;
state.parallelFallback = logical(getField(study,'parallelFallback',false));
end

function state = normalizeState(state,study,eventParams)
state.study = study;
if ~isfield(state,'eventParams'), state.eventParams = eventParams; end
if ~isfield(state,'rampData'), state.rampData = cell(height(state.designTable),1); end
if ~isfield(state,'pointData'), state.pointData = cell(height(state.designTable),1); end
if ~isfield(state,'caseStatus'), state.caseStatus = repmat("complete",height(state.designTable),1); end
if ~isfield(state,'rampBackfillDiagnostics')
    state.rampBackfillDiagnostics = emptyRampBackfillDiagnostics();
end
if ~isfield(state,'selectionHistory'), state.selectionHistory = emptySelectionTable(); end
if ~isfield(state,'pendingU'), state.pendingU = zeros(0,size(state.U,2)); end
if ~isfield(state,'pendingSelection'), state.pendingSelection = emptySelectionTable(); end
if ~isfield(state,'batchNumber'), state.batchNumber = 0; end
if ~isfield(state,'elapsed'), state.elapsed = 0; end
if ~isfield(state,'randomState'), state.randomState = seededState(study.randomSeed); end
if ~isfield(state,'parallelFallback'), state.parallelFallback = logical(getField(study,'parallelFallback',false)); end
end

function validateResumeStudy(savedStudy,requestedStudy)
allowed = ["mode","maxCases","numWorkers","batchSize","objective","ramps"];
coveredBySignature = ["parameters","toPhysical","signature"];
runtimeOnly = ["resume","parallelFallback"];
protected = setdiff(string(fieldnames(requestedStudy)), ...
    [allowed,coveredBySignature,runtimeOnly]);
for fieldName = protected'
    name = char(fieldName);
    if ~isfield(savedStudy,name) || ~isequaln(savedStudy.(name),requestedStudy.(name))
        error('runAdaptiveDOE:resumeStudyMismatch', ...
            'Cannot change study.%s while resuming; only mode, maxCases, numWorkers, batchSize, objective, and ramps may change.', ...
            name)
    end
end
end

function state = backfillRamps(state,checkpointPath)
if ~state.resolvedStudy.ramps.enabled || isempty(state.carCell), return, end
missing = false(height(state.designTable),1);
for i = 1:height(state.designTable)
    missing(i) = state.caseStatus(i) == "complete" && ...
        (i > numel(state.rampData) || isempty(state.rampData{i}));
end
if ~any(missing), return, end

indices = find(missing);
results = doeRunBatch(state.carCell(indices,:),state.eventParams, ...
    state.resolvedStudy,indices,"rampOnly");
successful = zeros(0,1);
for j = 1:numel(indices)
    i = indices(j);
    result = results{j};
    if result.status == "complete"
        state = mergeSuccessfulRamp(state,i,result);
        successful(end+1,1) = i; %#ok<AGROW>
    else
        state.rampBackfillDiagnostics = [state.rampBackfillDiagnostics; ...
            rampBackfillDiagnostic(i,result)];
    end
end
state = recomputeScores(state,successful);
saveDOECheckpoint(checkpointPath,state)
end

function state = mergeSuccessfulRamp(state,index,result)
state.rampData{index,1} = rampCache(result);
rampFields = ["understeer_gradient_10_deg_per_g", ...
    "understeer_gradient_25_deg_per_g","rebalance_speed_mps"];
for fieldName = rampFields
    name = char(fieldName);
    if ismember(name,result.metricRow.Properties.VariableNames) && ...
            ismember(name,state.metricTable.Properties.VariableNames)
        state.metricTable.(name)(index) = result.metricRow.(name);
    end
end
end

function row = rampBackfillDiagnostic(caseIndex,result)
row = table(caseIndex,string(result.status),string(result.errorIdentifier), ...
    string(result.errorMessage),result.elapsed, ...
    'VariableNames',{'case_index','status','error_identifier', ...
    'error_message','elapsed_s'});
end

function state = scheduleNextBatch(state,checkpointPath)
remaining = state.resolvedStudy.maxCases-height(state.designTable);
n = min(state.resolvedStudy.batchSize,remaining);
saveDOECheckpoint(checkpointPath,state)

validCount = nnz(state.metricTable.valid);
if validCount <= size(state.U,2)
    [state.pendingU,state.pendingSelection,state.randomState] = ...
        spaceFillingBatch(state,n);
else
    selectionStudy = state.resolvedStudy;
    selectionStudy.batchSize = n;
    [state.pendingU,state.pendingSelection] = ...
        doeSelectAdaptiveBatch(state,selectionStudy);
    nextState = state.pendingSelection.Properties.UserData.randomState;
    if ~isempty(nextState), state.randomState = nextState; end
    state.pendingSelection.Properties.UserData = [];
end
saveDOECheckpoint(checkpointPath,state)
end

function state = appendBatch(state,batchU,batchTable,batchSelection,results)
n = numel(results);
batchCars = cell(n,2);
rows = cell(n,1);
ramps = cell(n,1);
points = cell(n,1);
status = strings(n,1);
for j = 1:n
    result = results{j};
    batchCars(j,:) = {result.car,result.accelCar};
    rows{j} = result.metricRow;
    ramps{j} = rampCache(result);
    points{j} = result.points;
    status(j) = result.status;
end

state.U = [state.U; batchU];
state.designTable = [state.designTable; batchTable];
state.carCell = [state.carCell; batchCars];
state.metricTable = [state.metricTable; vertcat(rows{:})];
state.rampData = [state.rampData; ramps];
state.pointData = [state.pointData; points];
state.caseStatus = [state.caseStatus; status];
state.selectionHistory = [state.selectionHistory; batchSelection];
state = recomputeScores(state,(height(state.designTable)-n+1:height(state.designTable))');
end

function state = recomputeScores(state,indices)
for k = 1:numel(indices)
    i = indices(k);
    points = cachedPoints(state,i);
    [~,breakdown] = doeScoreCase(points,state.metricTable(i,:), ...
        state.resolvedStudy.objective);
    names = fieldnames(breakdown);
    for j = 1:numel(names)
        state.metricTable.(names{j})(i) = breakdown.(names{j});
    end
end
end

function points = cachedPoints(state,index)
points = emptyPoints();
if index <= numel(state.pointData) && ~isempty(state.pointData{index})
    points = state.pointData{index};
    return
end
try
    points = state.carCell{index,1}.comp.points;
catch
end
end

function cache = rampCache(result)
cache = struct('summary',result.rampSummary,'points',result.rampPoints);
if isempty(cache.summary) && isempty(cache.points), cache = []; end
end

function [U,selection,nextRandomState] = spaceFillingBatch(state,n)
candidateCount = max(1000,20*n);
previous = rng;
restore = onCleanup(@() rng(previous)); %#ok<NASGU>
rng(state.randomState)
candidates = rand(candidateCount,size(state.U,2));
nextRandomState = rng;

U = zeros(n,size(state.U,2));
nearest = zeros(n,1);
reference = state.U;
for j = 1:n
    distances = nearestDistances(candidates,reference);
    [nearest(j),pick] = max(distances);
    U(j,:) = candidates(pick,:);
    reference = [reference; U(j,:)]; %#ok<AGROW>
    candidates(pick,:) = [];
end
selection = table(repmat(string(state.resolvedStudy.mode),n,1), ...
    repmat("space_filling",n,1),nearest,nearest,nan(n,1), ...
    'VariableNames',{'mode','source','acquisition_value', ...
    'nearest_existing_distance','candidate_index'});
end

function distances = nearestDistances(points,reference)
if isempty(reference)
    distances = inf(size(points,1),1);
    return
end
distances = inf(size(points,1),1);
for i = 1:size(points,1)
    distances(i) = min(sqrt(mean((reference-points(i,:)).^2,2)));
end
end

function selection = initialSelection(n)
selection = table(repmat("initial",n,1),repmat("initial",n,1), ...
    nan(n,1),inf(n,1),nan(n,1), ...
    'VariableNames',{'mode','source','acquisition_value', ...
    'nearest_existing_distance','candidate_index'});
end

function T = emptySelectionTable()
T = table(strings(0,1),strings(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
    'VariableNames',{'mode','source','acquisition_value', ...
    'nearest_existing_distance','candidate_index'});
end

function T = emptyRampBackfillDiagnostics()
T = table(zeros(0,1),strings(0,1),strings(0,1),strings(0,1),zeros(0,1), ...
    'VariableNames',{'case_index','status','error_identifier', ...
    'error_message','elapsed_s'});
end

function T = emptyMetricTable()
[row,~] = doeCaseMetrics([],NaN,[]);
score = struct('points_skidpad',NaN,'points_accel',NaN, ...
    'points_autocross',NaN,'points_endurance',NaN, ...
    'modeled_dynamic_points',NaN,'penalty_invalid_case',NaN, ...
    'penalty_solve_failure',NaN,'penalty_wheel_lift',NaN, ...
    'penalty_energy',NaN,'penalty_understeer',NaN, ...
    'penalty_rebalance',NaN,'total_penalty_points',NaN, ...
    'objective_score',NaN);
T = [row struct2table(score)];
T(1,:) = [];
end

function resolvedStudy = attachStudySettings(resolved,study)
resolvedStudy = study;
resolvedStudy.parameters = resolved.parameters;
resolvedStudy.toPhysical = resolved.toPhysical;
resolvedStudy.signature = resolved.signature;
end

function study = applyParallelGuard(study)
study.parallelFallback = false;
if study.numWorkers == 0 || hasParallelToolbox(), return, end
if study.allowSerialFallback
    study.numWorkers = 0;
    study.parallelFallback = true;
else
    error('runAdaptiveDOE:noParallelToolbox', ...
        'Parallel Computing Toolbox is required when study.numWorkers is greater than zero.')
end
end

function tf = hasParallelToolbox()
try
    tf = license('test','Distrib_Computing_Toolbox');
catch
    tf = false;
end
end

function requireFitrgp()
if exist('fitrgp','file') == 0
    error('runAdaptiveDOE:noFitrgp', ...
        'fitrgp is required before an adaptive DOE study can run.')
end
end

function state = seededState(seed)
previous = rng;
restore = onCleanup(@() rng(previous)); %#ok<NASGU>
rng(seed,'twister')
state = rng;
end

function value = getField(s,name,default)
if isfield(s,name), value = s.(name); else, value = default; end
end

function points = emptyPoints()
points = struct('skidpad',NaN,'accel',NaN,'autocross',NaN, ...
    'endurance',NaN,'total',NaN);
end

function writeResults(path,state)
carCell = state.carCell; %#ok<NASGU>
designTable = state.designTable; %#ok<NASGU>
eventParams = state.eventParams; %#ok<NASGU>
metricTable = state.metricTable; %#ok<NASGU>
rampData = state.rampData; %#ok<NASGU>
study = state.study; %#ok<NASGU>
selectionHistory = state.selectionHistory; %#ok<NASGU>
rampBackfillDiagnostics = state.rampBackfillDiagnostics; %#ok<NASGU>
save(char(path),'carCell','designTable','eventParams','metricTable', ...
    'rampData','study','selectionHistory','rampBackfillDiagnostics','-v7.3')
end

function finishLog(job,state,resultsPath,status)
validCount = nnz(state.metricTable.valid);
invalidCount = height(state.metricTable)-validCount;
backfillFailures = height(state.rampBackfillDiagnostics);
details = sprintf(['status=%s; mode=%s; batches=%d; valid=%d; invalid=%d; ' ...
    'ramp_backfill_failures=%d; checkpoint=%s; results=%s'], ...
    status,string(state.resolvedStudy.mode),state.batchNumber,validCount, ...
    invalidCount,backfillFailures, ...
    fullfile(string(state.resolvedStudy.output.directory), ...
    string(state.resolvedStudy.output.checkpoint)),resultsPath);
simLog.finish(job,'events',cellstr(string(state.resolvedStudy.events)), ...
    'workers',state.resolvedStudy.numWorkers,'nCases',height(state.designTable), ...
    'details',details)
end
