function [study,events] = runStudy(cars,cases,request,callbacks)
%RUNSTUDY Execute an ordered, checkpointed ramp-speed study.

if nargin < 2 || isempty(cases)
    cases = repmat(emptyCase(),0,1);
end
if nargin < 3 || isempty(request)
    request = struct();
end
if nargin < 4 || isempty(callbacks)
    callbacks = struct();
end
validateInputs(cases,request,callbacks);
request = normalizeRequest(request);
callbacks = normalizeCallbacks(callbacks);
cases = normalizeCases(cases);
type = normalizeType(request.rampType);

study = rampSpeed.makeStudy(request.appVersion);
study.cases = cases;
study.status = "running";
started = datetime('now');
study.runMeta = makeStudyMeta(started,type,request);
study.runs = initializeRuns(cases,type,request.settings);
events = repmat(emptyEvent(),0,1);
totalCases = numel(cases);
workerCancelFiles = strings(1,totalCases);
completedCases = 0;
effectiveWorkers = 0;
fallbackReason = "";
activeFutures = [];
activeIndices = [];

if totalCases == 0
    study.status = "complete";
    study.runMeta.status = study.status;
    study.runMeta.completed = datetime('now');
    saveCheckpoint();
    return
end

[effectiveWorkers,fallbackReason,useParallel,pool] = ...
    chooseExecution(request);
study.runMeta.effectiveWorkers = effectiveWorkers;
study.runMeta.parallelFallbackReason = fallbackReason;
applyExecutionMetadata();
if strlength(fallbackReason) > 0
    addStudyWarning(fallbackReason);
    emitEvent(makeEvent("study","",NaN,NaN,completedCases,totalCases, ...
        fallbackReason));
end

if useParallel
    try
        runParallelCases(pool);
    catch ME
        interruptParallelCases(ME);
        reason = "Parallel execution failed; falling back to serial execution: " + ...
            string(ME.message);
        effectiveWorkers = 0;
        fallbackReason = reason;
        study.runMeta.effectiveWorkers = 0;
        study.runMeta.parallelFallbackReason = reason;
        applyExecutionMetadata();
        addStudyWarning(reason);
        emitEvent(makeEvent("study","",NaN,NaN,completedCases,totalCases, ...
            reason));
        runSerialCases(findPendingCases());
    end
else
    runSerialCases(1:totalCases);
end

if study.status == "running"
    statusValues = strings(numel(study.runs),1);
    for i = 1:numel(study.runs)
        statusValues(i) = string(study.runs(i).status);
    end
    if any(statusValues == "cancelled")
        study.status = "cancelled";
    elseif all(statusValues == "complete")
        study.status = "complete";
    elseif any(statusValues == "complete")
        study.status = "partial";
    else
        study.status = "failed";
    end
end
study.runMeta.completed = datetime('now');
study.runMeta.status = study.status;
saveCheckpoint();
cleanupWorkerCancellationFiles();

    function runSerialCases(indices)
        indices = double(indices(:)).';
        for index = indices
            if isCancelledRequested()
                cancelRemainingCases(index:totalCases);
                study.status = "cancelled";
                emitEvent(makeEvent("study","",NaN,NaN,completedCases, ...
                    totalCases,"study cancelled"));
                saveCheckpoint();
                return
            end

            emitEvent(makeEvent("case",caseId(index),NaN,NaN, ...
                completedCases,totalCases,"starting " + caseLabel(index)));
            run = executeOneCase(index);
            study.runs(index) = finalizeRun(run,index,study.runs(index), ...
                type,request.settings,cases(index),effectiveWorkers, ...
                fallbackReason,started,request);
            completedCases = completedCases + 1;
            emitEvent(makeEvent("case",caseId(index),NaN,NaN, ...
                completedCases,totalCases, ...
                string(study.runs(index).status) + ": " + caseLabel(index)));
            saveCheckpoint();

            if string(study.runs(index).status) == "cancelled" || ...
                    isCancelledRequested()
                cancelRemainingCases(index+1:totalCases);
                study.status = "cancelled";
                emitEvent(makeEvent("study","",NaN,NaN,completedCases, ...
                    totalCases,"study cancelled"));
                saveCheckpoint();
                return
            end
        end
    end

    function runParallelCases(activePool)
        queue = request.progressQueue;
        if isempty(queue)
            queue = parallel.pool.DataQueue;
            afterEach(queue,@handleQueuedProgress);
        end
        workerRequest = request;
        workerRequest.progressQueue = [];
        workerRequest.totalCases = totalCases;
        workerRequest.cancelFile = "";
        parallelCars = cell(1,totalCases);
        parallelRequests = cell(1,totalCases);
        probeCallbacks = struct("onProgress",@noopProgress, ...
            "isCancelled",@alwaysFalse);
        prepared = 0;
        for index = 1:totalCases
            if isCancelledRequested()
                cancelRemainingCases(index:totalCases);
                study.status = "cancelled";
                emitEvent(makeEvent("study","",NaN,NaN,completedCases, ...
                    totalCases,"study cancelled"));
                break
            end
            car = selectCarForCase(cases(index),index);
            caseRequest = workerRequest;
            caseRequest.cancelFile = string(tempname) + ".cancel";
            workerCancelFiles(index) = caseRequest.cancelFile;
            parallelCars{index} = car;
            parallelRequests{index} = caseRequest;
            prepared = index;
        end
        if prepared < totalCases
            return
        end
        preflightParallelInputs(activePool,request.runCaseFcn,parallelCars, ...
            cases,parallelRequests,probeCallbacks,queue);
        futures = parallel.FevalFuture.empty(0,totalCases);
        submitted = 0;
        for index = 1:prepared
            futures(index) = parfeval(activePool,@runCaseWorker,1, ...
                parallelCars{index},cases(index), ...
                parallelRequests{index},queue);
            submitted = index;
            activeIndices = 1:submitted;
            activeFutures = futures(1:submitted);
        end
        active = futures(1:submitted);
        activeIndices = 1:submitted;
        activeFutures = active;
        while ~isempty(active)
            [completedIndex,run] = fetchNext(active);
            index = activeIndices(completedIndex);
            active(completedIndex) = [];
            activeIndices(completedIndex) = [];
            activeFutures = active;
            study.runs(index) = finalizeRun(run,index,study.runs(index), ...
                type,request.settings,cases(index),effectiveWorkers, ...
                fallbackReason,started,request);
            completedCases = completedCases + 1;
            emitEvent(makeEvent("case",caseId(index),NaN,NaN, ...
                completedCases,totalCases, ...
                string(study.runs(index).status) + ": " + caseLabel(index)));
            saveCheckpoint();
            if isCancelledRequested()
                signalWorkerCancellation();
                cancelCause = MException("rampSpeed:cancelled", ...
                    "study cancelled");
                interruptParallelCases(cancelCause);
                cancelRemainingCases(findPendingCases());
                study.status = "cancelled";
                emitEvent(makeEvent("study","",NaN,NaN,completedCases, ...
                    totalCases,"study cancelled"));
                return
            end
        end
        activeFutures = [];
        activeIndices = [];
    end

    function interruptParallelCases(cause)
        futures = activeFutures;
        indices = activeIndices;
        activeFutures = [];
        activeIndices = [];
        if isempty(futures)
            return
        end
        signalWorkerCancellation();
        settleDeadline = tic;
        while ~allFuturesSettled(futures) && toc(settleDeadline) < 2
            pause(0.05);
        end
        try
            cancel(futures);
        catch
        end
        for k = 1:numel(futures)
            try
                wait(futures(k));
            catch
            end
        end
        message = "parallel execution interrupted: " + string(cause.message);
        for k = 1:numel(indices)
            index = indices(k);
            if index < 1 || index > totalCases || ...
                    ~any(string(study.runs(index).status) == ["pending","running"])
                continue
            end
            workerError = [];
            state = "";
            try
                state = string(futures(k).State);
                candidate = futures(k).Error;
                if ~isempty(candidate)
                    workerError = candidate;
                end
            catch
            end
            if state == "failed" && ~isCancellationError(cause)
                if isempty(workerError)
                    workerError = cause;
                end
                run = failureRun(cases(index),workerError,type,request.settings);
            else
                run = cancelledRun(cases(index),type,request.settings, ...
                    message);
            end
            study.runs(index) = finalizeRun(run,index,study.runs(index), ...
                type,request.settings,cases(index),effectiveWorkers, ...
                fallbackReason,started,request);
        end
    end
    function value = allFuturesSettled(futures)
        value = true;
        for futureIndex = 1:numel(futures)
            try
                state = string(futures(futureIndex).State);
            catch
                value = false;
                return
            end
            if ~any(state == ["finished","failed","cancelled", ...
                    "unavailable"])
                value = false;
                return
            end
        end
    end
    function run = executeOneCase(index)
        caseInfo = cases(index);
        try
            car = selectCarForCase(caseInfo,index);
            caseCallbacks = callbacks;
            caseCallbacks.onProgress = @(rawEvent) ...
                handleCaseProgress(rawEvent,index);
            caseCallbacks.isCancelled = @()isCancelledRequested();
            run = feval(request.runCaseFcn,car,caseInfo,request, ...
                caseCallbacks);
        catch ME
            if isCancellationError(ME)
                run = cancelledRun(caseInfo,type,request.settings, ...
                    "study cancelled during case: " + string(ME.message));
            else
                run = failureRun(caseInfo,ME,type,request.settings);
            end
        end
    end

    function car = selectCarForCase(caseInfo,index)
        role = selectedRole(caseInfo,type,request);
        column = 1;
        if role == "acceleration"
            column = 2;
        end
        sourceIndex = numericField(caseInfo,"sourceIndex", ...
            numericField(caseInfo,"designRow",index));
        if ~isfinite(sourceIndex) || sourceIndex < 1
            sourceIndex = index;
        end
        sourceIndex = round(sourceIndex);
        if iscell(cars)
            [nRows,nCols] = size(cars);
            if nRows >= sourceIndex && nCols >= column
                car = cars{sourceIndex,column};
                return
            end
            if nRows == 1 && nCols >= column
                car = cars{1,column};
                return
            end
            if isvector(cars) && numel(cars) >= index
                car = cars{index};
                return
            end
        elseif isstruct(cars)
            if isscalar(cars)
                car = cars;
                return
            elseif numel(cars) >= sourceIndex
                car = cars(sourceIndex);
                return
            end
        elseif isscalar(cars)
            car = cars;
            return
        end
        error("rampSpeed:missingCar", ...
            "No car is available for case %s.",caseId(index));
    end

    function handleCaseProgress(rawEvent,index)
        emitEvent(makeEvent("case",caseId(index), ...
            numericField(rawEvent,"speedIndex",NaN), ...
            numericField(rawEvent,"speed_mps",NaN),completedCases, ...
            totalCases,stringField(rawEvent,"message", ...
            "running " + caseLabel(index))));
    end

    function handleQueuedProgress(rawEvent)
        if ~isstruct(rawEvent) || ~isscalar(rawEvent)
            return
        end
        emitEvent(makeEvent("case",stringField(rawEvent,"caseId",""), ...
            numericField(rawEvent,"speedIndex",NaN), ...
            numericField(rawEvent,"speed_mps",NaN), ...
            numericField(rawEvent,"completedCases",completedCases), ...
            numericField(rawEvent,"totalCases",totalCases), ...
            stringField(rawEvent,"message","")));
        if isCancelledRequested()
            signalWorkerCancellation();
        end
    end

    function cleanupWorkerCancellationFiles()
        for k = 1:numel(workerCancelFiles)
            fileName = workerCancelFiles(k);
            if strlength(fileName) > 0 && isfile(char(fileName))
                try
                    delete(char(fileName));
                catch
                end
            end
        end
    end
    function signalWorkerCancellation()
        for k = 1:numel(workerCancelFiles)
            fileName = workerCancelFiles(k);
            if strlength(fileName) == 0
                continue
            end
            try
                fid = fopen(char(fileName),"w");
                if fid >= 0
                    fclose(fid);
                end
            catch
            end
        end
    end
    function emitEvent(event)
        event = normalizeEvent(event);
        events(end+1,1) = event;
        if ~isempty(callbacks.onProgress)
            try
                callbacks.onProgress(event);
            catch ME
                addStudyWarning("Progress callback failed: " + string(ME.message));
            end
        end
        if ~isempty(request.progressQueue)
            sendProgress(request.progressQueue,event);
        end
    end

    function saveCheckpoint()
        if strlength(request.checkpointPath) == 0
            return
        end
        try
            rampSpeed.saveStudy(request.checkpointPath,study);
        catch ME
            addStudyWarning("Checkpoint save failed: " + string(ME.message));
        end
    end

    function addStudyWarning(message)
        message = string(message);
        if ~any(string(study.runMeta.warnings) == message)
            study.runMeta.warnings(end+1,1) = message;
        end
    end

    function applyExecutionMetadata()
        for j = 1:numel(study.runs)
            study.runs(j).runMeta.effectiveWorkers = effectiveWorkers;
            study.runs(j).runMeta.parallelFallbackReason = string(fallbackReason);
            study.runs(j).runMeta.parallelRequested = request.parallelRequested;
            study.runs(j).runMeta.requestedWorkers = request.numWorkers;
            study.runs(j).runMeta.checkpointPath = request.checkpointPath;
            if strlength(string(fallbackReason)) > 0 && ...
                    ~any(string(study.runs(j).runMeta.warnings) == ...
                    string(fallbackReason))
                study.runs(j).runMeta.warnings(end+1,1) = ...
                    string(fallbackReason);
            end
        end
    end

    function cancelRemainingCases(indices)
        indices = double(indices(:)).';
        for j = indices
            if j < 1 || j > totalCases || ...
                    string(study.runs(j).status) == "cancelled"
                continue
            end
            run = cancelledRun(cases(j),type,request.settings, ...
                "study cancelled before case was started");
            study.runs(j) = finalizeRun(run,j,study.runs(j),type, ...
                request.settings,cases(j),effectiveWorkers,fallbackReason, ...
                started,request);
        end
    end

    function indices = findPendingCases()
        indices = [];
        for j = 1:totalCases
            if any(string(study.runs(j).status) == ["pending","running"])
                indices(end+1) = j; %#ok<AGROW>
            end
        end
    end

    function value = isCancelledRequested()
        value = false;
        if ~isempty(callbacks.isCancelled)
            value = value || readCancellation(callbacks.isCancelled);
        end
        if isfield(request,'isCancelled') && ~isempty(request.isCancelled)
            value = value || readCancellation(request.isCancelled);
        end
        if isfield(request,'cancelToken') && ~isempty(request.cancelToken)
            value = value || readCancellation(request.cancelToken);
        end
        if isfield(callbacks,'cancelToken') && ~isempty(callbacks.cancelToken)
            value = value || readCancellation(callbacks.cancelToken);
        end
    end

    function value = caseId(index)
        value = string(cases(index).id);
    end

    function value = caseLabel(index)
        value = stringField(cases(index),'label',caseId(index));
    end
end

function run = runCaseWorker(car,caseInfo,request,queue)
cancelFile = "";
if isfield(request,'cancelFile') && ~isempty(request.cancelFile)
    cancelFile = string(request.cancelFile);
end
workerCancelled = false;
callbacks = struct("onProgress",[],"isCancelled",@readWorkerCancellation);
if ~isempty(queue)
    callbacks.onProgress = @(event)sendProgress(queue, ...
        workerProgressEvent(event,caseInfo,request));
end
try
    run = feval(request.runCaseFcn,car,caseInfo,request,callbacks);
catch ME
    type = normalizeType(stringField(request,"rampType","lateral"));
    if isCancellationError(ME)
        run = cancelledRun(caseInfo,type,request.settings, ...
            "worker case cancelled: " + string(ME.message));
    else
        run = failureRun(caseInfo,ME,type,request.settings);
    end
end

    function value = readWorkerCancellation()
        if ~workerCancelled
            workerCancelled = strlength(cancelFile) > 0 && isfile(char(cancelFile));
        end
        value = workerCancelled;
    end
end


function event = workerProgressEvent(rawEvent,caseInfo,request)
event = makeEvent("case",stringField(caseInfo,"id",""), ...
    numericField(rawEvent,"speedIndex",NaN), ...
    numericField(rawEvent,"speed_mps",NaN), ...
    numericField(rawEvent,"completedCases",0), ...
    numericField(request,"totalCases",0), ...
    stringField(rawEvent,"message",""));
end

function preflightParallelInputs(pool,runCaseFcn,cars,cases,requests,callbacks,queue)
futures = parallel.FevalFuture.empty(0,numel(cases));
submitted = 0;
active = futures;
try
    for index = 1:numel(cases)
        futures(index) = parfeval(pool,@parallelSerializationProbe,1, ...
            runCaseFcn,cars{index},cases(index),requests{index}, ...
            callbacks,queue);
        submitted = index;
    end
    active = futures(1:submitted);
    while ~isempty(active)
        [completedIndex,~] = fetchNext(active);
        active(completedIndex) = [];
    end
catch ME
    if submitted > 0
        active = futures(1:submitted);
    end
    cancelAndWait(active);
    error("rampSpeed:parallelSerialization", ...
        "Parallel inputs are not serializable; falling back to serial execution: %s", ...
        string(ME.message));
end
end

function cancelAndWait(futures)
if isempty(futures)
    return
end
try
    cancel(futures);
catch
end
for k = 1:numel(futures)
    try
        wait(futures(k));
    catch
    end
end
end

function value = parallelSerializationProbe(~,~,~,~,~,queue)
if isa(queue,'parallel.pool.PollableDataQueue')
    poll(queue,0);
end
value = true;
end

function noopProgress(~)
end

function value = alwaysFalse()
value = false;
end

function validateInputs(cases,request,callbacks)
if ~isstruct(cases)
    error("rampSpeed:invalidCases","cases must be a struct array.");
end
if ~isstruct(request) || ~isscalar(request)
    error("rampSpeed:invalidRequest","request must be a scalar struct.");
end
if ~isstruct(callbacks) || ~isscalar(callbacks)
    error("rampSpeed:invalidCallbacks", ...
        "callbacks must be a scalar struct.");
end
end

function request = normalizeRequest(request)
if ~isfield(request,'rampType') || isempty(request.rampType)
    request.rampType = "lateral";
end
request.rampType = normalizeType(request.rampType);
if ~isfield(request,'settings') || isempty(request.settings)
    request.settings = struct();
end
if ~isstruct(request.settings) || ~isscalar(request.settings)
    error("rampSpeed:invalidSettings", ...
        "request.settings must be a scalar struct.");
end
if request.rampType == "lateral" && ...
        (~isfield(request.settings,'mode') || isempty(request.settings.mode) || ...
        strlength(string(request.settings.mode)) == 0)
    request.settings.mode = "coast";
end
if ~isfield(request,'parallelRequested') || isempty(request.parallelRequested)
    request.parallelRequested = false;
end
if ~isscalar(request.parallelRequested)
    error("rampSpeed:invalidParallelRequest", ...
        "parallelRequested must be scalar logical.");
end
request.parallelRequested = logical(request.parallelRequested);
if ~isfield(request,'numWorkers') || isempty(request.numWorkers)
    request.numWorkers = 0;
end
if ~isscalar(request.numWorkers) || ~isnumeric(request.numWorkers) || ...
        ~isfinite(request.numWorkers) || request.numWorkers < 0
    error("rampSpeed:invalidWorkerCount", ...
        "numWorkers must be a nonnegative finite scalar.");
end
request.numWorkers = floor(double(request.numWorkers));
if ~isfield(request,'checkpointPath') || isempty(request.checkpointPath)
    request.checkpointPath = "";
end
request.checkpointPath = string(request.checkpointPath);
if ~isscalar(request.checkpointPath)
    error("rampSpeed:invalidCheckpointPath", ...
        "checkpointPath must be scalar text.");
end
if ~isfield(request,'appVersion') || isempty(request.appVersion)
    request.appVersion = "dev";
else
    request.appVersion = string(request.appVersion);
end
if ~isfield(request,'runCaseFcn') || isempty(request.runCaseFcn)
    request.runCaseFcn = @rampSpeed.runCase;
elseif ~isa(request.runCaseFcn,'function_handle')
    error("rampSpeed:invalidRunCaseFcn", ...
        "runCaseFcn must be a function handle.");
end
if ~isfield(request,'progressQueue')
    request.progressQueue = [];
end
end

function callbacks = normalizeCallbacks(callbacks)
if ~isfield(callbacks,'onProgress')
    callbacks.onProgress = [];
end
if ~isfield(callbacks,'isCancelled')
    callbacks.isCancelled = [];
end
if ~isempty(callbacks.onProgress) && ...
        ~isa(callbacks.onProgress,'function_handle')
    error("rampSpeed:invalidCallback", ...
        "callbacks.onProgress must be a function handle.");
end
if ~isempty(callbacks.isCancelled) && ...
        ~isa(callbacks.isCancelled,'function_handle')
    error("rampSpeed:invalidCallback", ...
        "callbacks.isCancelled must be a function handle.");
end
end

function cases = normalizeCases(cases)
if isempty(cases)
    cases = repmat(emptyCase(),0,1);
    return
end
ids = strings(numel(cases),1);
for i = 1:numel(cases)
    if ~isfield(cases,'id') || strlength(string(cases(i).id)) == 0
        cases(i).id = "car-" + compose("%03d",i);
    else
        cases(i).id = string(cases(i).id);
    end
    ids(i) = cases(i).id;
    if ~isfield(cases,'label') || strlength(string(cases(i).label)) == 0
        cases(i).label = cases(i).id;
    else
        cases(i).label = string(cases(i).label);
    end
    if ~isfield(cases,'source')
        cases(i).source = "unknown";
    end
    if ~isfield(cases,'designRow') || isempty(cases(i).designRow)
        cases(i).designRow = i;
    end
    if ~isfield(cases,'sourceIndex') || isempty(cases(i).sourceIndex)
        cases(i).sourceIndex = cases(i).designRow;
    end
    if ~isfield(cases,'carRole') || isempty(cases(i).carRole)
        cases(i).carRole = "auto";
    else
        cases(i).carRole = normalizeRole(cases(i).carRole);
    end
    if ~isfield(cases,'carColumn') || isempty(cases(i).carColumn)
        cases(i).carColumn = roleColumn(cases(i).carRole);
    end
end
if numel(unique(ids)) ~= numel(ids)
    error("rampSpeed:duplicateCaseId", ...
        "cases must have unique IDs before execution.");
end
end

function [workers,reason,useParallel,pool] = chooseExecution(request)
workers = 0;
reason = "";
useParallel = false;
pool = [];
if ~request.parallelRequested
    return
end
if request.numWorkers < 2
    reason = "Parallel execution requested with fewer than two workers; falling back to serial execution.";
    return
end
try
    hasToolbox = license('test','Distrib_Computing_Toolbox') && ...
        (exist('parfeval','file') == 2 || exist('parfeval','builtin') == 5);
catch
    hasToolbox = false;
end
if ~hasToolbox
    reason = "Parallel Computing Toolbox unavailable; falling back to serial execution.";
    return
end
try
    pool = gcp('nocreate');
    if isempty(pool)
        pool = parpool('local',request.numWorkers);
    end
    workers = min(request.numWorkers,pool.NumWorkers);
    if workers < 1
        reason = "No parallel workers are available; falling back to serial execution.";
        workers = 0;
        pool = [];
        return
    end
    useParallel = true;
catch ME
    reason = "Parallel worker pool unavailable; falling back to serial execution: " + ...
        string(ME.message);
    workers = 0;
    pool = [];
end
end

function runs = initializeRuns(cases,type,settings)
if isempty(cases)
    runs = rampSpeed.makeStudy().runs;
    return
end
mode = "";
if type == "lateral"
    mode = "coast";
    if isfield(settings,'mode') && ~isempty(settings.mode)
        mode = string(settings.mode);
    end
end
runs = repmat(rampSpeed.makeRun(type,mode,settings,cases(1)), ...
    numel(cases),1);
for i = 1:numel(cases)
    runs(i) = rampSpeed.makeRun(type,mode,settings,cases(i));
end
end

function meta = makeStudyMeta(started,type,request)
meta = struct('source',"rampSpeed.runStudy",'rampType',type, ...
    'started',started,'completed',datetime.empty,'status',"running", ...
    'effectiveWorkers',0,'parallelFallbackReason',"", ...
    'warnings',strings(0,1),'errors',strings(0,1), ...
    'checkpointPath',request.checkpointPath);
end

function run = finalizeRun(run,index,template,type,settings,caseInfo, ...
        effectiveWorkers,fallbackReason,started,request)
if ~isstruct(run) || ~isscalar(run)
    invalidRunError = MException('rampSpeed:invalidRun', ...
        'runCaseFcn must return a scalar run struct.');
    run = failureRun(caseInfo,invalidRunError,type,settings);
end
if ~hasCanonicalRunFields(run)
    try
        run = rampSpeed.normalizeRampResult(run,type,settings,caseInfo, ...
            struct('source',"rampSpeed.runStudy"));
    catch ME
        run = failureRun(caseInfo,ME,type,settings);
    end
end
run = matchRunFields(run,template);
run.caseId = string(caseInfo.id);
run.type = type;
if type == "longitudinal"
    run.mode = "";
end
if ~isfield(run,'status') || isempty(run.status)
    run.status = "complete";
end
status = lower(string(run.status));
if any(status == ["completed","complete","pending","running"])
    status = "complete";
end
run.status = status;
if ~isfield(run,'runMeta') || ~isstruct(run.runMeta) || ...
        ~isscalar(run.runMeta)
    run.runMeta = struct();
end
run.runMeta.caseInfo = caseInfo;
run.runMeta.executionOrder = index;
run.runMeta.effectiveWorkers = effectiveWorkers;
run.runMeta.parallelFallbackReason = string(fallbackReason);
run.runMeta.parallelRequested = request.parallelRequested;
run.runMeta.requestedWorkers = request.numWorkers;
run.runMeta.checkpointPath = request.checkpointPath;
if ~isfield(run.runMeta,'started') || isempty(run.runMeta.started)
    run.runMeta.started = started;
end
if ~isfield(run.runMeta,'completed') || isempty(run.runMeta.completed)
    run.runMeta.completed = datetime('now');
end
run.runMeta.status = status;
if ~isfield(run.runMeta,'warnings') || isempty(run.runMeta.warnings)
    run.runMeta.warnings = strings(0,1);
else
    run.runMeta.warnings = string(run.runMeta.warnings(:));
end
if strlength(string(fallbackReason)) > 0 && ...
        ~any(run.runMeta.warnings == string(fallbackReason))
    run.runMeta.warnings(end+1,1) = string(fallbackReason);
end
if ~isfield(run.runMeta,'errors') || isempty(run.runMeta.errors)
    run.runMeta.errors = strings(0,1);
else
    run.runMeta.errors = string(run.runMeta.errors(:));
end
validationStudy = rampSpeed.makeStudy(request.appVersion);
validationStudy.cases = caseInfo;
validationStudy.runs = run;
[isValid,issues] = rampSpeed.validateStudy(validationStudy);
if ~isValid
    try
        error("rampSpeed:invalidRun", ...
            "Invalid canonical run output: %s",strjoin(issues,"; "));
    catch ME
        run = failureRun(caseInfo,ME,type,settings);
    end
    run = finalizeRun(run,index,template,type,settings,caseInfo, ...
        effectiveWorkers,fallbackReason,started,request);
end
end

function tf = hasCanonicalRunFields(run)
required = ["schemaVersion","caseId","type","mode","settings", ...
    "perSpeed","points","runMeta","status","raw"];
tf = all(isfield(run,cellstr(required)));
end

function run = matchRunFields(run,template)
names = fieldnames(template);
for i = 1:numel(names)
    if ~isfield(run,names{i})
        run.(names{i}) = template.(names{i});
    end
end
extra = setdiff(fieldnames(run),names);
if ~isempty(extra)
    run = rmfield(run,extra);
end
run = orderfields(run,names);
end

function run = failureRun(caseInfo,ME,type,settings)
mode = "";
if type == "lateral"
    mode = "coast";
    if isfield(settings,'mode') && ~isempty(settings.mode)
        mode = string(settings.mode);
    end
end
run = rampSpeed.makeRun(type,mode,settings,caseInfo);
message = string(ME.message);
identifier = string(ME.identifier);
run.status = "failed";
if height(run.perSpeed) > 0
    run.perSpeed.valid(:) = false;
    run.perSpeed.status(:) = "failed";
    run.perSpeed.reason(:) = "case failed: " + message;
end
run.runMeta.status = "failed";
run.runMeta.errors = message;
run.runMeta.error = struct("identifier",identifier, ...
    "message",message,"stack",ME.stack);
run.runMeta.errorIdentifier = identifier;
run.runMeta.errorStack = ME.stack;
run.runMeta.completed = datetime('now');
run.raw = struct('error',struct('identifier',identifier, ...
    'message',message,'stack',ME.stack),'status',"failed");
end

function run = cancelledRun(caseInfo,type,settings,reason)
mode = "";
if type == "lateral"
    mode = "coast";
    if isfield(settings,'mode') && ~isempty(settings.mode)
        mode = string(settings.mode);
    end
end
run = rampSpeed.makeRun(type,mode,settings,caseInfo);
run.status = "cancelled";
if height(run.perSpeed) > 0
    run.perSpeed.valid(:) = false;
    run.perSpeed.status(:) = "cancelled";
    run.perSpeed.reason(:) = string(reason);
end
run.runMeta.status = "cancelled";
run.runMeta.cancellationReason = string(reason);
run.runMeta.completed = datetime('now');
end

function event = makeEvent(phase,caseId,speedIndex,speed,completed,total,message)
event = struct('phase',string(phase),'caseId',string(caseId), ...
    'speedIndex',speedIndex,'speed_mps',speed, ...
    'completedCases',completed,'totalCases',total, ...
    'message',string(message));
end

function event = normalizeEvent(event)
event = makeEvent(stringField(event,'phase','case'), ...
    stringField(event,'caseId',''),numericField(event,'speedIndex',NaN), ...
    numericField(event,'speed_mps',NaN),numericField(event,'completedCases',0), ...
    numericField(event,'totalCases',0),stringField(event,'message',''));
end

function event = emptyEvent()
event = makeEvent("case","",NaN,NaN,0,0,"");
end

function value = selectedRole(caseInfo,type,request)
value = "auto";
if isfield(request,'carRole') && ~isempty(request.carRole)
    value = normalizeRole(request.carRole);
elseif isfield(caseInfo,'carRole') && ~isempty(caseInfo.carRole)
    value = normalizeRole(caseInfo.carRole);
end
if value == "auto"
    if type == "longitudinal"
        value = "acceleration";
    else
        value = "lap";
    end
end
end

function value = normalizeRole(value)
value = lower(strtrim(string(value)));
if ~isscalar(value)
    error("rampSpeed:invalidCarRole", ...
        "carRole must be scalar text.");
end
switch value
    case {"","auto","default"}
        value = "auto";
    case {"lap","lateral"}
        value = "lap";
    case {"accel","acceleration","longitudinal","acceleration-car"}
        value = "acceleration";
    otherwise
        error("rampSpeed:invalidCarRole", ...
            "carRole must be auto, lap, or acceleration.");
end
end

function value = roleColumn(role)
switch normalizeRole(role)
    case "lap"
        value = 1;
    case "acceleration"
        value = 2;
    otherwise
        value = 0;
end
end

function value = normalizeType(value)
value = lower(strtrim(string(value)));
if ~isscalar(value)
    error("rampSpeed:unsupportedType", ...
        "rampType must be scalar text.");
end
switch value
    case {"lateral","lateral-limit","lateral limit"}
        value = "lateral";
    case {"longitudinal","pure-longitudinal", ...
            "pure longitudinal","pure_longitudinal"}
        value = "longitudinal";
    otherwise
        error("rampSpeed:unsupportedType", ...
            "Unsupported ramp type: %s",value);
end
end

function value = numericField(S,name,default)
value = default;
if isstruct(S) && isfield(S,name) && ~isempty(S.(name))
    candidate = double(S.(name));
    if isscalar(candidate)
        value = candidate;
    end
end
end

function value = stringField(S,name,default)
value = string(default);
if isstruct(S) && isfield(S,name) && ~isempty(S.(name))
    candidate = string(S.(name));
    if isscalar(candidate)
        value = candidate;
    end
end
end

function value = readCancellation(token)
if isa(token,'function_handle')
    value = token();
elseif islogical(token) || isnumeric(token)
    value = token;
elseif isobject(token) && isprop(token,'Cancelled')
    value = token.Cancelled;
elseif isobject(token) && ismethod(token,'isCancelled')
    value = token.isCancelled();
else
    error("rampSpeed:invalidCancellation", ...
        "Cancellation token must be callable or scalar logical.");
end
if ~isscalar(value)
    error("rampSpeed:invalidCancellation", ...
        "Cancellation token must return a scalar logical.");
end
value = logical(value);
end

function value = isCancellationError(ME)
value = contains(lower(string(ME.identifier)),"cancel") || ...
    contains(lower(string(ME.message)),"cancel");
end

function sendProgress(queue,event)
try
    if isa(queue,'function_handle')
        queue(event);
    else
        send(queue,event);
    end
catch
    % Progress is advisory; a disconnected UI queue must not abort a solve.
end
end

function value = emptyCase()
value = struct('id',"",'label',"",'source',"", ...
    'designRow',NaN,'sourceIndex',NaN,'carRole',"auto",'carColumn',0);
end
