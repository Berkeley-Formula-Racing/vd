classdef StudyExecutor
    %STUDYEXECUTOR Own serial/parallel study lifecycle and cancellation.

    methods (Static)
        function job = start(cars,cases,request,callbacks)
            if nargin < 3 || isempty(request)
                request = struct();
            end
            if nargin < 4 || isempty(callbacks)
                callbacks = struct();
            end
            if ~isstruct(request) || ~isscalar(request)
                error("rampSpeed:invalidRequest", ...
                    "request must be a scalar struct.");
            end
            if ~isstruct(callbacks) || ~isscalar(callbacks)
                error("rampSpeed:invalidCallbacks", ...
                    "callbacks must be a scalar struct.");
            end

            job = rampSpeed.StudyJob();
            job.state = "running";
            if isfield(callbacks,"onJobCreated") && ...
                    ~isempty(callbacks.onJobCreated)
                callbacks.onJobCreated(job);
            end

            [useOuterFuture,pool] = canUseOuterFuture(request);
            if useOuterFuture
                rampSpeed.StudyExecutor.startParallel( ...
                    job,cars,cases,request,callbacks,pool);
                return
            end

            coreRequest = request;
            coreCallbacks = makeCoreCallbacks(job,callbacks);
            try
                [study,events] = rampSpeed.runStudy( ...
                    cars,cases,coreRequest,coreCallbacks);
                job.study = study;
                job.progressEvents = events;
                job.state = mapStudyState(study.status);
                job.completed = datetime("now");
            catch ME
                job.error = ME;
                job.state = "failed";
                job.completed = datetime("now");
            end
        end

        function cancel(job)
            validateJob(job);
            if ~any(string(job.state) == ["running","queued"])
                return
            end
            job.requestCancel();
            if isParallelScheduler(job)
                scheduler = job.scheduler;
                if isnan(scheduler.cancelRequestedAt)
                    scheduler.cancelRequestedAt = schedulerClock();
                    job.scheduler = scheduler;
                end
            end
        end

        function job = poll(job)
            validateJob(job);
            if isParallelScheduler(job)
                rampSpeed.StudyExecutor.pollParallel(job);
                return
            end
            if isempty(job.future) || ...
                    ~any(string(job.state) == ["running","queued"])
                return
            end
            try
                futureState = string(job.future.State);
            catch ME
                job.error = ME;
                job.state = "failed";
                job.completed = datetime("now");
                return
            end

            switch futureState
                case "finished"
                    try
                        [job.study,job.progressEvents] = fetchOutputs(job.future);
                        job.state = mapStudyState(job.study.status);
                    catch ME
                        job.error = ME;
                        job.state = "failed";
                    end
                    job.completed = datetime("now");
                case {"failed","cancelled","unavailable"}
                    try
                        job.error = job.future.Error;
                    catch
                        job.error = MException( ...
                            "rampSpeed:studyFutureFailed", ...
                            "The study worker ended in state %s.",futureState);
                    end
                    if futureState == "cancelled"
                        job.state = "cancelled";
                    else
                        job.state = "failed";
                    end
                    job.completed = datetime("now");
            end
        end

        function job = wait(job,timeoutSeconds)
            validateJob(job);
            if nargin < 2 || isempty(timeoutSeconds)
                timeoutSeconds = Inf;
            end
            if isempty(job.future) && ~isParallelScheduler(job)
                return
            end
            deadline = tic;
            while any(string(job.state) == ["running","queued"])
                rampSpeed.StudyExecutor.poll(job);
                if ~any(string(job.state) == ["running","queued"])
                    break
                end
                if toc(deadline) >= timeoutSeconds
                    break
                end
                pause(0.05);
            end
            rampSpeed.StudyExecutor.poll(job);
        end

        function [study,events] = runWorker(cars,cases,request,queue)
            request.parallelRequested = false;
            request.numWorkers = 0;
            request.progressQueue = queue;
            [study,events] = rampSpeed.runStudy(cars,cases,request,struct());
        end
    end

    methods (Static, Access=private)
        function startParallel(job,cars,cases,request,callbacks,pool)
            try
                queue = parallel.pool.DataQueue;
                afterEach(queue,@(event) ...
                    rampSpeed.StudyExecutor.forwardProgress(job,callbacks,event));
                startScheduler(job,cars,cases,request,callbacks,pool,queue);
            catch ME
                job.error = ME;
                job.state = "failed";
                job.completed = datetime("now");
            end
        end

        function pollParallel(job)
            pollScheduler(job);
        end

        function forwardProgress(job,callbacks,event)
            forwardSchedulerProgress(job,callbacks,event);
        end
    end
end

function startScheduler(job,cars,cases,request,callbacks,pool,queue)
cancelFile = string(tempname) + ".cancel";
job.cancelFile = cancelFile;

shellRequest = request;
shellRequest.parallelRequested = false;
shellRequest.numWorkers = 0;
shellRequest.progressQueue = [];
shellRequest.cancelFile = "";
shellRequest.checkpointPath = "";
shellRequest.isCancelled = [];
[study,~] = rampSpeed.runStudy(cars,cases,shellRequest, ...
    struct("isCancelled",@alwaysTrueSchedulerCancellation));

totalCases = numel(study.cases);
requestedWorkers = requestedWorkerCount(request);
poolWorkers = double(pool.NumWorkers);
if totalCases == 0
    effectiveWorkers = 0;
else
    effectiveWorkers = min([requestedWorkers,poolWorkers,totalCases]);
end
if effectiveWorkers < 1 && totalCases > 0
    error("rampSpeed:parallelPoolUnavailable", ...
        "The parallel pool has no usable workers.");
end

study = resetSchedulerStudy(study,job,request,effectiveWorkers);
s = struct();
s.parallel = true;
s.pool = pool;
s.queue = queue;
s.cars = cars;
s.cases = study.cases;
s.request = request;
s.callbacks = callbacks;
s.totalCases = totalCases;
s.effectiveWorkers = effectiveWorkers;
s.requestedWorkers = requestedWorkers;
s.nextIndex = 1;
s.completedCases = 0;
s.activeFutures = [];
s.activeIndices = zeros(0,1);
s.submitted = false(totalCases,1);
s.settled = false(totalCases,1);
s.cancelRequestedAt = NaN;
s.cancelGraceSeconds = 2;
s.futuresCancelled = false;
s.finished = false;
s.started = job.started;
s.type = schedulerType(request,study);
s.settings = schedulerSettings(request,study);

job.scheduler = s;
job.study = study;
job.future = [];
if totalCases == 0
    finishScheduler(job);
else
    pollScheduler(job);
end
end

function pollScheduler(job)
if ~isParallelScheduler(job) || job.scheduler.finished
    return
end

if schedulerCancellationRequested(job)
    s = job.scheduler;
    if isnan(s.cancelRequestedAt)
        s.cancelRequestedAt = schedulerClock();
        job.scheduler = s;
    end
end

s = job.scheduler;
if ~isnan(s.cancelRequestedAt) && ...
        ~s.futuresCancelled && ~isempty(s.activeFutures) && ...
        schedulerClock() - s.cancelRequestedAt >= s.cancelGraceSeconds
    try
        cancel(s.activeFutures);
    catch
    end
    s.futuresCancelled = true;
    job.scheduler = s;
end

indexInActive = 1;
while true
    s = job.scheduler;
    if indexInActive > numel(s.activeFutures)
        break
    end
    future = s.activeFutures(indexInActive);
    index = s.activeIndices(indexInActive);
    workerError = [];
    try
        futureState = string(future.State);
    catch ME
        futureState = "failed";
        workerError = ME;
    end
    if ~any(futureState == ["finished","failed","cancelled","unavailable"])
        indexInActive = indexInActive + 1;
        continue
    end

    if indexInActive == 1
        s.activeFutures = s.activeFutures(2:end);
        s.activeIndices = s.activeIndices(2:end);
    elseif indexInActive == numel(s.activeFutures)
        s.activeFutures = s.activeFutures(1:end-1);
        s.activeIndices = s.activeIndices(1:end-1);
    else
        s.activeFutures = [s.activeFutures(1:indexInActive-1), ...
            s.activeFutures(indexInActive+1:end)];
        s.activeIndices = [s.activeIndices(1:indexInActive-1); ...
            s.activeIndices(indexInActive+1:end)];
    end
    job.scheduler = s;
    job.future = s.activeFutures;

    run = [];
    if futureState == "finished"
        try
            workerStudy = fetchOutputs(future);
            run = schedulerRunFromStudy(workerStudy,job,index);
        catch ME
            workerError = ME;
        end
    else
        try
            workerError = future.Error;
        catch
            workerError = MException("rampSpeed:workerFutureFailed", ...
                sprintf("Setup worker ended in state %s.",futureState));
        end
    end

    if isempty(run)
        if schedulerCancellationRequested(job) || futureState == "cancelled"
            run = schedulerCancelledRun(job,index, ...
                "study cancelled while setup worker was running");
        else
            run = schedulerFailureRun(job,index,workerError);
        end
    end
    finalizeSchedulerCase(job,index,run);
    indexInActive = 1;
end

s = job.scheduler;
if ~isnan(s.cancelRequestedAt)
    if isempty(s.activeFutures)
        cancelSchedulerPending(job);
        finishScheduler(job);
    end
    return
end

scheduleSchedulerCases(job);
s = job.scheduler;
if isempty(s.activeFutures) && s.nextIndex > s.totalCases && ...
        all(s.settled)
    finishScheduler(job);
end
end

function scheduleSchedulerCases(job)
while true
    s = job.scheduler;
    if s.nextIndex > s.totalCases || ...
            numel(s.activeFutures) >= s.effectiveWorkers || ...
            schedulerCancellationRequested(job)
        return
    end
    index = s.nextIndex;
    s.nextIndex = s.nextIndex + 1;
    s.submitted(index) = true;
    job.scheduler = s;

    try
        [car,sourceIndex,column] = selectSchedulerCar( ...
            s.cars,s.cases(index),index,s.request);
        workerCars = singleSchedulerCar(car,sourceIndex,column);
        workerRequest = schedulerWorkerRequest(s.request,job.cancelFile);
        future = parfeval(s.pool,@runCaseWorker,1,workerCars, ...
            s.cases(index),workerRequest,s.queue);
        s = job.scheduler;
        if isempty(s.activeFutures)
            s.activeFutures = future;
        else
            s.activeFutures(end+1) = future;
        end
        if isempty(s.activeIndices)
            s.activeIndices = index;
        else
            s.activeIndices = [s.activeIndices(:); index];
        end
        job.scheduler = s;
        job.future = s.activeFutures;
        emitSchedulerEvent(job,schedulerEvent("case",s.cases(index), ...
            s.completedCases,s.totalCases, ...
            "starting " + schedulerCaseLabel(s.cases(index))));
    catch ME
        run = schedulerFailureRun(job,index,ME);
        finalizeSchedulerCase(job,index,run);
    end
end
end

function finalizeSchedulerCase(job,index,run)
s = job.scheduler;
if index < 1 || index > s.totalCases || s.settled(index)
    return
end
run = patchSchedulerRun(run,job,index);
study = job.study;
try
    study.runs(index) = run;
catch ME
    run = schedulerFailureRun(job,index,ME);
    study.runs(index) = run;
end
s.settled(index) = true;
s.completedCases = s.completedCases + 1;
study = appendSchedulerRunMetadata(study,run,index);
job.study = study;
job.scheduler = s;
emitSchedulerEvent(job,schedulerEvent("case",s.cases(index), ...
    s.completedCases,s.totalCases, ...
    string(run.status) + ": " + schedulerCaseLabel(s.cases(index))));
saveSchedulerCheckpoint(job);
end

function cancelSchedulerPending(job)
s = job.scheduler;
for index = 1:s.totalCases
    if s.settled(index) || s.submitted(index)
        continue
    end
    run = schedulerCancelledRun(job,index, ...
        "study cancelled before setup was started");
    finalizeSchedulerCase(job,index,run);
end
end

function finishScheduler(job)
s = job.scheduler;
if s.finished || ~isempty(s.activeFutures)
    return
end
if ~isnan(s.cancelRequestedAt)
    cancelSchedulerPending(job);
    studyStatus = "cancelled";
else
    studyStatus = schedulerStudyStatus(job.study);
end
study = job.study;
study.status = studyStatus;
if ~isfield(study,"runMeta") || ~isstruct(study.runMeta)
    study.runMeta = struct();
end
study.runMeta.completed = datetime("now");
study.runMeta.status = studyStatus;
job.study = study;
saveSchedulerCheckpoint(job);
s = job.scheduler;
s.finished = true;
job.scheduler = s;
job.future = [];
job.state = mapStudyState(studyStatus);
job.completed = datetime("now");
end

function value = schedulerCancellationRequested(job)
value = job.isCancelled();
s = job.scheduler;
if value
    return
end
tokens = {s.request,s.callbacks};
for k = 1:numel(tokens)
    tokenRecord = tokens{k};
    for name = ["isCancelled","cancelToken"]
        if isstruct(tokenRecord) && isfield(tokenRecord,char(name)) && ...
                ~isempty(tokenRecord.(char(name)))
            try
                value = value || readSchedulerCancellation( ...
                    tokenRecord.(char(name)));
            catch
            end
        end
        if value
            job.requestCancel();
            return
        end
    end
end
end

function value = readSchedulerCancellation(token)
if isa(token,"function_handle")
    value = token();
elseif islogical(token) || isnumeric(token)
    value = token;
elseif isobject(token) && isprop(token,"Cancelled")
    value = token.Cancelled;
elseif isobject(token) && ismethod(token,"isCancelled")
    value = token.isCancelled();
else
    value = false;
end
value = isscalar(value) && logical(value);
end

function value = alwaysTrueSchedulerCancellation()
value = true;
end

function emitSchedulerEvent(job,event)
if ~isstruct(event) || ~isscalar(event)
    return
end
appendEvent(job,event);
s = job.scheduler;
callbacks = s.callbacks;
if isfield(callbacks,"onProgress") && ~isempty(callbacks.onProgress)
    try
        callbacks.onProgress(event);
    catch
        % UI progress is advisory.
    end
end
end

function saveSchedulerCheckpoint(job)
s = job.scheduler;
if ~isfield(s.request,"checkpointPath") || ...
        isempty(s.request.checkpointPath)
    return
end
path = string(s.request.checkpointPath);
if ~isscalar(path) || strlength(strtrim(path)) == 0
    return
end
try
    rampSpeed.saveStudy(path,job.study);
catch ME
    study = job.study;
    if ~isfield(study,"runMeta") || ~isstruct(study.runMeta)
        study.runMeta = struct();
    end
    if ~isfield(study.runMeta,"warnings") || isempty(study.runMeta.warnings)
        study.runMeta.warnings = strings(0,1);
    else
        study.runMeta.warnings = string(study.runMeta.warnings(:));
    end
    message = "Checkpoint save failed: " + string(ME.message);
    if ~any(study.runMeta.warnings == message)
        study.runMeta.warnings(end+1,1) = message;
    end
    job.study = study;
end
end

function study = resetSchedulerStudy(study,job,request,effectiveWorkers)
if ~isfield(study,"runMeta") || ~isstruct(study.runMeta)
    study.runMeta = struct();
end
study.status = "running";
study.runMeta.status = "running";
study.runMeta.started = job.started;
study.runMeta.completed = datetime.empty;
study.runMeta.effectiveWorkers = effectiveWorkers;
study.runMeta.parallelFallbackReason = "";
study.runMeta.parallelRequested = logicalField(request,"parallelRequested",true);
study.runMeta.requestedWorkers = requestedWorkerCount(request);
study.runMeta.checkpointPath = stringField(request,"checkpointPath","");
if ~isfield(study.runMeta,"warnings") || isempty(study.runMeta.warnings)
    study.runMeta.warnings = strings(0,1);
else
    study.runMeta.warnings = string(study.runMeta.warnings(:));
end
if ~isfield(study.runMeta,"errors") || isempty(study.runMeta.errors)
    study.runMeta.errors = strings(0,1);
else
    study.runMeta.errors = string(study.runMeta.errors(:));
end
for index = 1:numel(study.runs)
    run = study.runs(index);
    run.status = "pending";
    run.raw = struct();
    if istable(run.perSpeed) && height(run.perSpeed) > 0
        run.perSpeed.valid(:) = false;
        run.perSpeed.status(:) = "pending";
        run.perSpeed.reason(:) = "";
    end
    if istable(run.points) && height(run.points) > 0
        run.points.valid(:) = false;
        run.points.status(:) = "pending";
    end
    if ~isfield(run,"runMeta") || ~isstruct(run.runMeta)
        run.runMeta = struct();
    end
    run.runMeta.status = "pending";
    run.runMeta.started = datetime.empty;
    run.runMeta.completed = datetime.empty;
    run.runMeta.effectiveWorkers = effectiveWorkers;
    run.runMeta.parallelFallbackReason = "";
    run.runMeta.parallelRequested = logicalField(request,"parallelRequested",true);
    run.runMeta.requestedWorkers = requestedWorkerCount(request);
    run.runMeta.checkpointPath = stringField(request,"checkpointPath","");
    if isfield(run.runMeta,"cancellationReason")
        run.runMeta = rmfield(run.runMeta,"cancellationReason");
    end
    if isfield(run.runMeta,"error")
        run.runMeta = rmfield(run.runMeta,"error");
    end
    run.runMeta.errors = strings(0,1);
    study.runs(index) = run;
end
end

function callbacks = makeCoreCallbacks(job,outerCallbacks)
callbacks = struct();
callbacks.onProgress = @(event)forwardCoreProgress(job,outerCallbacks,event);
callbacks.isCancelled = @()job.isCancelled();
end

function run = patchSchedulerRun(run,job,index)
s = job.scheduler;
caseInfo = s.cases(index);
requiredFields = {"schemaVersion","caseId","type","mode", ...
    "settings","perSpeed","points","runMeta","status","raw"};
if ~isstruct(run) || ~isscalar(run)
    run = schedulerFailureRun(job,index,MException( ...
        "rampSpeed:invalidWorkerResult", ...
        "Setup worker did not return a scalar run struct."));
    return
end
missingFields = requiredFields(~cellfun(@(name)isfield(run,name),requiredFields));
if ~isempty(missingFields)
    run = schedulerFailureRun(job,index,MException( ...
        "rampSpeed:invalidWorkerResult",sprintf( ...
        "Setup worker run is missing fields: %s.",strjoin(missingFields,", "))));
    return
end
run.caseId = string(caseInfo.id);
run.type = s.type;
if ~isfield(run,"runMeta") || ~isstruct(run.runMeta)
    run.runMeta = struct();
end
run.runMeta.caseInfo = caseInfo;
if isfield(caseInfo,"setupSpec")
    run.runMeta.setupSpec = caseInfo.setupSpec;
end
if isfield(caseInfo,"derived")
    run.runMeta.setupDerived = caseInfo.derived;
end
status = lower(string(run.status));
if any(status == ["completed","complete"])
    status = "complete";
end
run.status = status;
run.runMeta.status = status;
run.runMeta.executionOrder = index;
run.runMeta.effectiveWorkers = s.effectiveWorkers;
run.runMeta.parallelFallbackReason = "";
run.runMeta.parallelRequested = logicalField(s.request,"parallelRequested",true);
run.runMeta.requestedWorkers = s.requestedWorkers;
run.runMeta.checkpointPath = stringField(s.request,"checkpointPath","");
if ~isfield(run.runMeta,"started") || isempty(run.runMeta.started)
    run.runMeta.started = s.started;
end
if ~isfield(run.runMeta,"completed") || isempty(run.runMeta.completed)
    run.runMeta.completed = datetime("now");
end
if ~isfield(run.runMeta,"warnings") || isempty(run.runMeta.warnings)
    run.runMeta.warnings = strings(0,1);
else
    run.runMeta.warnings = string(run.runMeta.warnings(:));
end
if ~isfield(run.runMeta,"errors") || isempty(run.runMeta.errors)
    run.runMeta.errors = strings(0,1);
else
    run.runMeta.errors = string(run.runMeta.errors(:));
end
end

function study = appendSchedulerRunMetadata(study,run,index)
if ~isfield(study,"runMeta") || ~isstruct(study.runMeta)
    study.runMeta = struct();
end
if ~isfield(study.runMeta,"warnings") || isempty(study.runMeta.warnings)
    study.runMeta.warnings = strings(0,1);
else
    study.runMeta.warnings = string(study.runMeta.warnings(:));
end
if ~isfield(study.runMeta,"errors") || isempty(study.runMeta.errors)
    study.runMeta.errors = strings(0,1);
else
    study.runMeta.errors = string(study.runMeta.errors(:));
end
if isfield(run.runMeta,"warnings")
    for message = string(run.runMeta.warnings(:).')
        if strlength(strtrim(message)) > 0 && ...
                ~any(study.runMeta.warnings == message)
            study.runMeta.warnings(end+1,1) = message;
        end
    end
end
if string(run.status) == "failed"
    message = schedulerRunErrorMessage(run);
    if strlength(strtrim(message)) > 0
        caseMessage = string(run.caseId) + ": " + message;
        if ~any(study.runMeta.errors == caseMessage)
            study.runMeta.errors(end+1,1) = caseMessage;
        end
    end
end
if isfield(study.runMeta,"status") && string(study.runMeta.status) == "running"
    study.runMeta.lastCompletedCase = index;
end
end

function message = schedulerRunErrorMessage(run)
message = "";
if isfield(run.runMeta,"error") && isstruct(run.runMeta.error) && ...
        isfield(run.runMeta.error,"message")
    message = string(run.runMeta.error.message);
elseif isfield(run.runMeta,"errors") && ~isempty(run.runMeta.errors)
    message = string(run.runMeta.errors(1));
end
end

function run = schedulerRunFromStudy(workerStudy,job,index)
run = [];
if ~isstruct(workerStudy) || ~isscalar(workerStudy)
    return
end
if isfield(workerStudy,"runs")
    if isempty(workerStudy.runs)
        return
    end
    try
        run = workerStudy.runs(1);
    catch
        run = [];
    end
else
    requiredFields = {"schemaVersion","caseId","type","mode", ...
        "settings","perSpeed","points","runMeta","status","raw"};
    if all(cellfun(@(name)isfield(workerStudy,name),requiredFields))
        run = workerStudy;
    end
end
if ~isempty(run) && isfield(workerStudy,"runMeta") && ...
        isstruct(workerStudy.runMeta) && isfield(workerStudy.runMeta,"warnings")
    study = job.study;
    if ~isfield(study,"runMeta") || ~isstruct(study.runMeta)
        study.runMeta = struct();
    end
    if ~isfield(study.runMeta,"warnings") || isempty(study.runMeta.warnings)
        study.runMeta.warnings = strings(0,1);
    else
        study.runMeta.warnings = string(study.runMeta.warnings(:));
    end
    for message = string(workerStudy.runMeta.warnings(:).')
        if strlength(strtrim(message)) > 0 && ...
                ~any(study.runMeta.warnings == message)
            study.runMeta.warnings(end+1,1) = message;
        end
    end
    job.study = study;
end
if isempty(run)
    run = schedulerFailureRun(job,index,MException( ...
        "rampSpeed:invalidWorkerResult", ...
        "Setup worker returned no run."));
end
end

function run = schedulerFailureRun(job,index,ME)
s = job.scheduler;
if isempty(ME) || ~isa(ME,"MException")
    ME = MException("rampSpeed:workerFailed","Setup worker failed.");
end
caseInfo = s.cases(index);
run = rampSpeed.makeRun(s.type,schedulerMode(s),s.settings,caseInfo);
run.status = "failed";
message = string(ME.message);
identifier = string(ME.identifier);
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
run.runMeta.completed = datetime("now");
run.raw = struct("error",struct("identifier",identifier, ...
    "message",message,"stack",ME.stack),"status","failed");
end

function run = schedulerCancelledRun(job,index,reason)
s = job.scheduler;
caseInfo = s.cases(index);
run = rampSpeed.makeRun(s.type,schedulerMode(s),s.settings,caseInfo);
run.status = "cancelled";
if height(run.perSpeed) > 0
    run.perSpeed.valid(:) = false;
    run.perSpeed.status(:) = "cancelled";
    run.perSpeed.reason(:) = string(reason);
end
run.runMeta.status = "cancelled";
run.runMeta.cancellationReason = string(reason);
run.runMeta.completed = datetime("now");
end

function request = schedulerWorkerRequest(request,cancelFile)
request.parallelRequested = false;
request.numWorkers = 0;
request.progressQueue = [];
request.checkpointPath = "";
request.cancelFile = string(cancelFile);
request.totalCases = 1;
end

function value = schedulerStudyStatus(study)
if ~isfield(study,"runs") || isempty(study.runs)
    value = "complete";
    return
end
statuses = strings(numel(study.runs),1);
for index = 1:numel(study.runs)
    statuses(index) = string(study.runs(index).status);
end
if any(statuses == "cancelled")
    value = "cancelled";
elseif all(statuses == "complete")
    value = "complete";
elseif all(ismember(statuses,["complete","partial","warning","completed"]))
    if any(ismember(statuses,["partial","warning"]))
        value = "partial";
    else
        value = "complete";
    end
elseif any(ismember(statuses,["complete","partial","warning","completed"]))
    value = "partial";
else
    value = "failed";
end
end

function value = schedulerMode(s)
value = "";
if s.type == "lateral"
    value = "coast";
    if isfield(s.settings,"mode") && ~isempty(s.settings.mode)
        value = string(s.settings.mode);
    end
end
end

function value = schedulerType(request,study)
value = "lateral";
if isfield(study,"runMeta") && isfield(study.runMeta,"rampType") && ...
        ~isempty(study.runMeta.rampType)
    value = string(study.runMeta.rampType);
elseif isfield(request,"rampType") && ~isempty(request.rampType)
    value = lower(strtrim(string(request.rampType)));
end
if value == "lateral-limit" || value == "lateral limit"
    value = "lateral";
elseif any(value == ["longitudinal","pure-longitudinal","pure longitudinal", ...
        "pure_longitudinal"])
    value = "longitudinal";
end
end

function value = schedulerSettings(request,study)
value = struct();
if isfield(request,"settings") && isstruct(request.settings) && ...
        isscalar(request.settings)
    value = request.settings;
end
if isfield(study,"runs") && ~isempty(study.runs) && ...
        isfield(study.runs(1),"settings")
    value = study.runs(1).settings;
end
end

function value = requestedWorkerCount(request)
value = 0;
if isfield(request,"numWorkers") && ~isempty(request.numWorkers)
    value = floor(double(request.numWorkers));
end
end

function value = logicalField(record,name,default)
value = default;
if isstruct(record) && isfield(record,name) && ~isempty(record.(name))
    value = logical(record.(name));
end
end

function value = stringField(record,name,default)
value = string(default);
if isstruct(record) && isfield(record,name) && ~isempty(record.(name))
    candidate = string(record.(name));
    if isscalar(candidate)
        value = candidate;
    end
end
end

function value = schedulerCaseLabel(caseInfo)
value = stringField(caseInfo,"label",stringField(caseInfo,"id",""));
end

function forwardCoreProgress(job,outerCallbacks,event)
appendEvent(job,event);
if isfield(outerCallbacks,"onProgress") && ...
        ~isempty(outerCallbacks.onProgress)
    try
        outerCallbacks.onProgress(event);
    catch
        % Progress callbacks are advisory and must not abort a solve.
    end
end
end

function [car,sourceIndex,column] = selectSchedulerCar(cars,caseInfo,index,request)
type = schedulerType(request,struct());
role = "auto";
if isfield(request,"carRole") && ~isempty(request.carRole)
    role = schedulerRole(request.carRole);
elseif isfield(caseInfo,"carRole") && ~isempty(caseInfo.carRole)
    role = schedulerRole(caseInfo.carRole);
end
if role == "auto"
    if type == "longitudinal"
        role = "acceleration";
    else
        role = "lap";
    end
end
column = 1;
if role == "acceleration"
    column = 2;
end
sourceIndex = numericSchedulerField(caseInfo,"sourceIndex", ...
    numericSchedulerField(caseInfo,"designRow",index));
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
    if nRows == 1 && isvector(cars) && sourceIndex > 1
        if numel(cars) >= sourceIndex
            car = cars{sourceIndex};
            return
        end
        error("rampSpeed:carCellDesignMismatch", ...
            "Row-vector cars do not contain sourceIndex %d.",sourceIndex);
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
    "No car is available for case %s.",string(caseInfo.id));
end

function value = singleSchedulerCar(car,sourceIndex,column)
rows = max(1,sourceIndex);
columns = max(1,column);
value = cell(rows,columns);
value{sourceIndex,column} = car;
end

function value = numericSchedulerField(record,name,default)
value = default;
if isstruct(record) && isfield(record,name) && ~isempty(record.(name))
    candidate = double(record.(name));
    if isscalar(candidate)
        value = candidate;
    end
end
end

function value = schedulerRole(value)
value = lower(strtrim(string(value)));
if ~isscalar(value) || any(value == ["","auto","default"])
    if isscalar(value)
        value = "auto";
    else
        error("rampSpeed:invalidCarRole","carRole must be scalar text.");
    end
    return
end
switch value
    case {"lap","lateral"}
        value = "lap";
    case {"accel","acceleration","longitudinal","acceleration-car"}
        value = "acceleration";
    otherwise
        error("rampSpeed:invalidCarRole", ...
            "carRole must be auto, lap, or acceleration.");
end
end

function event = schedulerEvent(phase,caseInfo,completed,total,message)
event = struct("phase",string(phase), ...
    "caseId",stringField(caseInfo,"id",""), ...
    "speedIndex",NaN,"speed_mps",NaN, ...
    "completedCases",completed,"totalCases",total, ...
    "message",string(message));
end

function sendSchedulerWorkerProgress(queue,rawEvent,caseInfo)
if ~isstruct(rawEvent) || ~isscalar(rawEvent)
    return
end
event = rawEvent;
if ~isfield(event,"caseId") || strlength(string(event.caseId)) == 0
    event.caseId = stringField(caseInfo,"id","");
end
try
    send(queue,event);
catch
    % A disconnected progress queue must not abort a setup.
end
end

function run = runCaseWorker(singleCars,caseInfo,request,queue)
cancelFile = stringField(request,"cancelFile","");
callbacks = struct("onProgress",[], ...
    "isCancelled",@readWorkerCancellation);
if ~isempty(queue)
    callbacks.onProgress = @(event)sendSchedulerWorkerProgress( ...
        queue,event,caseInfo);
end
[workerStudy,~] = rampSpeed.runStudy(singleCars,caseInfo,request,callbacks);
run = [];
if isstruct(workerStudy) && isscalar(workerStudy) && ...
        isfield(workerStudy,"runs") && ~isempty(workerStudy.runs)
    run = workerStudy.runs(1);
end

    function value = readWorkerCancellation()
        value = strlength(cancelFile) > 0 && isfile(char(cancelFile));
    end
end

function forwardSchedulerProgress(job,callbacks,event)
if ~isstruct(event) || ~isscalar(event)
    return
end
if isParallelScheduler(job)
    s = job.scheduler;
    phase = stringField(event,"phase","case");
    if phase == "study"
        return
    end
    message = stringField(event,"message","");
    trimmed = lower(strtrim(message));
    if phase == "case" && ...
            (startsWith(trimmed,"starting ") || ...
            startsWith(trimmed,"complete:") || ...
            startsWith(trimmed,"failed:") || ...
            startsWith(trimmed,"cancelled:") || ...
            startsWith(trimmed,"partial:") || startsWith(trimmed,"warning:"))
        return
    end
    event.phase = phase;
    event.caseId = stringField(event,"caseId","");
    event.speedIndex = numericSchedulerField(event,"speedIndex",NaN);
    event.speed_mps = numericSchedulerField(event,"speed_mps",NaN);
    event.completedCases = s.completedCases;
    event.totalCases = s.totalCases;
    event.message = message;
end
appendEvent(job,event);
if isfield(callbacks,"onProgress") && ~isempty(callbacks.onProgress)
    try
        callbacks.onProgress(event);
    catch
        % UI progress is advisory.
    end
end
end

function value = isParallelScheduler(job)
value = false;
try
    value = isstruct(job.scheduler) && isfield(job.scheduler,"parallel") && ...
        logical(job.scheduler.parallel);
catch
end
end

function value = schedulerClock()
value = posixtime(datetime("now"));
end

function appendEvent(job,event)
if isempty(job.progressEvents)
    job.progressEvents = event;
    return
end
try
    job.progressEvents(end+1,1) = event;
catch
    % Keep a valid event collection even if an injected callback changes fields.
    job.progressEvents = repmat(event,numel(job.progressEvents)+1,1);
end
end

function [tf,pool] = canUseOuterFuture(request)
tf = false;
pool = [];
if ~isfield(request,"parallelRequested") || ...
        ~logical(request.parallelRequested)
    return
end
workers = 0;
if isfield(request,"numWorkers") && ~isempty(request.numWorkers)
    workers = double(request.numWorkers);
end
if ~isscalar(workers) || workers < 2
    return
end
try
    hasToolbox = license("test","Distrib_Computing_Toolbox") && ...
        (exist("parfeval","file") == 2 || exist("parfeval","builtin") == 5);
catch
    hasToolbox = false;
end
if ~hasToolbox
    return
end
try
    pool = gcp("nocreate");
    if isempty(pool)
        pool = parpool("local",floor(workers));
    end
    tf = ~isempty(pool);
catch
    pool = [];
    tf = false;
end
end

function value = mapStudyState(status)
status = lower(string(status));
switch status
    case "cancelled"
        value = "cancelled";
    case "failed"
        value = "failed";
    case {"complete","completed","partial","warning"}
        value = "completed";
    otherwise
        value = "failed";
end
end

function validateJob(job)
if ~isa(job,"rampSpeed.StudyJob") || ~isscalar(job)
    error("rampSpeed:invalidStudyJob", ...
        "job must be a scalar rampSpeed.StudyJob.");
end
end
