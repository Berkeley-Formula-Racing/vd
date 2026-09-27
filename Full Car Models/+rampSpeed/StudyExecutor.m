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
        end

        function job = poll(job)
            validateJob(job);
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
            if isempty(job.future)
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
            queue = parallel.pool.DataQueue;
            afterEach(queue,@(event) ...
                rampSpeed.StudyExecutor.forwardProgress(job,callbacks,event));
            request.parallelRequested = false;
            request.numWorkers = 0;
            request.progressQueue = queue;
            request.cancelFile = string(tempname) + ".cancel";
            job.cancelFile = request.cancelFile;
            try
                job.future = parfeval(pool, ...
                    @rampSpeed.StudyExecutor.runWorker,2, ...
                    cars,cases,request,queue);
            catch ME
                job.error = ME;
                job.state = "failed";
                job.completed = datetime("now");
            end
        end

        function forwardProgress(job,callbacks,event)
            if ~isstruct(event) || ~isscalar(event)
                return
            end
            appendEvent(job,event);
            if isfield(callbacks,"onProgress") && ...
                    ~isempty(callbacks.onProgress)
                try
                    callbacks.onProgress(event);
                catch
                    % UI progress is advisory.
                end
            end
        end
    end
end

function callbacks = makeCoreCallbacks(job,outerCallbacks)
callbacks = struct();
callbacks.onProgress = @(event)forwardCoreProgress(job,outerCallbacks,event);
callbacks.isCancelled = @()job.isCancelled();
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
