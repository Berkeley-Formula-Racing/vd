function tests = test_studyExecutorParallelScheduler
%TEST_STUDYEXECUTORPARALLELSCHEDULER Regression tests for setup-level jobs.

tests = functiontests(localfunctions);
end

function testParallelExecutionOverlapsIndependentSetups(testCase)
assumeTrue(testCase,hasParallelRuntime(), ...
    "Parallel Computing Toolbox is required for this scheduler test.");
cleanupState = containers.Map('KeyType','char','ValueType','any');
cleanupState('poolBefore') = gcp("nocreate");
cleanupState('job') = [];
fixture = makeRampFixture();
cases = makeFourCases(fixture.cases);
logDirectory = string(tempname);
mkdir(char(logDirectory));
cleanup = onCleanup(@()removeDirectory(logDirectory)); %#ok<NASGU>
cleanupResources = onCleanup(@()cleanupSchedulerTestResources(cleanupState)); %#ok<NASGU>

request = parallelRequest(@recordTimedCase);
request.logDirectory = logDirectory;
cleanupState('job') = rampSpeed.StudyExecutor.start( ...
    fixture.cars,cases,request,struct());
job = cleanupState('job');
rampSpeed.StudyExecutor.wait(job,Inf);

verifyEqual(testCase,job.state,"completed");
verifyEqual(testCase,string({job.study.runs.caseId}),string({cases.id}));
verifyTrue(testCase,all(string({job.study.runs.status}) == "complete"));
verifyEqual(testCase,job.study.runMeta.effectiveWorkers,2);
executionOrder = arrayfun(@(run)run.runMeta.executionOrder,job.study.runs);
verifyEqual(testCase,executionOrder(:),(1:4).');

intervals = readIntervals(logDirectory,cases);
verifyEqual(testCase,size(intervals),[numel(cases),2]);
overlap = false;
for i = 1:size(intervals,1)
    for j = i+1:size(intervals,1)
        overlap = overlap || intervals(i,1) < intervals(j,2) && ...
            intervals(j,1) < intervals(i,2);
    end
end
verifyTrue(testCase,overlap, ...
    "At least two setup executions should overlap in time.");
verifyTrue(testCase,all(countLogStarts(logDirectory,cases) == 1), ...
    "Each setup must execute exactly once.");
end

function testWorkerFailureRetainsOtherRunsAndCatalogOrder(testCase)
assumeTrue(testCase,hasParallelRuntime(), ...
    "Parallel Computing Toolbox is required for this scheduler test.");
cleanupState = containers.Map('KeyType','char','ValueType','any');
cleanupState('poolBefore') = gcp("nocreate");
cleanupState('job') = [];
fixture = makeRampFixture();
cases = makeFourCases(fixture.cases);
cases(2).id = "fail";
logDirectory = string(tempname);
mkdir(char(logDirectory));
cleanup = onCleanup(@()removeDirectory(logDirectory)); %#ok<NASGU>
cleanupResources = onCleanup(@()cleanupSchedulerTestResources(cleanupState)); %#ok<NASGU>

request = parallelRequest(@recordFailureCase);
request.logDirectory = logDirectory;
cleanupState('job') = rampSpeed.StudyExecutor.start( ...
    fixture.cars,cases,request,struct());
job = cleanupState('job');
rampSpeed.StudyExecutor.wait(job,Inf);

verifyEqual(testCase,job.state,"completed");
verifyEqual(testCase,job.study.status,"partial");
verifyEqual(testCase,string({job.study.runs.caseId}),string({cases.id}));
verifyEqual(testCase,string({job.study.runs.status}), ...
    ["complete","failed","complete","complete"]);
verifyEqual(testCase,job.study.runs(2).runMeta.error.identifier, ...
    "test:parallelWorkerFailure");
verifyTrue(testCase,isempty(job.error));
verifyTrue(testCase,all(countLogStarts(logDirectory,cases) == 1), ...
    "A failed setup must not be retried or duplicated.");
end

function testCancellationCreatesTerminalRowsAndCheckpoint(testCase)
assumeTrue(testCase,hasParallelRuntime(), ...
    "Parallel Computing Toolbox is required for this scheduler test.");
cleanupState = containers.Map('KeyType','char','ValueType','any');
cleanupState('poolBefore') = gcp("nocreate");
cleanupState('job') = [];
fixture = makeRampFixture();
cases = makeFourCases(fixture.cases);
checkpointPath = string(tempname) + ".mat";
logDirectory = string(tempname);
mkdir(char(logDirectory));
cleanupCheckpoint = onCleanup(@()deleteIfPresent(checkpointPath)); %#ok<NASGU>
cleanupDirectory = onCleanup(@()removeDirectory(logDirectory)); %#ok<NASGU>
jobReference = [];
cleanupResources = onCleanup(@()cleanupSchedulerTestResources(cleanupState)); %#ok<NASGU>

request = parallelRequest(@recordCancellableCase);
request.logDirectory = logDirectory;
request.checkpointPath = checkpointPath;
callbacks = struct("onJobCreated",@rememberJob, ...
    "onProgress",@cancelAfterWorkerStarts);
cleanupState('job') = rampSpeed.StudyExecutor.start( ...
    fixture.cars,cases,request,callbacks);
job = cleanupState('job');
rampSpeed.StudyExecutor.wait(job,Inf);

verifyEqual(testCase,job.state,"cancelled");
verifyEqual(testCase,job.study.status,"cancelled");
verifyTrue(testCase,isfile(checkpointPath));
statuses = string({job.study.runs.status});
verifyTrue(testCase,all(ismember(statuses,["complete","failed","cancelled"])));
verifyFalse(testCase,any(ismember(statuses,["pending","running"])));
checkpoint = rampSpeed.loadStudy(checkpointPath,"scheduler-test");
verifyEqual(testCase,checkpoint.status,"cancelled");
verifyEqual(testCase,string({checkpoint.runs.status}),statuses);
verifyTrue(testCase,all(countLogStarts(logDirectory,cases) <= 1), ...
    "Cancellation must not cause setup retries.");

    function rememberJob(createdJob)
        jobReference = createdJob;
    end

    function cancelAfterWorkerStarts(event)
        if isempty(jobReference) || ~isstruct(event) || ...
                ~isfield(event,"message")
            return
        end
        if string(event.message) == "worker-started"
            rampSpeed.StudyExecutor.cancel(jobReference);
        end
    end
end

function request = parallelRequest(runCaseFcn)
request = struct("rampType","lateral", ...
    "settings",struct("speeds",5,"nRamp",1,"nBisect",0,"mode","coast"), ...
    "parallelRequested",true,"numWorkers",2, ...
    "checkpointPath","","appVersion","scheduler-test", ...
    "runCaseFcn",runCaseFcn);
end

function cases = makeFourCases(twoCases)
cases = repmat(twoCases(1),4,1);
cases(1) = twoCases(1);
cases(1).id = "first";
cases(1).label = "first";
cases(1).designRow = 1;
cases(1).sourceIndex = 1;
cases(2).carRole = twoCases(2).carRole;
cases(2).id = "second";
cases(2).label = "second";
cases(2).designRow = 2;
cases(2).sourceIndex = 2;
cases(3).carRole = twoCases(1).carRole;
cases(3).id = "third";
cases(3).label = "third";
cases(3).designRow = 1;
cases(3).sourceIndex = 1;
cases(4).carRole = twoCases(2).carRole;
cases(4).id = "fourth";
cases(4).label = "fourth";
cases(4).designRow = 2;
cases(4).sourceIndex = 2;
end

function run = recordTimedCase(~,caseInfo,request,~)
path = logPath(request,caseInfo);
writeLog(path,"start");
pause(0.45);
writeLog(path,"finish");
run = makeFixtureRun(caseInfo);
end

function run = recordFailureCase(~,caseInfo,request,~)
path = logPath(request,caseInfo);
writeLog(path,"start");
pause(0.1);
if string(caseInfo.id) == "fail"
    error("test:parallelWorkerFailure","injected worker failure");
end
writeLog(path,"finish");
run = makeFixtureRun(caseInfo);
end

function run = recordCancellableCase(~,caseInfo,request,callbacks)
path = logPath(request,caseInfo);
writeLog(path,"start");
if ~isempty(callbacks.onProgress)
    callbacks.onProgress(struct("phase","speed","speedIndex",1, ...
        "speed_mps",5,"message","worker-started"));
end
for i = 1:100
    if callbacks.isCancelled()
        writeLog(path,"cancelled");
        run = rampSpeed.makeRun("lateral","coast",request.settings,caseInfo);
        run.status = "cancelled";
        run.runMeta.status = "cancelled";
        run.runMeta.cancellationReason = "test cancellation";
        run.perSpeed.status(:) = "cancelled";
        run.perSpeed.reason(:) = "test cancellation";
        return
    end
    pause(0.03);
end
writeLog(path,"finish");
run = makeFixtureRun(caseInfo);
end

function path = logPath(request,caseInfo)
path = fullfile(char(request.logDirectory),char(string(caseInfo.id) + ".log"));
end

function writeLog(path,label)
fid = fopen(path,"a");
if fid < 0
    error("test:logOpen","Could not open worker log.");
end
cleanup = onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,"%s %.17g\n",label,posixtime(datetime("now")));
end

function intervals = readIntervals(directory,cases)
intervals = NaN(numel(cases),2);
for i = 1:numel(cases)
    path = logPath(struct("logDirectory",directory),cases(i));
    lines = splitlines(strtrim(string(fileread(path))));
    intervals(i,1) = sscanf(char(lines(1)),"start %f");
    intervals(i,2) = sscanf(char(lines(2)),"finish %f");
end
end

function counts = countLogStarts(directory,cases)
counts = zeros(numel(cases),1);
for i = 1:numel(cases)
    path = logPath(struct("logDirectory",directory),cases(i));
    if isfile(path)
        lines = splitlines(strtrim(string(fileread(path))));
        counts(i) = sum(startsWith(lines,"start "));
    end
end
end

function value = hasParallelRuntime()
try
    value = license("test","Distrib_Computing_Toolbox") && ...
        (exist("parfeval","file") == 2 || exist("parfeval","builtin") == 5);
catch
    value = false;
end
end

function deleteIfPresent(path)
if isfile(path)
    delete(path);
end
end

function removeDirectory(path)
if isfolder(path)
    rmdir(path,"s");
end
end

function cleanupSchedulerTestResources(state)
% Release scheduler-owned references before optionally closing a test pool.
job = state('job');
poolBefore = state('poolBefore');
if isa(job,"rampSpeed.StudyJob") && isscalar(job)
    try
        rampSpeed.StudyExecutor.wait(job,Inf);
    catch
    end
    try
        scheduler = job.scheduler;
        if isstruct(scheduler) && isscalar(scheduler)
            if isfield(scheduler,"activeFutures")
                scheduler.activeFutures = [];
            end
            if isfield(scheduler,"activeIndices")
                scheduler.activeIndices = zeros(0,1);
            end
            if isfield(scheduler,"queue")
                scheduler.queue = [];
            end
            if isfield(scheduler,"pool")
                scheduler.pool = [];
            end
            job.scheduler = scheduler;
        end
        job.future = [];
        delete(job);
    catch
    end
end
state('job') = [];

% StudyExecutor reuses caller-owned pools. Close only a pool created by this
% test file, after its job/future references have been released.
if isempty(poolBefore)
    try
        pool = gcp("nocreate");
        if ~isempty(pool)
            delete(pool);
        end
    catch
    end
end
end
