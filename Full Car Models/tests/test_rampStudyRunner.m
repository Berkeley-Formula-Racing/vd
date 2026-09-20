function tests = test_rampStudyRunner
tests = functiontests(localfunctions);
end

function testRunnerRetainsCaseMetadata(testCase)
fixture = makeRampFixture();
cases = fixture.cases;
request = struct("rampType","lateral","settings",struct(), ...
    "parallelRequested",false,"numWorkers",0, ...
    "checkpointPath","", "appVersion","test", ...
    "runCaseFcn",@(car,caseInfo,request,callbacks) ...
        makeFixtureRun(caseInfo));
[study,events] = rampSpeed.runStudy(fixture.cars,cases,request,struct());
verifyEqual(testCase,numel(study.runs),2);
verifyEqual(testCase,study.runs(1).caseId,cases(1).id);
verifyEqual(testCase,study.runs(1).status,"complete");
verifyGreaterThanOrEqual(testCase,numel(events),2);
end

function testRunnerExecutesInOrderAndRetainsCaseFailure(testCase)
fixture = makeRampFixture();
callOrder = strings(0,1);
request = struct("rampType","lateral","settings",struct(), ...
    "parallelRequested",false,"numWorkers",0, ...
    "checkpointPath","", "appVersion","test", ...
    "runCaseFcn",@runInjectedCase);
[study,~] = rampSpeed.runStudy(fixture.cars,fixture.cases,request,struct());
verifyEqual(testCase,callOrder,["baseline";"accel"]);
verifyEqual(testCase,study.runs(1).status,"complete");
verifyEqual(testCase,study.runs(2).status,"failed");
verifyThat(testCase,study.runs(2).runMeta.error.message, ...
    matlab.unittest.constraints.ContainsSubstring("injected failure"));
verifyThat(testCase,study.runs(2).runMeta.error.identifier, ...
    matlab.unittest.constraints.ContainsSubstring("test:caseFailure"));

    function run = runInjectedCase(~,caseInfo,~,~)
        callOrder(end+1,1) = string(caseInfo.id);
        if string(caseInfo.id) == "accel"
            error("test:caseFailure","injected failure");
        end
        run = makeFixtureRun(caseInfo);
    end
end

function testRunnerPublishesRequiredProgressFields(testCase)
fixture = makeRampFixture();
events = struct.empty;
request = struct("rampType","lateral","settings",struct(), ...
    "parallelRequested",false,"numWorkers",0, ...
    "checkpointPath","", "appVersion","test", ...
    "runCaseFcn",@(car,caseInfo,request,callbacks) ...
        makeFixtureRun(caseInfo));
callbacks = struct("onProgress",@capture);
[~,returnedEvents] = rampSpeed.runStudy(fixture.cars,fixture.cases, ...
    request,callbacks);
verifyGreaterThanOrEqual(testCase,numel(events),2);
verifyEqual(testCase,returnedEvents,events);
required = ["phase","caseId","speedIndex","speed_mps", ...
    "completedCases","totalCases","message"];
verifyTrue(testCase,all(isfield(events,cellstr(required))));
verifyEqual(testCase,events(1).totalCases,2);

    function capture(event)
        if isempty(events)
            events = event;
        else
            events(end+1,1) = event;
        end
    end
end

function testParallelRequestWithNoWorkerCountFallsBackToSerial(testCase)
fixture = makeRampFixture();
request = struct("rampType","lateral","settings",struct(), ...
    "parallelRequested",true,"numWorkers",0, ...
    "checkpointPath","", "appVersion","test", ...
    "runCaseFcn",@(car,caseInfo,request,callbacks) ...
        makeFixtureRun(caseInfo));
[study,~] = rampSpeed.runStudy(fixture.cars,fixture.cases,request,struct());
effectiveWorkers = arrayfun(@(run)run.runMeta.effectiveWorkers,study.runs);
fallbackReasons = arrayfun(@(run)string(run.runMeta.parallelFallbackReason), ...
    study.runs);
statuses = arrayfun(@(run)string(run.status),study.runs);
verifyEqual(testCase,effectiveWorkers,[0;0]);
verifyTrue(testCase,all(strlength(fallbackReasons) > 0));
verifyTrue(testCase,all(statuses == "complete"));
end

function testSetupCatalogHasStableRoleAwareMetadata(testCase)
cars = {struct("name","lap-1"),struct("name","accel-1"); ...
    struct("name","lap-2"),struct("name","accel-2")};
designTable = table([10;20],[1;2], ...
    'VariableNames',{'mass','index'});
lapCases = rampSpeed.setupCaseCatalog(cars,designTable,"lap");
accelCases = rampSpeed.setupCaseCatalog(cars,designTable,"acceleration");
lapCasesAgain = rampSpeed.setupCaseCatalog(cars,designTable,"lap");
verifyEqual(testCase,numel(lapCases),2);
verifyEqual(testCase,string({lapCases.id}),string({lapCasesAgain.id}));
verifyEqual(testCase,[lapCases.designRow],[1 2]);
verifyEqual(testCase,[lapCases.sourceIndex],[1 2]);
verifyEqual(testCase,string({lapCases.carRole}),["lap","lap"]);
verifyEqual(testCase,[lapCases.carColumn],[1 1]);
verifyEqual(testCase,string({accelCases.carRole}),["acceleration","acceleration"]);
verifyEqual(testCase,[accelCases.carColumn],[2 2]);
verifyTrue(testCase,all(strlength(string({lapCases.label})) > 0));
end

function testMalformedInjectedRunIsRetainedAndStudyContinues(testCase)
fixture = makeRampFixture();
checkpointPath = fullfile(tempdir,"ramp-study-malformed-run-test.mat");
if isfile(checkpointPath)
    delete(checkpointPath);
end
cleanup = onCleanup(@()deleteIfPresent(checkpointPath));
request = struct("rampType","lateral","settings",struct(), ...
    "parallelRequested",false,"numWorkers",0, ...
    "checkpointPath",checkpointPath,"appVersion","test", ...
    "runCaseFcn",@returnMalformed);
[study,~] = rampSpeed.runStudy(fixture.cars,fixture.cases,request,struct());

verifyEqual(testCase,study.runs(1).status,"failed");
verifyEqual(testCase,study.runs(1).runMeta.error.identifier, ...
    "rampSpeed:invalidRun");
verifyThat(testCase,study.runs(1).runMeta.error.message, ...
    matlab.unittest.constraints.ContainsSubstring("scalar run struct"));
verifyEqual(testCase,study.runs(2).status,"complete");
verifyTrue(testCase,isfile(checkpointPath));

    function run = returnMalformed(~,caseInfo,~,~)
        if string(caseInfo.id) == "baseline"
            run = [];
        else
            run = makeFixtureRun(caseInfo);
        end
    end
end

function deleteIfPresent(fileName)
if isfile(fileName)
    delete(fileName);
end
end
function testSetupCatalogRejectsDuplicateIds(testCase)
cars = {struct("name","lap-1");struct("name","lap-2")};
designTable = table(["duplicate";"duplicate"], ...
    'VariableNames',{'id'});
verifyError(testCase,@()rampSpeed.setupCaseCatalog(cars,designTable), ...
    "rampSpeed:duplicateCaseId");
end

function testSetupCatalogRejectsVectorDesignRowMismatch(testCase)
cars = {struct("name","lap-1"),struct("name","accel-1")};
oneDesignRow = table("one",'VariableNames',{'id'});
oneSetup = rampSpeed.setupCaseCatalog(cars,oneDesignRow);
verifyEqual(testCase,numel(oneSetup),1);
twoDesignRows = table(["one";"two"],'VariableNames',{'id'});
verifyError(testCase,@()rampSpeed.setupCaseCatalog(cars,twoDesignRows), ...
    "rampSpeed:carCellDesignMismatch");
end

function testRunnerRejectsDuplicateCaseIdsBeforeExecution(testCase)
fixture = makeRampFixture();
cases = fixture.cases;
cases(2).id = cases(1).id;
request = struct("rampType","lateral","settings",struct(), ...
    "parallelRequested",false,"numWorkers",0, ...
    "checkpointPath","","appVersion","test", ...
    "runCaseFcn",@unexpectedExecution);
verifyError(testCase,@()rampSpeed.runStudy(fixture.cars,cases,request,struct()), ...
    "rampSpeed:duplicateCaseId");

    function unexpectedExecution(varargin)
        error("test:unexpectedExecution","duplicate cases must be rejected first.");
    end
end

function testDirectRunCaseProgressCallbackIsAdvisory(testCase)
[cars,~] = carConfig();
settings = struct("speeds",5,"nRamp",1,"nBisect",0, ...
    "mode","coast","verbose",false);
request = struct("rampType","lateral","settings",settings, ...
    "progressQueue",@captureQueue);
callbacks = struct("onProgress",@throwProgress);
run = rampSpeed.runCase(cars{1,1}, ...
    struct("id","progress-callback","label","progress-callback", ...
    "carRole","lap"),request,callbacks);
verifyTrue(testCase,isstruct(run));
verifyEqual(testCase,run.caseId,"progress-callback");

    function throwProgress(~)
        error("test:progressCallback","advisory callback failure");
    end

    function captureQueue(~)
    end
end
function testParallelWorkerFailureDoesNotRerunActiveCases(testCase)
fixture = makeRampFixture();
cases = fixture.cases;
cases(1).id = "fail";
cases(1).label = "fail";
cases(3) = cases(2);
cases(3).id = "third";
cases(3).label = "third";
logPath = string(tempname) + ".log";
cleanup = onCleanup(@()deleteIfPresent(logPath));
request = struct("rampType","lateral","settings",struct(), ...
    "parallelRequested",true,"numWorkers",2, ...
    "checkpointPath","","appVersion","test", ...
    "runCaseFcn",@parallelInjectedWorkerFailure, ...
    "workerLogPath",logPath);
[study,~] = rampSpeed.runStudy(fixture.cars,cases,request,struct());
logLines = splitlines(strtrim(string(fileread(logPath))));
verifyEqual(testCase,sum(logLines == "fail"),1);
verifyEqual(testCase,study.runs(1).status,"failed");
verifyEqual(testCase,study.runs(1).runMeta.error.identifier, ...
    "test:workerFailure");
verifyThat(testCase,study.runs(1).runMeta.error.message, ...
    matlab.unittest.constraints.ContainsSubstring("injected worker failure"));
verifyTrue(testCase,~isempty(study.runs(1).runMeta.error.stack));
verifyEqual(testCase,study.runs(2).status,"complete");
verifyEqual(testCase,study.runs(3).status,"complete");
verifyGreaterThan(testCase,study.runs(2).runMeta.effectiveWorkers,0);
end

function run = parallelInjectedWorkerFailure(~,caseInfo,request,~)
logPath = char(request.workerLogPath);
fid = fopen(logPath,"a");
if fid < 0
    error("test:logOpen","Could not open worker log.");
end
fprintf(fid,"%s\n",char(string(caseInfo.id)));
fclose(fid);
if string(caseInfo.id) == "fail"
    error("test:workerFailure","injected worker failure");
end
run = rampSpeed.makeRun("lateral","coast",request.settings,caseInfo);
end

function testParallelUnserializableRequestFallsBackToSerial(testCase)
fixture = makeRampFixture();
badQueue = parallel.pool.PollableDataQueue;
cleanup = onCleanup(@()close(badQueue));
request = struct("rampType","lateral","settings",struct(), ...
    "parallelRequested",true,"numWorkers",2, ...
    "checkpointPath","","appVersion","test", ...
    "runCaseFcn",@serializableFallbackRun, ...
    "progressQueue",badQueue);
[study,~] = rampSpeed.runStudy(fixture.cars,fixture.cases,request,struct());
effectiveWorkers = arrayfun(@(run)run.runMeta.effectiveWorkers,study.runs);
fallbackReasons = arrayfun(@(run)string(run.runMeta.parallelFallbackReason), ...
    study.runs);
verifyEqual(testCase,effectiveWorkers,[0;0]);
verifyTrue(testCase,all(contains(fallbackReasons,"serializable")));
verifyEqual(testCase,string({study.runs.status}),["complete","complete"]);

end

function run = serializableFallbackRun(~,caseInfo,request,~)
run = rampSpeed.makeRun("lateral","coast",request.settings,caseInfo);
end

function testInvalidCanonicalInjectedRunIsRetainedAndStudyContinues(testCase)
fixture = makeRampFixture();
checkpointPath = fullfile(tempdir,"ramp-study-invalid-canonical-test.mat");
cleanup = onCleanup(@()deleteIfPresent(checkpointPath));
request = struct("rampType","lateral","settings",struct(), ...
    "parallelRequested",false,"numWorkers",0, ...
    "checkpointPath",checkpointPath,"appVersion","test", ...
    "runCaseFcn",@returnInvalidCanonical);
[study,~] = rampSpeed.runStudy(fixture.cars,fixture.cases,request,struct());
verifyEqual(testCase,study.runs(1).status,"failed");
verifyEqual(testCase,study.runs(1).runMeta.error.identifier, ...
    "rampSpeed:invalidRun");
verifyThat(testCase,study.runs(1).runMeta.error.message, ...
    matlab.unittest.constraints.ContainsSubstring("Invalid canonical run"));
verifyEqual(testCase,study.runs(2).status,"complete");
verifyTrue(testCase,isfile(checkpointPath));

    function run = returnInvalidCanonical(~,caseInfo,~,~)
        run = makeFixtureRun(caseInfo);
        if string(caseInfo.id) == "baseline"
            run.perSpeed.injected_bad = false(height(run.perSpeed),1);
        end
    end
end

function testMissingCanonicalColumnsAreRejectedAndStudyContinues(testCase)
fixture = makeRampFixture();
checkpointPath = fullfile(tempdir,"ramp-study-missing-canonical-columns-test.mat");
cleanup = onCleanup(@()deleteIfPresent(checkpointPath));
request = struct("rampType","lateral","settings",struct(), ...
    "parallelRequested",false,"numWorkers",0, ...
    "checkpointPath",checkpointPath,"appVersion","test", ...
    "runCaseFcn",@returnMissingCanonical);
[study,~] = rampSpeed.runStudy(fixture.cars,fixture.cases,request,struct());
verifyEqual(testCase,study.runs(1).status,"failed");
verifyTrue(testCase,isfield(study.runs(1).runMeta,'error'));
if isfield(study.runs(1).runMeta,'error')
    verifyEqual(testCase,study.runs(1).runMeta.error.identifier, ...
        "rampSpeed:invalidRun");
    verifyThat(testCase,study.runs(1).runMeta.error.message, ...
        matlab.unittest.constraints.ContainsSubstring("canonical"));
end
verifyEqual(testCase,study.runs(2).status,"complete");
verifyTrue(testCase,isfile(checkpointPath));

    function run = returnMissingCanonical(~,caseInfo,~,~)
        run = makeFixtureRun(caseInfo);
        if string(caseInfo.id) == "baseline"
            run.perSpeed = run.perSpeed(:,"speed_mps");
            run.points = run.points(:,"speed_mps");
        end
    end
end

function testCanonicalColumnClassesAreRejectedAndStudyContinues(testCase)
fixture = makeRampFixture();
checkpointPath = fullfile(tempdir,"ramp-study-canonical-class-test.mat");
cleanup = onCleanup(@()deleteIfPresent(checkpointPath));
request = struct("rampType","lateral","settings",struct(), ...
    "parallelRequested",false,"numWorkers",0, ...
    "checkpointPath",checkpointPath,"appVersion","test", ...
    "runCaseFcn",@returnWrongCanonicalClasses);
[study,~] = rampSpeed.runStudy(fixture.cars,fixture.cases,request,struct());
verifyEqual(testCase,study.runs(1).status,"failed");
verifyTrue(testCase,isfield(study.runs(1).runMeta,'error'));
if isfield(study.runs(1).runMeta,'error')
    verifyEqual(testCase,study.runs(1).runMeta.error.identifier, ...
        "rampSpeed:invalidRun");
    verifyThat(testCase,study.runs(1).runMeta.error.message, ...
        matlab.unittest.constraints.ContainsSubstring("valid"));
end
verifyEqual(testCase,study.runs(2).status,"complete");
verifyTrue(testCase,isfile(checkpointPath));

    function run = returnWrongCanonicalClasses(~,caseInfo,~,~)
        run = makeFixtureRun(caseInfo);
        if string(caseInfo.id) == "baseline"
            run.perSpeed.valid = double(run.perSpeed.valid);
        end
    end
end

function testRunnerSelectsDistinctCarsFromRowVector(testCase)
cars = {struct("marker","one"),struct("marker","two")};
cases(1) = struct('id',"one",'label',"one",'source',"test", ...
    'designRow',1,'sourceIndex',1,'carRole',"lap",'carColumn',1);
cases(2) = struct('id',"two",'label',"two",'source',"test", ...
    'designRow',2,'sourceIndex',2,'carRole',"lap",'carColumn',1);
markers = strings(0,1);
request = struct("rampType","lateral","settings",struct(), ...
    "parallelRequested",false,"numWorkers",0, ...
    "checkpointPath","","appVersion","test", ...
    "runCaseFcn",@recordCar);
[study,~] = rampSpeed.runStudy(cars,cases,request,struct());
verifyEqual(testCase,markers,["one";"two"]);
verifyEqual(testCase,string({study.runs.status}),["complete","complete"]);

    function run = recordCar(car,caseInfo,~,~)
        markers(end+1,1) = string(car.marker);
        run = makeFixtureRun(caseInfo);
    end
end

function testParallelInterruptRetainsFinishedFutureResult(testCase)
cars = {struct("marker","fast");struct("marker","late")};
cases(1) = struct('id',"fast",'label',"fast",'source',"test", ...
    'designRow',1,'sourceIndex',1,'carRole',"lap",'carColumn',1);
cases(2) = struct('id',"late",'label',"late",'source',"test", ...
    'designRow',2,'sourceIndex',2,'carRole',"lap",'carColumn',1);
cancelled = false;
callbacks = struct("onProgress",@captureProgress,"isCancelled",@isCancelled);
request = struct("rampType","lateral","settings",struct(), ...
    "parallelRequested",true,"numWorkers",2, ...
    "checkpointPath","","appVersion","test", ...
    "runCaseFcn",@runInterruptedCase);
[study,~] = rampSpeed.runStudy(cars,cases,request,callbacks);
verifyEqual(testCase,study.status,"cancelled");
verifyEqual(testCase,study.runs(1).status,"complete");
verifyEqual(testCase,study.runs(2).status,"failed");
if isfield(study.runs(2).runMeta,'error')
    verifyEqual(testCase,study.runs(2).runMeta.error.identifier, ...
        "test:lateFailure");
end

    function captureProgress(event)
        if string(event.caseId) == "fast" && ...
                contains(string(event.message),"complete")
            cancelled = true;
        end
    end

    function value = isCancelled()
        value = cancelled;
    end

    function run = runInterruptedCase(~,caseInfo,~,~)
        if string(caseInfo.id) == "late"
            pause(1.0);
            run = makeFixtureRun(caseInfo);
            run.status = "failed";
            run.runMeta.status = "failed";
            run.runMeta.errors = "late future failure";
            run.runMeta.error = struct("identifier","test:lateFailure", ...
                "message","late future failure","stack",struct.empty(0,1));
            return
        end
        run = makeFixtureRun(caseInfo);
    end
end
