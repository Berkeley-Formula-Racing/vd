function tests = test_rampStudyCancellation
tests = functiontests(localfunctions);
end

function testCancellationRetainsCompletedCaseAndCheckpoint(testCase)
fixture = makeRampFixture();
checkpointPath = fullfile(tempdir,"ramp-study-cancellation-test.mat");
if isfile(checkpointPath)
    delete(checkpointPath);
end
cleanup = onCleanup(@()deleteIfPresent(checkpointPath));

cancelled = false;
callbacks = struct("isCancelled",@isCancelled);
request = struct("rampType","lateral","settings",struct(), ...
    "parallelRequested",false,"numWorkers",0, ...
    "checkpointPath",checkpointPath,"appVersion","test", ...
    "runCaseFcn",@runCase);
[study,events] = rampSpeed.runStudy(fixture.cars,fixture.cases, ...
    request,callbacks);

verifyEqual(testCase,study.runs(1).status,"complete");
verifyEqual(testCase,study.runs(2).status,"cancelled");
verifyEqual(testCase,study.status,"cancelled");
verifyTrue(testCase,isfile(checkpointPath));
verifyTrue(testCase,any(strcmp(string({events.message}), ...
    "study cancelled")));

    function run = runCase(~,caseInfo,~,~)
        run = makeFixtureRun(caseInfo);
        cancelled = true;
    end

    function value = isCancelled()
        value = cancelled;
    end
end

function testParallelWorkerCancellationUsesWorkerQueue(testCase)
fixture = makeRampFixture();
cases = fixture.cases;
cases(1).id = "first";
cases(1).label = "first";
cases(3) = cases(2);
cases(3).id = "third";
cases(3).label = "third";
checkpointPath = fullfile(tempdir,"ramp-study-parallel-cancellation-test.mat");
logPath = string(tempname) + ".log";
cleanupCheckpoint = onCleanup(@()deleteIfPresent(checkpointPath));
cleanupLog = onCleanup(@()deleteIfPresent(logPath));
cancelled = false;
callbacks = struct("onProgress",@captureProgress,"isCancelled",@isCancelled);
request = struct("rampType","lateral","settings",struct(), ...
    "parallelRequested",true,"numWorkers",2, ...
    "checkpointPath",checkpointPath,"appVersion","test", ...
    "runCaseFcn",@parallelCancellationCase,"workerLogPath",logPath);
[study,~] = rampSpeed.runStudy(fixture.cars,cases,request,callbacks);
logLines = splitlines(strtrim(string(fileread(logPath))));
verifyEqual(testCase,study.status,"cancelled");
verifyEqual(testCase,study.runs(1).status,"cancelled");
verifyFalse(testCase,any(logLines == "first:complete"));
verifyTrue(testCase,any(logLines == "first:cancelled"));
verifyTrue(testCase,isfile(checkpointPath));

    function captureProgress(event)
        if string(event.caseId) == "first"
            cancelled = true;
        end
    end

    function value = isCancelled()
        value = cancelled;
    end
end
function deleteIfPresent(fileName)
if isfile(fileName)
    delete(fileName);
end
end

function run = parallelCancellationCase(~,caseInfo,request,callbacks)
if ~isempty(callbacks.onProgress)
    callbacks.onProgress(struct("speedIndex",1,"speed_mps",5, ...
        "message","started"));
end
for k = 1:50
    if callbacks.isCancelled()
        fid = fopen(char(request.workerLogPath),"a");
        fprintf(fid,"%s\n",char(string(caseInfo.id) + ":cancelled"));
        fclose(fid);
        run = rampSpeed.makeRun("lateral","coast",struct(),caseInfo);
        run.status = "cancelled";
        run.runMeta.status = "cancelled";
        return
    end
    pause(0.05);
end
fid = fopen(char(request.workerLogPath),"a");
fprintf(fid,"%s\n",char(string(caseInfo.id) + ":complete"));
fclose(fid);
run = rampSpeed.makeRun("lateral","coast",struct(),caseInfo);
run.status = "complete";
end
