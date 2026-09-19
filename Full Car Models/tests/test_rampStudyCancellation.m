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

function deleteIfPresent(fileName)
if isfile(fileName)
    delete(fileName);
end
end
