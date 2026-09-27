function tests = test_rampSpeedSession
%TEST_RAMPSPEEDSESSION Task 6 lifecycle contracts.

tests = functiontests(localfunctions);
end

function testStudyExecutorSerialLifecycle(testCase)
fixture = makeRampFixture();
request = baseRequest(@fixtureRun);
job = rampSpeed.StudyExecutor.start(fixture.cars,fixture.cases, ...
    request,struct());

verifyEqual(testCase,job.state,"completed");
verifyTrue(testCase,isempty(job.future));
verifyEqual(testCase,job.study.status,"complete");
verifyEqual(testCase,numel(job.study.runs),2);
verifyGreaterThanOrEqual(testCase,numel(job.progressEvents),2);
end

function testStudyExecutorCancellationRetainsTerminalRows(testCase)
fixture = makeRampFixture();
jobRef = [];
request = baseRequest(@fixtureRun);
callbacks = struct( ...
    "onJobCreated",@rememberJob, ...
    "onProgress",@cancelAfterFirstCase);

job = rampSpeed.StudyExecutor.start(fixture.cars,fixture.cases, ...
    request,callbacks);

verifyEqual(testCase,job.state,"cancelled");
verifyEqual(testCase,job.study.status,"cancelled");
verifyEqual(testCase,string({job.study.runs.status}), ...
    ["complete","cancelled"]);

    function rememberJob(created)
        jobRef = created;
    end

    function cancelAfterFirstCase(event)
        if isempty(jobRef) || ~isstruct(event)
            return
        end
        if string(event.phase) == "case" && ...
                contains(string(event.message),"complete")
            rampSpeed.StudyExecutor.cancel(jobRef);
        end
    end
end

function testSessionSetupProtectionAndSelection(testCase)
fixture = makeRampFixture();
session = rampSpeed.RampSpeedSession(fixture.cars,fixture.cases, ...
    baseOptions(@fixtureRun));

verifyError(testCase,@()session.deleteSetup("baseline"), ...
    "rampSpeed:baselineProtected");
verifyError(testCase,@()session.editSetup("baseline", ...
    struct("label","changed")),"rampSpeed:baselineProtected");

session.duplicateSetup("accel","accel-copy","Acceleration copy");
model = session.viewModel();
verifyEqual(testCase,model.state,"idle");
verifyTrue(testCase,any(string(model.setupTable.id) == "accel-copy"));

session.selectCases(["baseline","accel-copy"]);
verifyEqual(testCase,session.viewModel().selectedCaseIds, ...
    ["baseline";"accel-copy"]);

session.deleteSetup("accel-copy");
verifyFalse(testCase,any(string(session.viewModel().setupTable.id) == ...
    "accel-copy"));
end

function testSessionLocksEditsDuringRun(testCase)
fixture = makeRampFixture();
busyIdentifier = "";
options = baseOptions(@fixtureRun);
options.onProgress = @captureProgress;
session = rampSpeed.RampSpeedSession(fixture.cars,fixture.cases,options);

job = session.start();
verifyEqual(testCase,job.state,"completed");
verifyEqual(testCase,session.state,"completed");
verifyEqual(testCase,busyIdentifier,"rampSpeed:sessionBusy");

    function captureProgress(event)
        if ~isstruct(event) || string(event.phase) ~= "case"
            return
        end
        if strlength(string(busyIdentifier)) > 0
            return
        end
        try
            session.editSetup("accel",struct("label","too soon"));
        catch ME
            busyIdentifier = string(ME.identifier);
        end
    end
end

function testSessionCancelTransitionsToCancelled(testCase)
fixture = makeRampFixture();
jobRef = [];
options = baseOptions(@fixtureRun);
options.onJobCreated = @rememberJob;
options.onProgress = @cancelFromProgress;
session = rampSpeed.RampSpeedSession(fixture.cars,fixture.cases,options);

job = session.start();
verifyEqual(testCase,job.state,"cancelled");
verifyEqual(testCase,session.state,"cancelled");
verifyEqual(testCase,session.viewModel().study.status,"cancelled");

    function rememberJob(created)
        jobRef = created;
    end

    function cancelFromProgress(event)
        if isempty(jobRef) || ~isstruct(event)
            return
        end
        if string(event.phase) == "case" && ...
                contains(string(event.message),"complete")
            session.cancel();
        end
    end
end

function testSessionClearResultsReturnsToIdle(testCase)
fixture = makeRampFixture();
session = rampSpeed.RampSpeedSession(fixture.cars,fixture.cases, ...
    baseOptions(@fixtureRun));
session.start();
verifyEqual(testCase,session.state,"completed");

session.clearResults();
model = session.viewModel();
verifyEqual(testCase,model.state,"idle");
verifyTrue(testCase,isempty(model.study.runs));
end

function request = baseRequest(runCaseFcn)
request = struct( ...
    "rampType","lateral", ...
    "settings",struct(), ...
    "parallelRequested",false, ...
    "numWorkers",0, ...
    "checkpointPath","", ...
    "appVersion","task6-test", ...
    "runCaseFcn",runCaseFcn);
end

function options = baseOptions(runCaseFcn)
options = baseRequest(runCaseFcn);
options.onProgress = [];
options.onJobCreated = [];
end

function run = fixtureRun(~,caseInfo,~,~)
run = makeFixtureRun(caseInfo);
end
