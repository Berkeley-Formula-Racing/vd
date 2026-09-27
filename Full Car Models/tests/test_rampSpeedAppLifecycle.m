function tests = test_rampSpeedAppLifecycle
tests = functiontests(localfunctions);
end

function testAppOwnsSessionController(testCase)
ensureRampSpeedAppPath();
app = RampSpeedApp("Visible","off");
cleanup = onCleanup(@()deleteIfValid(app)); %#ok<NASGU>

verifyTrue(testCase,isprop(app,"Session"));
verifyClass(testCase,app.Session,"rampSpeed.RampSpeedSession");
verifyEqual(testCase,string(app.Session.state),"idle");
verifyEqual(testCase,string(app.Session.rampType), ...
    string(app.RampTypeDropDown.Value));
end

function testFixtureRunUsesSessionLifecycle(testCase)
ensureRampSpeedAppPath();
fixture = makeRampFixture();
app = RampSpeedApp("Visible","off");
cleanup = onCleanup(@()deleteIfValid(app)); %#ok<NASGU>

fakeRunner = @(car,caseInfo,request,callbacks)makeFixtureRun(caseInfo);
[study,events] = app.runStudyForTest(fixture.cars,fixture.cases,fakeRunner);

verifyEqual(testCase,string(study.status),"complete");
verifyEqual(testCase,string(app.Session.state),"completed");
verifyTrue(testCase,any(string(app.Session.stateHistory) == "running"));
verifyTrue(testCase,any(string(app.Session.stateHistory) == "completed"));
verifyEqual(testCase,string(app.Session.study.status),"complete");
verifyGreaterThanOrEqual(testCase,numel(events),numel(fixture.cases));
verifyEqual(testCase,string(app.Session.progress.status),"completed");
end

function testCancelAndClearDelegateToSession(testCase)
ensureRampSpeedAppPath();
app = RampSpeedApp("Visible","off");
cleanup = onCleanup(@()deleteIfValid(app)); %#ok<NASGU>

app.CancelButtonPushed([],[]);
verifyEqual(testCase,string(app.Session.state),"idle");

app.ClearButtonPushed([],[]);
verifyEqual(testCase,string(app.Session.state),"idle");
verifyTrue(testCase,isempty(app.Session.study.runs));
end

function testSessionLocksSetupChangesWhileRunning(testCase)
ensureRampSpeedAppPath();
fixture = makeRampFixture();
app = RampSpeedApp("Visible","off");
cleanup = onCleanup(@()deleteIfValid(app)); %#ok<NASGU>

app.Session = rampSpeed.RampSpeedSession(fixture.cars,fixture.cases, ...
    struct("rampType","longitudinal","settings",struct(), ...
    "runCaseFcn",@(car,caseInfo,request,callbacks)makeFixtureRun(caseInfo)));
app.Session.start(struct("parallelRequested",false));
verifyEqual(testCase,string(app.Session.state),"completed");
verifyError(testCase,@()app.Session.editSetup("baseline",struct("label","x")), ...
    "rampSpeed:baselineProtected");
end

function testReadOnlySessionCannotStart(testCase)
ensureRampSpeedAppPath();
app = RampSpeedApp("Visible","off");
cleanup = onCleanup(@()deleteIfValid(app)); %#ok<NASGU>

app.Session.acceptResult(struct("status","complete","readOnly",true));
verifyTrue(testCase,app.Session.readOnly);
verifyError(testCase,@()app.Session.start(struct()), ...
    "rampSpeed:readOnlySession");
end

function ensureRampSpeedAppPath()
root = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(root));
sourceRoot = fullfile(root,".task7");
if isfolder(sourceRoot)
    addpath(sourceRoot);
end
end

function deleteIfValid(app)
if ~isempty(app) && isvalid(app)
    delete(app);
end
end
