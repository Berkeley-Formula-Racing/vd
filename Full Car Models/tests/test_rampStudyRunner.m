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
