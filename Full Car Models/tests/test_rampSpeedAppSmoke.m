function tests = test_rampSpeedAppSmoke
tests = functiontests(localfunctions);
end

function testAppExposesPersistentSetupAndAnalysisShell(testCase)
ensureRampSpeedAppPath();
app = RampSpeedApp("Visible","off");
cleanup = onCleanup(@()deleteIfValid(app));

required = ["SetupTable","CarRoleDropDown","RampTypeDropDown", ...
    "LateralModeDropDown","SpeedStartEditField","SpeedStopEditField", ...
    "SpeedStepEditField","SpeedVectorEditField","NRampEditField", ...
    "NBisectEditField","ResidualToleranceEditField", ...
    "ExecutionModeDropDown","WorkersEditField","UnitDropDown", ...
    "BaselineDropDown","RunButton","CancelButton","LoadButton", ...
    "SaveButton","ExportButton","ClearButton","ProgressTextArea", ...
    "MainTabGroup","CapabilityTab","BalanceTab","AeroLoadsTab", ...
    "SuspensionTab","RawRampTab","InspectorDataTab"];
verifyTrue(testCase,all(arrayfun(@(name)isprop(app,name),required)));

verifyEqual(testCase,string(app.RampTypeDropDown.Value),"lateral");
verifyEqual(testCase,string(app.LateralModeDropDown.Value),"coast");
verifyEqual(testCase,string(app.UnitDropDown.Value),"SI");
verifyEqual(testCase,string(app.CarRoleDropDown.Value),"lap");
verifyGreaterThan(testCase,height(app.SetupTable.Data),0);

tabTitles = string(arrayfun(@(child)child.Title, ...
    app.MainTabGroup.Children,"UniformOutput",false));
verifyTrue(testCase,all(ismember(["Capability","Balance","Aero & Loads", ...
    "Suspension","Raw Ramp","Inspector/Data"],tabTitles)));
end

function testFixtureRunUsesInjectedRunnerWithoutSolver(testCase)
ensureRampSpeedAppPath();
fixture = makeRampFixture();
app = RampSpeedApp("Visible","off");
cleanup = onCleanup(@()deleteIfValid(app));

fakeRunner = @(car,caseInfo,request,callbacks)makeFixtureRun(caseInfo);
[study,events] = app.runStudyForTest(fixture.cars,fixture.cases,fakeRunner);

verifyEqual(testCase,string(study.status),"complete");
verifyEqual(testCase,numel(study.runs),numel(fixture.cases));
verifyTrue(testCase,all(string({study.runs.status}) == "complete"));
verifyGreaterThanOrEqual(testCase,numel(events),numel(fixture.cases));
verifyEqual(testCase,string(app.ProgressTextArea.Value(end)), ...
    "Study complete.");
end

function ensureRampSpeedAppPath()
root = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(root));
end

function deleteIfValid(app)
if ~isempty(app) && isvalid(app)
    delete(app);
end
end
