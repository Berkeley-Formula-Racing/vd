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

function testRunCallbackUsesFunctionHandleBeforeOutputCount(testCase)
root = fileparts(fileparts(mfilename("fullpath")));
archive = fullfile(root,"RampSpeedApp.mlapp");
temporaryRoot = tempname;
mkdir(temporaryRoot);
cleanup = onCleanup(@()rmdir(temporaryRoot,"s")); %#ok<NASGU>
unzip(archive,temporaryRoot);

documentPath = fullfile(temporaryRoot,"matlab","document.xml");
verifyTrue(testCase,isfile(documentPath));
document = fileread(documentPath);
correctPattern = ['parfeval\s*\(\s*backgroundPool\s*,\s*' ...
    '@rampSpeed\.runStudy\s*,\s*2\s*,'];
reversedPattern = ['parfeval\s*\(\s*backgroundPool\s*,\s*2\s*,[\s\S]*?' ...
    '@rampSpeed\.runStudy'];

verifyNotEmpty(testCase,regexp(document,correctPattern,"once"), ...
    "The Run callback must pass the function handle before the output count.");
verifyEmpty(testCase,regexp(document,reversedPattern,"once"), ...
    "The Run callback still uses the reversed parfeval signature.");
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

function testCancelCallbackOwnsCancellationState(testCase)
ensureRampSpeedAppPath();
app = RampSpeedApp("Visible","off");
cleanup = onCleanup(@()deleteIfValid(app));

app.CancelButtonPushed([],[]);

verifyEqual(testCase,string(app.ProgressTextArea.Value(end)), ...
    "No study is running.");
end

function testFixtureRunPreservesOverlayLinesAndWarnings(testCase)
ensureRampSpeedAppPath();
fixture = makeRampFixture();
app = RampSpeedApp("Visible","off");
cleanup = onCleanup(@()deleteIfValid(app));

fakeRunner = @runRenderFixture;
app.runStudyForTest(fixture.cars,fixture.cases,fakeRunner);
lines = findall(app.CapabilityAxes,"Type","line");
markers = string(arrayfun(@(line)line.Marker,lines,"UniformOutput",false));

verifyGreaterThanOrEqual(testCase,numel(lines),2);
hasGap = any(arrayfun(@(line)any(isnan(double(line.YData))),lines));
verifyTrue(testCase,hasGap);
verifyTrue(testCase,any(markers == "o"));
end

function testBaselineSelectionRendersVariantMinusBaseline(testCase)
ensureRampSpeedAppPath();
fixture = makeRampFixture();
app = RampSpeedApp("Visible","off");
cleanup = onCleanup(@()deleteIfValid(app));

app.runStudyForTest(fixture.cars,fixture.cases, ...
    @(car,caseInfo,request,callbacks)makeFixtureRun(caseInfo));
app.BaselineDropDown.Value = "accel";
app.BaselineDropDown.ValueChangedFcn(app.BaselineDropDown,[]);
lines = findall(app.BalanceAxes,"Type","line");
labels = string(arrayfun(@(line)line.DisplayName,lines, ...
    "UniformOutput",false));

verifyGreaterThanOrEqual(testCase,numel(lines),1);
verifyTrue(testCase,any(contains(labels," - accel")));
verifyThat(testCase,string(app.BalanceAxes.YLabel.String), ...
    matlab.unittest.constraints.ContainsSubstring("Delta"));
end

function testUnitSelectionChangesPlotValuesAndLabels(testCase)
ensureRampSpeedAppPath();
fixture = makeRampFixture();
app = RampSpeedApp("Visible","off");
cleanup = onCleanup(@()deleteIfValid(app));

[study,~] = app.runStudyForTest(fixture.cars,fixture.cases, ...
    @(car,caseInfo,request,callbacks)makeFixtureRun(caseInfo));
siLines = findall(app.AeroLoadsAxes,"Type","line");
siLine = siLines(find(arrayfun(@(line)~isempty(line.YData),siLines),1));
siValues = double(siLine.YData);
siLabel = string(app.AeroLoadsAxes.YLabel.String);
app.UnitDropDown.Value = "imperial";
app.UnitDropDown.ValueChangedFcn(app.UnitDropDown,[]);
imperialLines = findall(app.AeroLoadsAxes,"Type","line");
imperialLine = imperialLines(find(arrayfun(@(line)~isempty(line.YData), ...
    imperialLines),1));
imperialValues = double(imperialLine.YData);
imperialLabel = string(app.AeroLoadsAxes.YLabel.String);

verifyNotEqual(testCase,imperialValues,siValues);
verifyThat(testCase,imperialLabel, ...
    matlab.unittest.constraints.ContainsSubstring("lbf"));
verifyNotEqual(testCase,imperialLabel,siLabel);
verifyEqual(testCase,study.runs(1).perSpeed.downforce_N,[200;220],"AbsTol",0);
end

function run = runRenderFixture(~,caseInfo,~,~)
run = makeFixtureRun(caseInfo);
if string(caseInfo.id) == "baseline"
    run.perSpeed.truncated(1) = true;
    run.perSpeed.valid(2) = false;
end
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
