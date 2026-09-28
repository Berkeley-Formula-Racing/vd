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
    "MainTabGroup","BalanceTab","AeroLoadsTab", ...
    "SuspensionTab","RawRampTab","InspectorDataTab"];
verifyTrue(testCase,all(arrayfun(@(name)isprop(app,name),required)));
verifyFalse(testCase,isprop(app,"CapabilityTab"));
verifyEqual(testCase,numel(app.BalanceAxes),4);
verifyEqual(testCase,numel(app.SuspensionAxes),4);
verifyEqual(testCase,string(app.SpeedGridModeDropDown.Items), ...
    ["Preview","Accurate","High accuracy"]);

verifyEqual(testCase,string(app.RampTypeDropDown.Value),"lateral");
verifyEqual(testCase,string(app.LateralModeDropDown.Value),"coast");
verifyEqual(testCase,string(app.UnitDropDown.Value),"SI");
verifyEqual(testCase,string(app.CarRoleDropDown.Value),"lap");
verifyGreaterThan(testCase,height(app.SetupTable.Data),0);

tabTitles = string(arrayfun(@(child)child.Title, ...
    app.MainTabGroup.Children,"UniformOutput",false));
verifyFalse(testCase,any(tabTitles == "Capability"));
verifyTrue(testCase,all(ismember(["Balance","Aero & Loads", ...
    "Suspension","Raw Ramp","Inspector/Data"],tabTitles)));
end

function testRunCallbackUsesSessionExecutor(testCase)
root = fileparts(fileparts(mfilename("fullpath")));
archive = fullfile(root,"RampSpeedApp.mlapp");
temporaryRoot = tempname;
mkdir(temporaryRoot);
cleanup = onCleanup(@()rmdir(temporaryRoot,"s")); %#ok<NASGU>
unzip(archive,temporaryRoot);

documentPath = fullfile(temporaryRoot,"matlab","document.xml");
verifyTrue(testCase,isfile(documentPath));
document = fileread(documentPath);
runSection = regexp(document, ...
    'function RunButtonPushed[\s\S]*?function CancelButtonPushed', ...
    'match','once');
cancelSection = regexp(document, ...
    'function CancelButtonPushed[\s\S]*?function LoadButtonPushed', ...
    'match','once');
verifyNotEmpty(testCase,runSection);
verifyNotEmpty(testCase,regexp(runSection,'app\.Session\.start',"once"));
verifyEmpty(testCase,regexp(runSection, ...
    'parfeval|parallel\.pool\.DataQueue',"once"));
verifyNotEmpty(testCase,regexp(cancelSection,'app\.Session\.cancel',"once"));
verifyEmpty(testCase,regexp(cancelSection, ...
    'cancel\s*\(\s*app\.Future',"once"));
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
lines = findall(app.BalanceAxes(1),"Type","line");
markers = string(arrayfun(@(line)line.Marker,lines,"UniformOutput",false));

verifyGreaterThanOrEqual(testCase,numel(lines),2);
hasGap = any(arrayfun(@(line)any(isnan(double(line.YData))),lines));
verifyTrue(testCase,hasGap);
verifyTrue(testCase,any(markers == "o"));
end

function testAllSelectedSetupsOverlayOnBalanceAndSuspensionTabs(testCase)
ensureRampSpeedAppPath();
fixture = makeRampFixture();
app = RampSpeedApp("Visible","off");
cleanup = onCleanup(@()deleteIfValid(app));

app.runStudyForTest(fixture.cars,fixture.cases, ...
    @(car,caseInfo,request,callbacks)makeFixtureRun(caseInfo));
balanceLines = gobjects(0,1);
for ax = app.BalanceAxes(:).'
    balanceLines = [balanceLines; findall(ax,"Type","line")]; %#ok<AGROW>
end
suspensionLines = gobjects(0,1);
for ax = app.SuspensionAxes(:).'
    suspensionLines = [suspensionLines; findall(ax,"Type","line")]; %#ok<AGROW>
end
balanceLabels = string(arrayfun(@(line)line.DisplayName,balanceLines, ...
    "UniformOutput",false));
suspensionLabels = string(arrayfun(@(line)line.DisplayName,suspensionLines, ...
    "UniformOutput",false));

verifyGreaterThanOrEqual(testCase,numel(balanceLines),2);
verifyGreaterThanOrEqual(testCase,numel(suspensionLines),2);
verifyTrue(testCase,any(contains(balanceLabels,"baseline")));
verifyTrue(testCase,any(contains(balanceLabels,"acceleration")));
verifyTrue(testCase,any(contains(suspensionLabels,"baseline")));
verifyTrue(testCase,any(contains(suspensionLabels,"acceleration")));
verifyEqual(testCase,string(app.BalanceAxes(1).Title.String), ...
    "Total balance");
verifyEqual(testCase,string(app.BalanceAxes(4).Title.String), ...
    "Steer angle");
verifyEqual(testCase,string(app.SuspensionAxes(3).Title.String), ...
    "Pitch angle");
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

function testClearRemovesAllPlotGraphicsIncludingHiddenWarnings(testCase)
ensureRampSpeedAppPath();
fixture = makeRampFixture();
app = RampSpeedApp("Visible","off");
cleanup = onCleanup(@()deleteIfValid(app)); %#ok<NASGU>

app.runStudyForTest(fixture.cars,fixture.cases,@runRenderFixture);
verifyNotEmpty(testCase,findall(app.AeroLoadsAxes,"Type","line"));

app.ClearButtonPushed([],[]);

axesList = [app.BalanceAxes(:);app.AeroLoadsAxes(:); ...
    app.SuspensionAxes(:);app.RawRampAxes(:)];
for ax = axesList(:).'
    verifyEmpty(testCase,findall(ax,"Type","line"));
    verifyEmpty(testCase,findall(ax,"Type","constantline"));
end
end


function testProgressShowsSpeedContext(testCase)
ensureRampSpeedAppPath();
fixture = makeRampFixture();
app = RampSpeedApp("Visible","off");
cleanup = onCleanup(@()deleteIfValid(app)); %#ok<NASGU>

app.runStudyForTest(fixture.cars,fixture.cases,@runWithSpeedProgress);
messages = string(app.ProgressTextArea.Value);
verifyTrue(testCase,any(contains(messages,"20")));

    function run = runWithSpeedProgress(~,caseInfo,~,callbacks)
        callbacks.onProgress(struct("phase","speed","speedIndex",2, ...
            "speed_mps",20,"completedSpeeds",1,"requestedSpeeds",2));
        run = makeFixtureRun(caseInfo);
    end
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
