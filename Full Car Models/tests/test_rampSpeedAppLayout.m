function tests = test_rampSpeedAppLayout
tests = functiontests(localfunctions);
end

function testResponsiveLayoutKeepsPlotsAndSetupReachable(testCase)
ensureRampSpeedAppPath();
app = RampSpeedApp("Visible","off");
cleanup = onCleanup(@()deleteIfValid(app)); %#ok<NASGU>

verifyEqual(testCase,string(app.SetupPanel.Scrollable),"on");
drawnow;
pause(1);
drawnow;
plotAxes = [app.BalanceAxes(:);app.AeroLoadsAxes(:); ...
    app.SuspensionAxes(:);app.RawRampAxes(:)];
windowSizes = [900 600; 1200 800; 700 450];
for sizeIndex = 1:size(windowSizes,1)
    app.UIFigure.Position = [100 100 windowSizes(sizeIndex,:)];
    drawnow;
    pause(2);
    drawnow;

    verifyGreaterThan(testCase,app.MainTabGroup.Position(3),1);
    verifyGreaterThan(testCase,app.MainTabGroup.Position(4),1);
    for axisIndex = 1:numel(plotAxes)
        plotAxis = plotAxes(axisIndex);
        parentPosition = plotAxis.Parent.Position;
        verifyLessThanOrEqual(testCase,plotAxis.Position(3), ...
            parentPosition(3));
        verifyLessThanOrEqual(testCase,plotAxis.Position(4), ...
            parentPosition(4));
    end
end
end

function testSetupControlsUseRunSetupTabAndCompactOverlayLegend(testCase)
ensureRampSpeedAppPath();
app = RampSpeedApp("Visible","off");
cleanup = onCleanup(@()deleteIfValid(app)); %#ok<NASGU>

verifyEqual(testCase,string(app.SetupTab.Title),"Setup");
verifyTrue(testCase,isvalid(app.SetupLegendPanel));
verifyEqual(testCase,string(app.SetupLegendPanel.Visible),"on");
verifyEqual(testCase,string(app.SetupLegendTable.ColumnName(:)), ...
    ["Show";"Setup"]);
verifyTrue(testCase,logical(app.SetupLegendTable.ColumnEditable(1)));

app.MainTabGroup.SelectedTab = app.SetupTab;
app.UIFigure.Position = [100 100 1200 800];
drawnow;
pause(0.2);
drawnow;

verifyEqual(testCase,string(app.SetupTable.Visible),"on");
verifyGreaterThanOrEqual(testCase,app.SetupTable.Position(4),100);

app.SetupLegendToggleButton.ButtonPushedFcn( ...
    app.SetupLegendToggleButton,[]);
verifyEqual(testCase,string(app.SetupLegendPanel.Visible),"off");
verifyEqual(testCase,app.RootGridLayout.ColumnWidth{1},0);

app.SetupLegendToggleButton.ButtonPushedFcn( ...
    app.SetupLegendToggleButton,[]);
verifyEqual(testCase,string(app.SetupLegendPanel.Visible),"on");
end

function testRunConfigurationLivesInSidebarAndSetupTabIsSetupOnly(testCase)
ensureRampSpeedAppPath();
app = RampSpeedApp("Visible","off");
cleanup = onCleanup(@()deleteIfValid(app)); %#ok<NASGU>

verifyEqual(testCase,string(app.SetupTab.Title),"Setup");
verifyTrue(testCase,isvalid(app.RunConfigPanel));
verifyTrue(testCase,isdescendantOf(app.RunConfigPanel,app.SetupLegendPanel));
verifyTrue(testCase,isdescendantOf(app.SetupTable,app.SetupPanel));

globalControls = [app.RampTypeDropDown;app.LateralModeDropDown; ...
    app.SpeedStartEditField;app.SpeedStopEditField;app.SpeedStepEditField; ...
    app.SpeedVectorEditField;app.NRampEditField;app.NBisectEditField; ...
    app.ResidualToleranceEditField;app.ExecutionModeDropDown; ...
    app.WorkersEditField;app.UnitDropDown;app.SolverProfileDropDown; ...
    app.SpeedGridModeDropDown;app.BaselineDropDown;app.RunButton; ...
    app.CancelButton;app.LoadButton;app.SaveButton;app.ExportButton; ...
    app.ClearButton;app.ProgressTextArea];
for control = globalControls(:).'
    verifyTrue(testCase,isdescendantOf(control,app.RunConfigPanel));
    verifyFalse(testCase,isdescendantOf(control,app.SetupPanel));
end
end

function testOverlayLegendControlsSelectedSetups(testCase)
ensureRampSpeedAppPath();
app = RampSpeedApp("Visible","off");
cleanup = onCleanup(@()deleteIfValid(app)); %#ok<NASGU>

app.duplicateSelectedSetup("legend-copy","Legend copy");
data = app.SetupLegendTable.Data;
verifyEqual(testCase,height(data),2);
data.show = [true;false];
app.SetupLegendTable.Data = data;
app.SetupLegendTable.CellEditCallback(app.SetupLegendTable, ...
    struct("Indices",[2 1]));

verifyEqual(testCase,app.Session.viewModel().selectedCaseIds, ...
    string(app.Cases(1).id));
data = app.SetupLegendTable.Data;
verifyTrue(testCase,data.show(1));
verifyFalse(testCase,data.show(2));
end

function ensureRampSpeedAppPath()
root = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(root));
end

function result = isdescendantOf(component,ancestor)
result = false;
parent = component.Parent;
while ~isempty(parent)
    if isequal(parent,ancestor)
        result = true;
        return;
    end
    parent = parent.Parent;
end
end

function deleteIfValid(app)
if ~isempty(app) && isvalid(app)
    delete(app);
end
end
