function tests = test_rampSpeedSetupEditorApp
tests = functiontests(localfunctions);
end

function testAppStartsWithRealBaselineAndExplicitSetupColumns(testCase)
ensureRampSpeedAppPath();
app = RampSpeedApp("Visible","off");
cleanup = onCleanup(@()deleteIfValid(app));
verifySize(testCase,app.Cars,[1 1]);
verifyClass(testCase,app.Cars{1},'Car');
verifyEqual(testCase,string(app.CarRoleDropDown.Visible),"off");
verifyTrue(testCase,app.Cars{1}.g > 0);
verifyEqual(testCase,string(app.SetupTable.ColumnName(:)), ...
    ["ID","Label","Rear ARB (N m/rad)","Front spring (lb/in)", ...
    "Rear spring (lb/in)","Front ride height (in)","Rear ride height (in)", ...
    "Driver weight (kg)","Rear weight dist. (%)","Aero map","R_sf","Source"]');
verifyEqual(testCase,app.Cases(1).setupSpec.frontRideHeight_in, ...
    3.86558910640757, ...
    'AbsTol',1e-12);
verifyEqual(testCase,app.Cases(1).setupSpec.rearRideHeight_in, ...
    5.76907823448674, ...
    'AbsTol',1e-12);
verifyEqual(testCase,app.FrontRideHeightEditField.Value, ...
    3.86558910640757,'AbsTol',1e-12);
verifyEqual(testCase,app.RearRideHeightEditField.Value, ...
    5.76907823448674,'AbsTol',1e-12);
verifyTrue(testCase,all(ismember( ...
    {'rearArbStiffness_NmPerRad','frontSpringRate_lb_in', ...
    'rearSpringRate_lb_in','frontRideHeight_in','rearRideHeight_in', ...
    'driverWeight_kg','rearWeightDistribution_percent','aeroMapId','R_sf'},app.SetupTable.Data.Properties.VariableNames)));
end

function testMultiRowTableSelectionRetainsAllOverlayCases(testCase)
ensureRampSpeedAppPath();
app = RampSpeedApp("Visible","off");
cleanup = onCleanup(@()deleteIfValid(app));
app.duplicateSelectedSetup("overlay-copy","Overlay copy");
ids = string({app.Cases.id}).';

% App Designer can report only the last clicked cell while the table still
% owns the complete multi-selection. The app must use that complete state.
app.SetupTable.Selection = [1 1; 2 1];
app.SetupTable.CellSelectionCallback(app.SetupTable, ...
    struct("Indices",[2 1]));

verifyEqual(testCase,app.Session.viewModel().selectedCaseIds,ids);
end

function testDuplicateEditAndDeletePreserveBaseline(testCase)
ensureRampSpeedAppPath();
app = RampSpeedApp("Visible","off");
cleanup = onCleanup(@()deleteIfValid(app));
baselineRsf = app.Cars{1}.R_sf;
copy = app.duplicateSelectedSetup("rear-high","Rear ARB high");
verifyEqual(testCase,string(copy.id),"rear-high");
verifySize(testCase,app.Cars,[2 1]);
verifyEqual(testCase,string(app.Cases(1).id),"ramp-baseline");
verifyTrue(testCase,app.Cases(1).isBaseline);
verifyFalse(testCase,app.Cases(2).isBaseline);

edited = app.Cases(2).setupSpec;
edited.rearArbStiffness_NmPerRad = ...
    app.BaselineConfig.options.rearArbStiffness_NmPerRad(end);
edited.driverWeight_kg = 80;
edited.rearWeightDistribution_percent = 54;
app.editSelectedSetup(edited);
verifyLessThan(testCase,app.Cars{2}.R_sf,baselineRsf);
verifyEqual(testCase,app.Cars{2}.M,242,'AbsTol',1e-12);
verifyEqual(testCase,app.Cars{2}.l_f/app.Cars{2}.W_b,0.54,'AbsTol',1e-12);
verifyEqual(testCase,app.Cars{1}.R_sf,baselineRsf,'AbsTol',1e-12);

app.deleteSelectedSetup();
verifySize(testCase,app.Cars,[1 1]);
verifyEqual(testCase,string(app.Cases(1).id),"ramp-baseline");
end

function testBaselineCannotBeDeleted(testCase)
ensureRampSpeedAppPath();
app = RampSpeedApp("Visible","off");
cleanup = onCleanup(@()deleteIfValid(app));
verifyError(testCase,@()app.deleteSelectedSetup(), ...
    'rampSpeed:baselineProtected');
end

function testRampTypeDoesNotSwitchSharedCar(testCase)
ensureRampSpeedAppPath();
app = RampSpeedApp("Visible","off");
cleanup = onCleanup(@()deleteIfValid(app));
app.CarRoleDropDown.Value = "lap";
app.RampTypeDropDown.Value = "longitudinal";
app.RampTypeDropDown.ValueChangedFcn(app.RampTypeDropDown,[]);
verifyEqual(testCase,string(app.CarRoleDropDown.Value),"lap");
verifySize(testCase,app.Cars,[1 1]);
end

function testReadOnlyLoadedStudyCannotRun(testCase)
ensureRampSpeedAppPath();
app = RampSpeedApp("Visible","off");
cleanup = onCleanup(@()deleteIfValid(app));
app.SetupReadOnly = true;
app.ClearButtonPushed([],[]);
verifyEqual(testCase,string(app.RunButton.Enable),"off");
app.RunButtonPushed([],[]);
verifyTrue(testCase,any(contains(string(app.ProgressTextArea.Value),"read-only")));
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
