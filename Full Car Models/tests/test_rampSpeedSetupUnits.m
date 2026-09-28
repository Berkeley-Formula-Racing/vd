function tests = test_rampSpeedSetupUnits
tests = functiontests(localfunctions);
end

function testRideHeightEditorAlwaysUsesInches(testCase)
root = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(root));
app = RampSpeedApp("Visible","off");
cleanup = onCleanup(@()deleteIfValid(app));
app.duplicateSelectedSetup("unit-test","Unit test");

verifyEqual(testCase,string(app.UnitDropDown.Value),"SI");
verifyEqual(testCase,app.FrontRideHeightEditField.Value,4,"AbsTol",1e-12);
verifyEqual(testCase,app.RearRideHeightEditField.Value,5.7,"AbsTol",1e-12);
app.FrontRideHeightEditField.Value = 2;
app.FrontRideHeightEditField.ValueChangedFcn( ...
    app.FrontRideHeightEditField,[]);
verifyEqual(testCase,app.Cases(2).setupSpec.frontRideHeight_in,2, ...
    'AbsTol',1e-12);

app.UnitDropDown.Value = "imperial";
app.UnitDropDown.ValueChangedFcn(app.UnitDropDown,[]);
verifyEqual(testCase,app.FrontRideHeightEditField.Value,2,'AbsTol',1e-12);
app.UnitDropDown.Value = "SI";
app.UnitDropDown.ValueChangedFcn(app.UnitDropDown,[]);
verifyEqual(testCase,app.FrontRideHeightEditField.Value,2,'AbsTol',1e-12);
verifyEqual(testCase,app.Cases(2).setupSpec.frontRideHeight_in,2, ...
    'AbsTol',1e-12);
end

function deleteIfValid(app)
if ~isempty(app) && isvalid(app)
    delete(app);
end
end
