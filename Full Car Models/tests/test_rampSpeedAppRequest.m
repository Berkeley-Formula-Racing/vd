function tests = test_rampSpeedAppRequest
tests = functiontests(localfunctions);
end

function testDefaultRequestUsesLateralCoastAndSI(testCase)
ensureRampSpeedAppPath();
app = RampSpeedApp("Visible","off");
cleanup = onCleanup(@()deleteIfValid(app));

request = app.buildRequestForTest();

verifyEqual(testCase,string(request.rampType),"lateral");
verifyEqual(testCase,string(request.settings.mode),"coast");
verifyEqual(testCase,string(request.carRole),"auto");
verifyEqual(testCase,string(request.displayUnits.speed),"m/s");
verifyGreaterThan(testCase,min(request.settings.speeds),0);
verifyGreaterThanOrEqual(testCase,request.settings.nRamp,2);
verifyGreaterThanOrEqual(testCase,request.settings.nBisect,0);
verifyGreaterThan(testCase,request.settings.residualTolerance,0);
end

function testLongitudinalRoleAndExplicitVectorPropagate(testCase)
ensureRampSpeedAppPath();
app = RampSpeedApp("Visible","off");
cleanup = onCleanup(@()deleteIfValid(app));

app.RampTypeDropDown.Value = "longitudinal";
app.CarRoleDropDown.Value = "acceleration";
app.SpeedVectorEditField.Value = "5, 7.5, 10";
app.NRampEditField.Value = 4;
app.NBisectEditField.Value = 3;
app.ResidualToleranceEditField.Value = 1e-7;

request = app.buildRequestForTest();

verifyEqual(testCase,string(request.rampType),"longitudinal");
verifyEqual(testCase,string(request.carRole),"auto");
verifyEqual(testCase,request.settings.speeds,[5 7.5 10],"AbsTol",0);
verifyEqual(testCase,request.settings.nRamp,4);
verifyEqual(testCase,request.settings.nBisect,3);
verifyEqual(testCase,request.settings.residualTolerance,1e-7,"AbsTol",0);
end

function testRequestValidationRejectsNonPhysicalValues(testCase)
ensureRampSpeedAppPath();
app = RampSpeedApp("Visible","off");
cleanup = onCleanup(@()deleteIfValid(app));

app.SpeedStartEditField.Value = 0;
verifyError(testCase,@()app.buildRequestForTest(), ...
    "rampSpeed:app:invalidRequest");
app.SpeedStartEditField.Value = 5;

app.NRampEditField.Value = 1;
verifyError(testCase,@()app.buildRequestForTest(), ...
    "rampSpeed:app:invalidRequest");
app.NRampEditField.Value = 3;

app.NBisectEditField.Value = -1;
verifyError(testCase,@()app.buildRequestForTest(), ...
    "rampSpeed:app:invalidRequest");
app.NBisectEditField.Value = 2;

app.ResidualToleranceEditField.Value = 0;
verifyError(testCase,@()app.buildRequestForTest(), ...
    "rampSpeed:app:invalidRequest");
app.ResidualToleranceEditField.Value = 1e-6;

app.SpeedVectorEditField.Value = "5, nope, 10";
verifyError(testCase,@()app.buildRequestForTest(), ...
    "rampSpeed:app:invalidRequest");
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
