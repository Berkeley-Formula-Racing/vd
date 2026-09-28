function tests = test_rampSpeedCompactModel
% Contract tests for the serializable ramp-speed model/cache boundary.
tests = functiontests(localfunctions);
end

function testBuildProducesSerializableModelWithCachedEnvelope(testCase)
ensureRampSpeedPath();
[car,config] = carConfigBaseline();
profile = rampSpeed.resolveSolverProfile("fastPreview");
speeds = [5 10 20 25];

model = rampSpeed.buildRampModel(car,config.defaultSetup,profile,speeds);

verifyEqual(testCase,string(model.modelKind),"rampSpeedLite");
verifyEqual(testCase,string(model.modelVersion),"ramp-speed-lite-1");
verifyEqual(testCase,string(model.setupId),string(config.defaultSetup.id));
verifyEqual(testCase,string(model.baselineVersion), ...
    string(config.defaultSetup.baselineVersion));
verifyEqual(testCase,string(model.powertrainModel),"continuousEnvelope");
verifyEqual(testCase,model.speed_mps,speeds(:),"AbsTol",0);
verifySize(testCase,model.continuousEnvelope,[numel(speeds) 1]);
verifyFalse(testCase,isfield(model,"car"));
verifyFalse(testCase,contains(jsonencode(model),"function_handle"));
verifyGreaterThan(testCase,strlength(string(jsonencode(model))),0);
verifyEqual(testCase,model.vehicle.mass_kg,car.M,"AbsTol",0);
verifyEqual(testCase,model.driver.weight_kg, ...
    config.defaultSetup.driverWeight_kg,"AbsTol",0);
verifyEqual(testCase,model.suspension.R_sf,car.R_sf,"AbsTol",0);
end

function testEnvelopeCacheIsUsedByContinuousPointSolver(testCase)
ensureRampSpeedPath();
[car,config] = carConfigBaseline();
profile = rampSpeed.resolveSolverProfile("fastPreview",struct( ...
    "maxFunctionEvaluations",350));
model = rampSpeed.buildRampModel(car,config.defaultSetup,profile,20);

result = rampSpeed.solveLongitudinalPoint(car,20,[],profile, ...
    struct("rampModel",model));

verifyTrue(testCase,ismember(string(result.status), ...
    ["converged","near_feasible"]));
verifyEqual(testCase,string(result.diagnostics.envelopeSource), ...
    "cachedRampModel");
verifyEqual(testCase,result.diagnostics.continuousEnvelope.speed_mps, ...
    20,"AbsTol",0);
end

function testEnvelopeCompilationRejectsInvalidSpeed(testCase)
ensureRampSpeedPath();
[car,~] = carConfigBaseline();
verifyError(testCase,@()rampSpeed.compilePowertrainEnvelope(car,[10 0]), ...
    "rampSpeed:invalidSpeed");
end

function testModelRetainsMixedSpeedFailuresWithoutAborting(testCase)
ensureRampSpeedPath();
[car,config] = carConfigBaseline();
model = rampSpeed.buildRampModel(car,config.defaultSetup, ...
    rampSpeed.resolveSolverProfile("fastPreview"),[5 1e6]);

verifyEqual(testCase,model.speed_mps,[5;1e6],"AbsTol",0);
verifySize(testCase,model.continuousEnvelope,[1 1]);
verifySize(testCase,model.continuousEnvelopeErrors,[1 1]);
verifyEqual(testCase,model.continuousEnvelopeErrors.speed_mps,1e6, ...
    "AbsTol",0);
end

function testCanonicalLongitudinalRunStoresAndUsesRampModel(testCase)
ensureRampSpeedPath();
[car,config] = carConfigBaseline();
settings = struct("speeds",20,"solverProfile","fastPreview", ...
    "speedGrid",struct("mode","fixed"),"verbose",false);
caseInfo = struct("id",config.defaultSetup.id, ...
    "label",config.defaultSetup.label,"setupSpec",config.defaultSetup);

run = rampSpeed.runLongitudinalRamp(car,settings,caseInfo,struct());

verifyEqual(testCase,string(run.runMeta.rampModel.modelKind), ...
    "rampSpeedLite");
verifyEqual(testCase,string(run.runMeta.rampModel.powertrainModel), ...
    "continuousEnvelope");
verifyEqual(testCase,run.runMeta.rampModel.speed_mps,20,"AbsTol",0);
verifyEqual(testCase,string(run.raw.diagnostics(1).diagnostics.envelopeSource), ...
    "cachedRampModel");
end

function testDuplicatedSetupModelsKeepDistinctProvenance(testCase)
ensureRampSpeedPath();
[~,config] = carConfigBaseline();
copy = rampSpeed.duplicateSetup(config.defaultSetup,"rear-high","Rear ARB high");
copy.rearArbStiffness_NmPerRad = config.options.rearArbStiffness_NmPerRad(end);
[cars,cases] = rampSpeed.buildSetupCatalog(config,[config.defaultSetup,copy]);
profile = rampSpeed.resolveSolverProfile("fastPreview");
first = rampSpeed.buildRampModel(cars{1},cases(1).setupSpec,profile,20);
second = rampSpeed.buildRampModel(cars{2},cases(2).setupSpec,profile,20);

verifyNotEqual(testCase,string(first.setupId),string(second.setupId));
verifyNotEqual(testCase, ...
    first.suspension.rearArbStiffness_NmPerRad, ...
    second.suspension.rearArbStiffness_NmPerRad);
verifyNotEqual(testCase,first.suspension.R_sf,second.suspension.R_sf);
verifyEqual(testCase,string(first.aero.aeroMapId),string(second.aero.aeroMapId));
verifyEqual(testCase,string(first.powertrain.powertrainModel), ...
    string(second.powertrain.powertrainModel));
end

function ensureRampSpeedPath()
root = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(root));
end
