function tests = test_vehicleConfigRampCompatibility
tests = functiontests(localfunctions);
end

function testRampAdapterMatchesFrozenConstructionSnapshot(testCase)
[baseline,config] = carConfigBaseline();
setup = config.defaultSetup;
assets = rampSpeed.loadRampAssets(config,setup);
reference = constructionReference(testCase);
spec = vehicleConfig.fromRampSetup(config.fixedParameters,setup, ...
    assets.map,assets);
verifyEqual(testCase,string(spec.profile),"ramp");
verifyEqual(testCase,string(spec.suspension.rollSplitPolicy),"derived");
verifyFalse(testCase,isfield(spec.parameters,'driver_weight'));
verifyFalse(testCase,isfield(spec.parameters,'cda'));
verifyFalse(testCase,isfield(spec.parameters,'cla'));
verifyFalse(testCase,isfield(spec.parameters,'R_sf'));
[actual,derived] = vehicleConfig.buildCar(spec,"accel",assets);
verifyEqual(testCase,carConfigurationSnapshot(actual),reference.ramp.snapshot);
verifyEqual(testCase,carConfigurationSnapshot(baseline), ...
    reference.ramp.snapshot);
verifyEqual(testCase,derived,reference.ramp.derived);
end

function testRampProfileRejectsLapRole(testCase)
[~,config] = carConfigBaseline();
setup = config.defaultSetup;
assets = rampSpeed.loadRampAssets(config,setup);
spec = vehicleConfig.fromRampSetup(config.fixedParameters,setup, ...
    assets.map,assets);
verifyError(testCase,@() vehicleConfig.buildCar(spec,"lap",assets), ...
    'vehicleConfig:badRole');
end

function testRampProfileRejectsQssOnlyParameterFields(testCase)
[~,config] = carConfigBaseline();
setup = config.defaultSetup;
assets = rampSpeed.loadRampAssets(config,setup);
spec = vehicleConfig.fromRampSetup(config.fixedParameters,setup, ...
    assets.map,assets);
spec.parameters.driver_weight = 64;
verifyError(testCase,@() vehicleConfig.buildCar(spec,"accel",assets), ...
    'vehicleConfig:unknownField');
end

function testRampTrackWidthMustBePositive(testCase)
[~,config] = carConfigBaseline();
source = config.defaultSetup;
for invalidWidth = [0 -1]
    setup = vehicleConfig.fromRampSetup(config.fixedParameters,source);
    setup.parameters.track_width = invalidWidth;
    verifyError(testCase,@() vehicleConfig.buildCar(setup,"accel"), ...
        'vehicleConfig:invalidParameter');
end
end

function testRampCanonicalIdentityMustBeNonempty(testCase)
[~,config] = carConfigBaseline();
source = config.defaultSetup;
setup = vehicleConfig.fromRampSetup(config.fixedParameters,source);
emptyValues = {''; ""};
identityFields = {'id','label'};
for fieldIndex = 1:numel(identityFields)
    for valueIndex = 1:numel(emptyValues)
        invalidSetup = setup;
        invalidSetup.(identityFields{fieldIndex}) = emptyValues{valueIndex};
        verifyError(testCase,@() vehicleConfig.buildCar(invalidSetup,"accel"), ...
            'vehicleConfig:badSchema');
    end
end
end

function testRampAdapterSuppliesDeterministicLegacyIdentity(testCase)
[~,config] = carConfigBaseline();
source = rmfield(config.defaultSetup,{'id','label'});
setup = vehicleConfig.fromRampSetup(config.fixedParameters,source);
verifyEqual(testCase,string(setup.id),"legacy-ramp");
verifyEqual(testCase,string(setup.label),"Legacy ramp setup");
car = vehicleConfig.buildCar(setup,"accel");
verifyGreaterThan(testCase,car.M,0);

source.id = "";
source.label = '';
setup = vehicleConfig.fromRampSetup(config.fixedParameters,source);
verifyEqual(testCase,string(setup.id),"legacy-ramp");
verifyEqual(testCase,string(setup.label),"Legacy ramp setup");
end

function testRampEditsRebuildDerivedRollRatesAndTransientFields(testCase)
[~,config] = carConfigBaseline();
setup = config.defaultSetup;
setup.frontSpringRate_lb_in = 350;
setup.rearSpringRate_lb_in = 300;
setup.rearArbStiffness_NmPerRad = 1493;
setup.driverWeight_kg = 70;
setup.rearWeightDistribution_percent = 54.2;
setup.frontRideHeight_in = setup.frontRideHeight_in + 0.05;
setup.rearRideHeight_in = setup.rearRideHeight_in - 0.04;
assets = rampSpeed.loadRampAssets(config,setup);
[car,derived] = buildSingleRampCar(config.fixedParameters,setup, ...
    assets.map,assets);
expected = rampSpeed.deriveRollStiffness(setup.frontSpringRate_lb_in, ...
    setup.rearSpringRate_lb_in,config.fixedParameters.vehicle.motion_ratio_front, ...
    config.fixedParameters.vehicle.motion_ratio_rear, ...
    config.fixedParameters.vehicle.track_width, ...
    config.fixedParameters.frontArbStiffness_NmPerRad, ...
    setup.rearArbStiffness_NmPerRad);

verifyEqual(testCase,car.M,232,'AbsTol',1e-12);
verifyEqual(testCase,car.R_sf,expected.R_sf,'AbsTol',1e-12);
verifyEqual(testCase,derived.R_sf,expected.R_sf,'AbsTol',1e-12);
verifyEqual(testCase,derived.wheelRateFront_Npm, ...
    expected.wheelRateFront_Npm,'AbsTol',1e-12);
verifyEqual(testCase,derived.wheelRateRear_Npm, ...
    expected.wheelRateRear_Npm,'AbsTol',1e-12);
verifyEqual(testCase,derived.frontWeightDistribution_percent,45.8, ...
    'AbsTol',1e-12);
verifyEqual(testCase,derived.rearWeightDistribution_percent,54.2, ...
    'AbsTol',1e-12);
verifyEqual(testCase,derived.driverWeight_kg,70,'AbsTol',1e-12);
verifyEqual(testCase,derived.staticRearLoad_N,car.M*car.g*0.542, ...
    'AbsTol',1e-12);
verifyEqual(testCase,derived.staticFrontLoad_N,car.M*car.g*0.458, ...
    'AbsTol',1e-12);
verifyEqual(testCase,car.rideHeightAero.static_front_ride_height_in, ...
    setup.frontRideHeight_in,'AbsTol',1e-12);
verifyEqual(testCase,car.rideHeightAero.static_rear_ride_height_in, ...
    setup.rearRideHeight_in,'AbsTol',1e-12);
verifyEqual(testCase,car.rideHeightAero.wheel_rate_front_Npm, ...
    expected.wheelRateFront_Npm,'AbsTol',1e-12);
verifyEqual(testCase,car.rideHeightAero.wheel_rate_rear_Npm, ...
    expected.wheelRateRear_Npm,'AbsTol',1e-12);
verifyEqual(testCase,car.k_rf,0,'AbsTol',1e-12);
verifyEqual(testCase,car.k_rr,setup.rearArbStiffness_NmPerRad,'AbsTol',1e-12);
verifyEqual(testCase,car.k_f_b,expected.wheelRateFront_Npm,'AbsTol',1e-12);
verifyEqual(testCase,car.k_r_b,expected.wheelRateRear_Npm,'AbsTol',1e-12);
verifyEqual(testCase,car.k_f_r,expected.wheelRateFront_Npm,'AbsTol',1e-12);
verifyEqual(testCase,car.k_r_r,expected.wheelRateRear_Npm,'AbsTol',1e-12);
verifyEqual(testCase,car.rs_total, ...
    (expected.frontRollStiffness_NmPerRad + ...
    expected.rearRollStiffness_NmPerRad)*pi/180,'AbsTol',1e-12);
verifyEqual(testCase,fieldnames(derived),{ ...
    'wheelRateFront_Npm';'wheelRateRear_Npm'; ...
    'springRollStiffnessFront_NmPerRad'; ...
    'springRollStiffnessRear_NmPerRad'; ...
    'frontRollStiffness_NmPerRad';'rearRollStiffness_NmPerRad';'R_sf'; ...
    'driverWeight_kg';'rearWeightDistribution_percent'; ...
    'frontWeightDistribution_percent';'totalMass_kg'; ...
    'staticRearLoad_N';'staticFrontLoad_N';'assetFingerprints'});
verifyEqual(testCase,derived.assetFingerprints,assets.fingerprints);
end

function testSelectedCorrectedMapAndStaticMapPolicy(testCase)
[~,config] = carConfigBaseline();
setup = config.defaultSetup;
assets = rampSpeed.loadRampAssets(config,setup);
selectedMap = assets.map.withCorrections(1.12,0.93,0.025);
[car,derived] = buildSingleRampCar(config.fixedParameters,setup, ...
    selectedMap,assets);
verifyTrue(testCase,car.aero.hasMap());
verifyEqual(testCase,carConfigurationSnapshot(car.aero.map), ...
    carConfigurationSnapshot(selectedMap));
verifyEqual(testCase,string(car.aero.map.sourcePath), ...
    string(selectedMap.sourcePath));
verifyEqual(testCase,car.aero.map.claScale,1.12,'AbsTol',1e-12);
verifyEqual(testCase,car.aero.map.cdaScale,0.93,'AbsTol',1e-12);
verifyEqual(testCase,car.aero.map.copOffset,0.025,'AbsTol',1e-12);
verifyEqual(testCase,car.rideHeightAero.map_reference_front_ride_height_in, ...
    config.fixedParameters.vehicle.map_reference_front_ride_height_in, ...
    'AbsTol',1e-12);
verifyEqual(testCase,car.rideHeightAero.map_reference_rear_ride_height_in, ...
    config.fixedParameters.vehicle.map_reference_rear_ride_height_in, ...
    'AbsTol',1e-12);
verifyEqual(testCase,derived.assetFingerprints,assets.fingerprints);

[staticCar,staticDerived] = buildSingleRampCar( ...
    config.fixedParameters,setup,[],assets);
verifyFalse(testCase,staticCar.aero.hasMap());
verifyEqual(testCase,staticDerived.assetFingerprints,assets.fingerprints);
[noAssetCar,noAssetDerived] = buildSingleRampCar( ...
    config.fixedParameters,setup);
verifyFalse(testCase,noAssetCar.aero.hasMap());
verifyFalse(testCase,isfield(noAssetDerived,'assetFingerprints'));
end

function testRampBaselineIsIndependentOfQssSetupChanges(testCase)
[baseline,config] = carConfigBaseline();
setup = config.defaultSetup;
rampAssets = rampSpeed.loadRampAssets(config,setup);
spec = vehicleConfig.fromRampSetup(config.fixedParameters,setup, ...
    rampAssets.map,rampAssets);
qss = vehicleConfig.defaultSetup();
qss.parameters.mass = qss.parameters.mass + 20;
qss.parameters.accel_cla = qss.parameters.accel_cla + 0.4;
qss.aeroMap.enabled = false;
qssAssets = vehicleConfig.loadAssets(qss);
vehicleConfig.buildCar(qss,"accel",qssAssets);
rebuilt = vehicleConfig.buildCar(spec,"accel",rampAssets);
verifyEqual(testCase,carConfigurationSnapshot(baseline), ...
    carConfigurationSnapshot(rebuilt));
end

function reference = constructionReference(testCase)
repo = repositoryRoot();
referencePath = fullfile(repo,'artifacts','car_configuration_refactor', ...
    '2026-10-04-reference','construction-reference.mat');
if strcmp(getenv('VD_TASK5_REQUIRE_REFERENCE'),'1') && ~isfile(referencePath)
    error('test_vehicleConfigRampCompatibility:missingReference', ...
        'VD_TASK5_REQUIRE_REFERENCE is set but the construction snapshot is missing.');
end
assumeTrue(testCase,isfile(referencePath), ...
    'The Task 1 construction reference is not present in this clone.');
loaded = load(referencePath,'reference');
reference = loaded.reference;
end

function root = repositoryRoot()
modelRoot = fileparts(fileparts(which('vehicleConfig.defaultSetup')));
root = fileparts(modelRoot);
end
