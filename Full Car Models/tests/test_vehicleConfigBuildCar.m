function tests = test_vehicleConfigBuildCar
%TEST_VEHICLECONFIGBUILDCAR Verify editable setup construction and shared assets.
tests = functiontests(localfunctions);
end

function testGeometryAndGearingAreRebuilt(testCase)
base = vehicleConfig.defaultSetup();
edit = base;
edit.parameters.wheelbase = 1.60;
edit.parameters.weight_dist = 0.55;
edit.parameters.wheel_radius = 0.205;
edit.parameters.final_drive = 3.2;
car = vehicleConfig.buildCar(edit,"lap");
verifyEqual(testCase,car.W_b,1.60,'AbsTol',1e-12);
verifyEqual(testCase,car.l_f,1.60*0.55,'AbsTol',1e-12);
verifyEqual(testCase,car.l_r,1.60*0.45,'AbsTol',1e-12);
verifyEqual(testCase,car.R,car.powertrain.wheel_radius);
pt = car.powertrain;
expected = pt.shift_point./pt.gears/pt.final_drive/ ...
    pt.primary_reduction*pi/30*car.R;
expected(end) = pt.redline/pt.gears(end)/pt.final_drive/ ...
    pt.primary_reduction*pi/30*car.R;
verifyEqual(testCase,pt.switch_gear_velocities,expected,'AbsTol',1e-12);
verifyEqual(testCase,base.parameters.weight_dist,0.512);
end

function testSpringEditUpdatesWheelRateWithoutRetuningQssSplit(testCase)
setup = vehicleConfig.defaultSetup();
setup.parameters.spring_rate_front_lb_in = 350;
assets = vehicleConfig.loadAssets(setup);
[car,derived] = vehicleConfig.buildCar(setup,"lap",assets);
verifyEqual(testCase,car.rideHeightAero.wheel_rate_front_Npm, ...
    350*175.126835*setup.parameters.motion_ratio_front^2,'AbsTol',1e-9);
verifyEqual(testCase,car.R_sf,setup.parameters.R_sf);
verifyEqual(testCase,derived.wheelRateFront_Npm, ...
    car.rideHeightAero.wheel_rate_front_Npm);
end

function testRequiredPhysicalBoundsRejectInvalidValues(testCase)
base = vehicleConfig.defaultSetup();
cases = { ...
    'wheelbase',0; ...
    'wheelbase',-1; ...
    'mass',0; ...
    'mass',-1; ...
    'track_width',0; ...
    'track_width',-1; ...
    'wheel_radius',0; ...
    'wheel_radius',-0.1; ...
    'weight_dist',0; ...
    'weight_dist',1; ...
    'weight_dist',-0.01; ...
    'weight_dist',1.01; ...
    'spring_rate_front_lb_in',0; ...
    'spring_rate_front_lb_in',-1; ...
    'spring_rate_rear_lb_in',0; ...
    'spring_rate_rear_lb_in',-1; ...
    'motion_ratio_front',0; ...
    'motion_ratio_front',-1; ...
    'motion_ratio_rear',0; ...
    'motion_ratio_rear',-1};
for k = 1:size(cases,1)
    setup = base;
    setup.parameters.(cases{k,1}) = cases{k,2};
    verifyError(testCase,@() vehicleConfig.buildCar(setup,"lap"), ...
        'vehicleConfig:invalidParameter');
end
end

function testDriverWeightsAreNonnegativeAndZeroRemainsValid(testCase)
setup = vehicleConfig.defaultSetup();
setup.parameters.driver_weight = -1;
verifyError(testCase,@() vehicleConfig.buildCar(setup,"lap"), ...
    'vehicleConfig:invalidParameter');
setup = vehicleConfig.defaultSetup();
setup.parameters.accel_driver_weight = -1;
verifyError(testCase,@() vehicleConfig.buildCar(setup,"accel"), ...
    'vehicleConfig:invalidParameter');

setup = vehicleConfig.defaultSetup();
setup.parameters.driver_weight = 0;
setup.parameters.accel_driver_weight = 0;
lap = vehicleConfig.buildCar(setup,"lap");
accel = vehicleConfig.buildCar(setup,"accel");
verifyEqual(testCase,lap.M,setup.parameters.mass,'AbsTol',1e-12);
verifyEqual(testCase,accel.M,setup.parameters.mass,'AbsTol',1e-12);
end

function testCanonicalSetupIdentityMustBeNonempty(testCase)
base = vehicleConfig.defaultSetup();
emptyValues = {''; ""};
identityFields = {'id','label'};
for fieldIndex = 1:numel(identityFields)
    for valueIndex = 1:numel(emptyValues)
        setup = base;
        setup.(identityFields{fieldIndex}) = emptyValues{valueIndex};
        verifyError(testCase,@() vehicleConfig.buildCar(setup,"lap"), ...
            'vehicleConfig:badSchema');
    end
end
end

function testNonScalarComplexAndNonFiniteSourcesAreRejected(testCase)
base = vehicleConfig.defaultSetup();
badValues = {[1 2],NaN,1+2i,"1.6"};
for k = 1:numel(badValues)
    setup = base;
    setup.parameters.wheelbase = badValues{k};
    verifyError(testCase,@() vehicleConfig.buildCar(setup,"lap"), ...
        'vehicleConfig:invalidParameter');
end
end

function testMalformedTorqueCurvesAreRejectedAtBuildTime(testCase)
base = vehicleConfig.defaultSetup();
badCurves = {42,ones(3,2),ones(2,1),[1 Inf;2 3], ...
    [1+1i 2;3 4],"42"};
for k = 1:numel(badCurves)
    setup = base;
    setup.fixed.torque_fn = badCurves{k};
    verifyError(testCase,@() vehicleConfig.buildCar(setup,"lap"), ...
        'vehicleConfig:invalidParameter');
end
end
function testUnknownFieldsRolesProfilesAndSchemasHaveSpecificErrors(testCase)
base = vehicleConfig.defaultSetup();
setup = base;
setup.parameters.futureValue = 1;
verifyError(testCase,@() vehicleConfig.buildCar(setup,"lap"), ...
    'vehicleConfig:unknownField');
setup = base;
setup.fixed.futureInput = 1;
verifyError(testCase,@() vehicleConfig.buildCar(setup,"lap"), ...
    'vehicleConfig:unknownField');
setup = base;
setup.futureSchemaKey = true;
verifyError(testCase,@() vehicleConfig.buildCar(setup,"lap"), ...
    'vehicleConfig:unknownField');
verifyError(testCase,@() vehicleConfig.buildCar(base,"ramp"), ...
    'vehicleConfig:badRole');
setup = base;
setup.profile = "unsupported-profile";
verifyError(testCase,@() vehicleConfig.buildCar(setup,"lap"), ...
    'vehicleConfig:badSchema');
setup = base;
setup.schemaVersion = 2;
verifyError(testCase,@() vehicleConfig.buildCar(setup,"lap"), ...
    'vehicleConfig:badSchema');
verifyError(testCase,@() vehicleConfig.buildCar(repmat(base,2,1),"lap"), ...
    'vehicleConfig:badSchema');
end

function testSignedCamberToeRideHeightsAndZeroTermsRemainValid(testCase)
setup = vehicleConfig.defaultSetup();
setup.parameters.gamma_f = -1.2;
setup.parameters.static_r_toe = -0.03;
setup.parameters.static_front_ride_height_in = -0.15;
setup.parameters.static_rear_ride_height_in = 0.08;
setup.parameters.map_reference_front_ride_height_in = -0.2;
setup.parameters.map_reference_rear_ride_height_in = 0.1;
setup.parameters.I_zz = 0;
setup.parameters.I_wheel = 0;
setup.parameters.I_driveline = 0;
setup.parameters.Crr = 0;
setup.parameters.camber_compliance_f = 0;
setup.parameters.camber_compliance_r = 0;
assets = vehicleConfig.loadAssets(setup);
car = vehicleConfig.buildCar(setup,"lap",assets);
verifyEqual(testCase,car.static_gamma_f,-1.2);
verifyEqual(testCase,car.static_r_toe,-0.03);
verifyEqual(testCase,car.rideHeightAero.static_front_ride_height_in,-0.15);
verifyEqual(testCase,car.rideHeightAero.static_rear_ride_height_in,0.08);
verifyEqual(testCase,car.rideHeightAero.map_reference_front_ride_height_in,-0.2);
verifyEqual(testCase,car.I_zz,0);
verifyEqual(testCase,car.I_wheel,0);
verifyEqual(testCase,car.I_driveline,0);
verifyEqual(testCase,car.Crr,0);
verifyEqual(testCase,car.camber_compliance_f,0);
verifyEqual(testCase,car.camber_compliance_r,0);
end

function testLapAndAccelerationRolesUseTheirOwnWeightsAndMapPolicy(testCase)
setup = vehicleConfig.defaultSetup();
assets = vehicleConfig.loadAssets(setup);
[lap,lapDerived] = vehicleConfig.buildCar(setup,"lap",assets);
[accel,accelDerived] = vehicleConfig.buildCar(setup,"accel",assets);
verifyEqual(testCase,lap.M,226,'AbsTol',1e-12);
verifyEqual(testCase,accel.M,221,'AbsTol',1e-12);
verifyEqual(testCase,lapDerived.driverWeight_kg,setup.parameters.driver_weight);
verifyEqual(testCase,accelDerived.driverWeight_kg,setup.parameters.accel_driver_weight);
verifyTrue(testCase,lap.aero.hasMap());
verifyFalse(testCase,accel.aero.hasMap());
verifyEqual(testCase,lap.R_sf,setup.parameters.R_sf);
verifyEqual(testCase,accel.R_sf,setup.parameters.R_sf);
verifyEqual(testCase,lap.aero.cda,setup.parameters.cda);
verifyEqual(testCase,lap.aero.cla,setup.parameters.cla);
verifyEqual(testCase,accel.aero.cda,setup.parameters.accel_cda);
verifyEqual(testCase,accel.aero.cla,setup.parameters.accel_cla);
end

function testDisabledMapDoesNotLoadOrAttachConfiguredPath(testCase)
setup = vehicleConfig.defaultSetup();
setup.aeroMap.enabled = false;
setup.aeroMap.path = fullfile(repositoryRoot(),'Full Car Models','does-not-exist.csv');
assets = vehicleConfig.loadAssets(setup);
verifyEmpty(testCase,assets.map);
[lap,~] = vehicleConfig.buildCar(setup,"lap",assets);
[accel,~] = vehicleConfig.buildCar(setup,"accel",assets);
verifyFalse(testCase,lap.aero.hasMap());
verifyFalse(testCase,accel.aero.hasMap());
verifyEqual(testCase,strlength(assets.fingerprints.aeroMap),0);
end

function testInjectedMapCorrectionsAndCamberAssetsArePreserved(testCase)
setup = vehicleConfig.defaultSetup();
assets = vehicleConfig.loadAssets(setup);
assets.map = assets.map.withCorrections(1.12,0.93,0.025);
assets.camberModelData.betaL = assets.camberModelData.betaL + 1e-4;
assets.camberRatioData.camber2ratio = ...
    assets.camberRatioData.camber2ratio + 1e-4;
expectedMap = assets.map;
expectedModelData = assets.camberModelData;
expectedRatioData = assets.camberRatioData;
car = vehicleConfig.buildCar(setup,"lap",assets);
verifyEqual(testCase,carConfigurationSnapshot(car.aero.map), ...
    carConfigurationSnapshot(expectedMap));
verifyEqual(testCase,car.camberKinematics.camberModelData,expectedModelData);
verifyEqual(testCase,car.tire.camber2ratio,expectedRatioData.camber2ratio);
verifyEqual(testCase,car.tire.camber4ratio,expectedRatioData.camber4ratio);
end

function testAssetContextCannotBeReusedForAnotherIdentity(testCase)
setup = vehicleConfig.defaultSetup();
assets = vehicleConfig.loadAssets(setup);
changed = setup;
changed.aeroMap.path = fullfile(repositoryRoot(),'Full Car Models', ...
    'different-map.csv');
verifyError(testCase,@() vehicleConfig.buildCar(changed,"lap",assets), ...
    'vehicleConfig:assetMismatch');
end

function testLoadedLc0ModelIsReusedAndKeepsUncertainty(testCase)
setup = vehicleConfig.defaultSetup("lc0_nd","high");
assets = vehicleConfig.loadAssets(setup);
assets.tireModel.task3Sentinel = "loaded-once";
[car,derived] = vehicleConfig.buildCar(setup,"accel",assets);
verifyClass(testCase,car.tire,'NondimensionalTire');
verifyEqual(testCase,car.tire.uncertainty_mode,"high");
verifyEqual(testCase,car.tire.model.task3Sentinel,"loaded-once");
verifyFalse(testCase,car.aero.hasMap());
verifyEqual(testCase,derived.assetFingerprints,assets.fingerprints);
end

function testRebuildStartsUnsolvedAndDoesNotMutateEarlierCar(testCase)
setup = vehicleConfig.defaultSetup();
assets = vehicleConfig.loadAssets(setup);
[first,~] = vehicleConfig.buildCar(setup,"lap",assets);
first.ggPoints = [1 2 3];
first.ss_info = ones(1,6);
first.longAccelLookup = [4 5];
first.longDecelLookup = [6 7];
first.comp = struct('sentinel',true);
firstBefore = carConfigurationSnapshot(first);
[rebuilt,~] = vehicleConfig.buildCar(setup,"lap",assets);
verifyEmpty(testCase,rebuilt.ggPoints);
verifyEmpty(testCase,rebuilt.ss_info);
verifyEmpty(testCase,rebuilt.accel_info);
verifyEmpty(testCase,rebuilt.decel_info);
verifyEmpty(testCase,rebuilt.longAccelLookup);
verifyEmpty(testCase,rebuilt.longDecelLookup);
verifyEmpty(testCase,rebuilt.comp);
verifyTrue(testCase,isstruct(rebuilt.ggMask));
verifyEmpty(testCase,fieldnames(rebuilt.ggMask));
[~,~] = vehicleConfig.buildCar(setup,"accel",assets);
verifyEqual(testCase,carConfigurationSnapshot(first),firstBefore);
end

function testQssRolesAndSelectionsMatchFrozenConstructionSnapshots(testCase)
reference = constructionReference(testCase);
base = vehicleConfig.defaultSetup();
assets = vehicleConfig.loadAssets(base);
[lap,~] = vehicleConfig.buildCar(base,"lap",assets);
[accel,~] = vehicleConfig.buildCar(base,"accel",assets);
verifyEqual(testCase,carConfigurationSnapshot(lap),reference.qss.default.lapCar);
verifyEqual(testCase,carConfigurationSnapshot(accel),reference.qss.default.accelCar);
verifyEqual(testCase,lap.M,226,'AbsTol',1e-12);
verifyEqual(testCase,accel.M,221,'AbsTol',1e-12);

explicit = reference.qss.explicit;
for row = 1:height(explicit.designTable)
    setup = applyDesignRow(base,explicit.designTable,row);
    [lap,~] = vehicleConfig.buildCar(setup,"lap",assets);
    [accel,~] = vehicleConfig.buildCar(setup,"accel",assets);
    verifyEqual(testCase,carConfigurationSnapshot(lap), ...
        explicit.lapCars{row});
    verifyEqual(testCase,carConfigurationSnapshot(accel), ...
        explicit.accelCars{row});
end

legacyAssets = vehicleConfig.loadAssets(vehicleConfig.defaultSetup("legacy","low"));
lc0Assets = vehicleConfig.loadAssets(vehicleConfig.defaultSetup("lc0_nd","low"));
selectionNames = fieldnames(reference.qss.tireSelections);
for k = 1:numel(selectionNames)
    item = reference.qss.tireSelections.(selectionNames{k});
    setup = vehicleConfig.defaultSetup(item.model,item.uncertainty);
    if string(item.model) == "legacy"
        selectionAssets = legacyAssets;
    else
        selectionAssets = lc0Assets;
    end
    [lap,~] = vehicleConfig.buildCar(setup,"lap",selectionAssets);
    [accel,~] = vehicleConfig.buildCar(setup,"accel",selectionAssets);
    verifyEqual(testCase,carConfigurationSnapshot(lap),item.lapCar);
    verifyEqual(testCase,carConfigurationSnapshot(accel),item.accelCar);
end
end

function testSavedValidStateHasMatchingDirectEquationResults(testCase)
reference = constructionReference(testCase);
numericalPath = fullfile(repositoryRoot(),'artifacts', ...
    'car_configuration_refactor','2026-10-04-reference-production-28mps', ...
    'numerical-reference.mat');
if strcmp(getenv('VD_TASK3_REQUIRE_NUMERICAL_REFERENCE'),'1') && ...
        ~isfile(numericalPath)
    error('test_vehicleConfigBuildCar:missingNumericalReference', ...
        'VD_TASK3_REQUIRE_NUMERICAL_REFERENCE is set but the numerical capture is missing.');
end
assumeTrue(testCase,isfile(numericalPath), ...
    'The Task 1 numerical reference is not present in this clone.');
numerical = load(numericalPath,'numerical').numerical;
gg = numerical.qss.ggParamSets;
validIndex = [];
for k = 1:numel(gg)
    candidate = gg(k);
    if candidate.maxLatValid && any(candidate.maxLatFlag == [1 2]) && ...
            numel(candidate.maxLatx) >= 12 && ...
            all(isfinite(candidate.maxLatx(4:12)))
        validIndex = k;
        break
    end
end
assumeTrue(testCase,~isempty(validIndex), ...
    'The Task 1 numerical reference contains no valid lateral state.');
candidate = gg(validIndex);
state = candidate.maxLatx(4:12);
verifySize(testCase,state,[1 9]);
verifyEqual(testCase,state(3),candidate.longVel,'AbsTol',1e-12);
verifyTrue(testCase,all(isfinite(state)));
verifyTrue(testCase,candidate.car.aero.hasMap());
verifyEqual(testCase,carConfigurationSnapshot(candidate.car), ...
    reference.qss.default.lapCar);

setup = vehicleConfig.defaultSetup();
assets = vehicleConfig.loadAssets(setup);
rebuilt = vehicleConfig.buildCar(setup,"lap",assets);
expected = equationOutputs(candidate.car,state);
actual = equationOutputs(rebuilt,state);
for k = 1:numel(expected)
    verifyEqual(testCase,isfinite(actual{k}),isfinite(expected{k}));
    verifyEqual(testCase,isnan(actual{k}),isnan(expected{k}));
    verifyEqual(testCase,isinf(actual{k}),isinf(expected{k}));
    finite = isfinite(expected{k});
    verifyEqual(testCase,actual{k}(finite),expected{k}(finite), ...
        'AbsTol',1e-12,'RelTol',1e-12);
end
end

function outputs = equationOutputs(car,state)
outputs = cell(1,15);
[outputs{:}] = car.equations(state);
end

function setup = applyDesignRow(setup,designTable,row)
names = designTable.Properties.VariableNames;
for k = 1:numel(names)
    setup.parameters.(names{k}) = designTable{row,names{k}};
end
end

function reference = constructionReference(testCase)
referencePath = fullfile(repositoryRoot(),'artifacts', ...
    'car_configuration_refactor','2026-10-04-reference', ...
    'construction-reference.mat');
if strcmp(getenv('VD_TASK3_REQUIRE_REFERENCE'),'1') && ~isfile(referencePath)
    error('test_vehicleConfigBuildCar:missingReference', ...
        'VD_TASK3_REQUIRE_REFERENCE is set but the construction snapshot is missing.');
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
