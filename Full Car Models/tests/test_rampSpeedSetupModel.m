function tests = test_rampSpeedSetupModel
tests = functiontests(localfunctions);
end

function testRollStiffnessMatchesSpecifiedFormula(testCase)
frontSpring = 300;
rearSpring = 250;
mrFront = 0.847;
mrRear = 0.984;
track = 47*0.0254;
frontArb = 0;
rearArb = 986;
derived = rampSpeed.deriveRollStiffness(frontSpring,rearSpring, ...
    mrFront,mrRear,track,frontArb,rearArb);
lbInToNpm = 175.126835;
frontWheel = frontSpring*lbInToNpm*mrFront^2;
rearWheel = rearSpring*lbInToNpm*mrRear^2;
frontTotal = frontWheel*track^2/2 + frontArb;
rearTotal = rearWheel*track^2/2 + rearArb;
verifyEqual(testCase,derived.R_sf,frontTotal/(frontTotal+rearTotal), ...
    'AbsTol',1e-12);
verifyEqual(testCase,derived.frontRollStiffness_NmPerRad,frontTotal, ...
    'AbsTol',1e-12);
verifyEqual(testCase,derived.rearRollStiffness_NmPerRad,rearTotal, ...
    'AbsTol',1e-12);
end

function testIncreasingRearArbReducesFrontRollSplit(testCase)
[~,config] = carConfigBaseline();
base = config.defaultSetup;
[~,~,low] = rampSpeed.buildCarFromSetup(config, ...
    setField(base,'rearArbStiffness_NmPerRad',config.options.rearArbStiffness_NmPerRad(2)));
[~,~,high] = rampSpeed.buildCarFromSetup(config, ...
    setField(base,'rearArbStiffness_NmPerRad',config.options.rearArbStiffness_NmPerRad(end)));
verifyLessThan(testCase,high.R_sf,low.R_sf);
end

function testDuplicateDoesNotMutateProtectedBaseline(testCase)
[~,config] = carConfigBaseline();
baseline = config.defaultSetup;
copy = rampSpeed.duplicateSetup(baseline,"rear-arb-high","Rear ARB high");
copy.rearArbStiffness_NmPerRad = config.options.rearArbStiffness_NmPerRad(end);
[car,normalized,derived] = rampSpeed.buildCarFromSetup(config,copy);
verifyEqual(testCase,baseline.rearArbStiffness_NmPerRad,986);
verifyEqual(testCase,config.defaultSetup.rearArbStiffness_NmPerRad,986);
verifyEqual(testCase,normalized.id,"rear-arb-high");
verifyEqual(testCase,normalized.label,"Rear ARB high");
verifyFalse(testCase,normalized.isBaseline);
verifyClass(testCase,car,'Car');
verifyEqual(testCase,car.R_sf,derived.R_sf,'AbsTol',1e-12);
end

function testOnlyBaselineDiscreteChoicesAreAccepted(testCase)
[~,config] = carConfigBaseline();
bad = config.defaultSetup;
bad.frontSpringRate_lb_in = 301;
verifyError(testCase,@()rampSpeed.buildCarFromSetup(config,bad), ...
    'rampSpeed:invalidSetupOption');
bad = config.defaultSetup;
bad.rearArbStiffness_NmPerRad = 1;
verifyError(testCase,@()rampSpeed.buildCarFromSetup(config,bad), ...
    'rampSpeed:invalidSetupOption');
end

function testRideHeightsReachCarAndWarnOutsideMapEnvelope(testCase)
[~,config] = carConfigBaseline();
setup = config.defaultSetup;
setup.frontRideHeight_in = 100;
setup.rearRideHeight_in = -100;
[car,normalized,derived] = rampSpeed.buildCarFromSetup(config,setup);
verifyEqual(testCase,car.rideHeightAero.static_front_ride_height_in,100);
verifyEqual(testCase,car.rideHeightAero.static_rear_ride_height_in,-100);
verifyEqual(testCase,normalized.aeroMapId,"b26");
verifyTrue(testCase,derived.aeroMapOutsideEnvelope);
end

function testAeroMapCatalogAndMalformedMapValidation(testCase)
catalog = rampSpeed.aeroMapCatalog();
[ok,issues] = rampSpeed.validateAeroMapCatalog(catalog);
verifyTrue(testCase,ok);
verifyEmpty(testCase,issues);
badPath = [tempname '.csv'];
cleanup = onCleanup(@()deleteIfPresent(badPath));
writetable(table([0;1;2],[0;1;2]),badPath);
badCatalog = catalog;
badCatalog(1).path = string(badPath);
[ok,issues] = rampSpeed.validateAeroMapCatalog(badCatalog);
verifyFalse(testCase,ok);
verifyNotEmpty(testCase,issues);
verifyError(testCase,@()rampSpeed.validateAeroMapFile(badPath), ...
    'rampSpeed:invalidAeroMap');
end

function testBaselineDoesNotCallMutableCarConfig(testCase)
source = fileread(which('carConfigBaseline'));
verifyFalse(testCase,contains(source,'carConfig('));
end

function out = setField(in,name,value)
out = in;
out.(name) = value;
end

function deleteIfPresent(path)
if isfile(path)
    delete(path);
end
end
