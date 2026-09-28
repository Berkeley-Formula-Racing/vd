function tests = test_rampSpeedDriverWeight
tests = functiontests(localfunctions);
end

function testBaselineExposesDriverAndRearDistribution(testCase)
[car,config] = carConfigBaseline();
verifyEqual(testCase,config.defaults.driverWeight_kg,59);
verifyEqual(testCase,config.defaults.rearWeightDistribution_percent,51.2);
verifyEqual(testCase,config.defaultSetup.driverWeight_kg,59);
verifyEqual(testCase,config.defaultSetup.rearWeightDistribution_percent,51.2);
verifyEqual(testCase,car.M,162 + 59,'AbsTol',1e-12);
verifyEqual(testCase,car.l_f/car.W_b,0.512,'AbsTol',1e-12);
end

function testDriverWeightAndDistributionReachGeneratedCar(testCase)
[~,config] = carConfigBaseline();
setup = config.defaultSetup;
setup.driverWeight_kg = 80;
setup.rearWeightDistribution_percent = 54;
[car,normalized,derived] = rampSpeed.buildCarFromSetup(config,setup);
verifyEqual(testCase,car.M,162 + 80,'AbsTol',1e-12);
verifyEqual(testCase,car.l_f/car.W_b,0.54,'AbsTol',1e-12);
verifyEqual(testCase,normalized.driverWeight_kg,80);
verifyEqual(testCase,normalized.rearWeightDistribution_percent,54);
verifyEqual(testCase,derived.totalMass_kg,242,'AbsTol',1e-12);
verifyEqual(testCase,derived.frontWeightDistribution_percent,46,'AbsTol',1e-12);
verifyEqual(testCase,derived.staticRearLoad_N,242*9.81*0.54,'AbsTol',1e-10);
verifyEqual(testCase,derived.staticFrontLoad_N,242*9.81*0.46,'AbsTol',1e-10);
end

function testDriverWeightAndDistributionValidation(testCase)
[~,config] = carConfigBaseline();
bad = config.defaultSetup;
bad.driverWeight_kg = -1;
verifyError(testCase,@()rampSpeed.buildCarFromSetup(config,bad), ...
    'rampSpeed:invalidSetup');
bad = config.defaultSetup;
bad.rearWeightDistribution_percent = 100;
verifyError(testCase,@()rampSpeed.buildCarFromSetup(config,bad), ...
    'rampSpeed:invalidSetup');
end

function testDuplicateRetainsIndependentWeightOverrides(testCase)
[~,config] = carConfigBaseline();
copy = rampSpeed.duplicateSetup(config.defaultSetup,"driver-heavy", ...
    "Driver heavy");
copy.driverWeight_kg = 90;
copy.rearWeightDistribution_percent = 53;
[car,normalized] = rampSpeed.buildCarFromSetup(config,copy);
verifyEqual(testCase,config.defaultSetup.driverWeight_kg,59);
verifyEqual(testCase,config.defaultSetup.rearWeightDistribution_percent,51.2);
verifyEqual(testCase,normalized.id,"driver-heavy");
verifyEqual(testCase,car.M,252,'AbsTol',1e-12);
verifyEqual(testCase,car.l_f/car.W_b,0.53,'AbsTol',1e-12);
end