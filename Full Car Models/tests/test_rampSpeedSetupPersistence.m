function tests = test_rampSpeedSetupPersistence
tests = functiontests(localfunctions);
end

function testSetupDefinitionsRoundTripWithStudy(testCase)
[~,config] = carConfigBaseline();
[cars,cases] = rampSpeed.buildSetupCatalog(config,config.defaultSetup);
request = struct('rampType',"longitudinal",'carRole',"auto", ...
    'settings',struct('speeds',5),'parallelRequested',false, ...
    'numWorkers',0,'checkpointPath',"",'appVersion',"test");
request.runCaseFcn = @captureRun;
[study,~] = rampSpeed.runStudy(cars,cases,request,struct('onProgress',[]));
path = [tempname '.mat'];
cleanup = onCleanup(@()deleteIfPresent(path));
rampSpeed.saveStudy(path,study);
loaded = rampSpeed.loadStudy(path,"test");
verifyFalse(testCase,loaded.readOnly);
verifyEqual(testCase,string(loaded.baselineVersion), ...
    string(config.baselineVersion));
verifyEqual(testCase,string(loaded.setupSpecifications(1).aeroMapId),"b26");
verifyEqual(testCase,loaded.setupSpecifications(1).rearArbStiffness_NmPerRad, ...
    config.defaultSetup.rearArbStiffness_NmPerRad);
verifyEqual(testCase,loaded.setupSpecifications(1).driverWeight_kg, ...
    config.defaultSetup.driverWeight_kg);
verifyEqual(testCase,loaded.setupSpecifications(1).rearWeightDistribution_percent, ...
    config.defaultSetup.rearWeightDistribution_percent);
verifyEqual(testCase,loaded.runs(1).runMeta.setupSpec.driverWeight_kg, ...
    config.defaultSetup.driverWeight_kg);
verifyEqual(testCase,loaded.runs(1).runMeta.setupSpec.rearWeightDistribution_percent, ...
    config.defaultSetup.rearWeightDistribution_percent);
verifyEqual(testCase,string(loaded.runs(1).runMeta.setupSpec.id), ...
    string(config.defaultSetup.id));
end

function testStudyWithoutSetupDefinitionsLoadsReadOnly(testCase)
study = rampSpeed.makeStudy("legacy-result");
study = rmfield(study,{'setupSpecifications','baselineVersion','readOnly'});
path = [tempname '.mat'];
cleanup = onCleanup(@()deleteIfPresent(path));
rampSpeed.saveStudy(path,study);
loaded = rampSpeed.loadStudy(path,"test");
verifyTrue(testCase,loaded.readOnly);
verifyTrue(testCase,~isfield(loaded,'setupSpecifications') || ...
    isempty(loaded.setupSpecifications));
end

function deleteIfPresent(path)
if isfile(path)
    delete(path);
end
end

function run = captureRun(car,caseInfo,request,~)
run = rampSpeed.makeRun(request.rampType,"",request.settings,caseInfo);
run.status = "complete";
run.runMeta.testCarRsf = car.R_sf;
end
