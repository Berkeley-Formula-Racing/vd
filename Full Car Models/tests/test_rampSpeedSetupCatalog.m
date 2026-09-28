function tests = test_rampSpeedSetupCatalog
tests = functiontests(localfunctions);
end

function testCatalogUsesOneCarPerSetupAndRetainsProvenance(testCase)
[~,config] = carConfigBaseline();
copy = rampSpeed.duplicateSetup(config.defaultSetup,"rear-high","Rear ARB high");
copy.rearArbStiffness_NmPerRad = config.options.rearArbStiffness_NmPerRad(end);
[cars,cases,designTable] = rampSpeed.buildSetupCatalog(config, ...
    [config.defaultSetup,copy]);
verifySize(testCase,cars,[2 1]);
verifySize(testCase,cases,[2 1]);
verifyEqual(testCase,height(designTable),2);
verifyEqual(testCase,string(cases(1).setupSpec.baselineVersion), ...
    string(config.baselineVersion));
verifyTrue(testCase,cases(1).isBaseline);
verifyFalse(testCase,cases(2).isBaseline);
verifyEqual(testCase,string(cases(2).setupSpec.id),"rear-high");
verifyTrue(testCase,all(ismember( ...
    {'driverWeight_kg','rearWeightDistribution_percent'}, ...
    designTable.Properties.VariableNames)));
verifyEqual(testCase,designTable.R_sf(2),cars{2}.R_sf,'AbsTol',1e-12);
end

function testSameSetupCarIsUsedForBothRampTypes(testCase)
[~,config] = carConfigBaseline();
[cars,cases] = rampSpeed.buildSetupCatalog(config,config.defaultSetup);
request = baseRequest("lateral");
request.runCaseFcn = @captureRun;
[lateral,~] = rampSpeed.runStudy(cars,cases,request,struct());
request.rampType = "longitudinal";
request.carRole = "auto";
[longitudinal,~] = rampSpeed.runStudy(cars,cases,request,struct());
verifyEqual(testCase,lateral.runs(1).runMeta.testCarRsf, ...
    longitudinal.runs(1).runMeta.testCarRsf,'AbsTol',1e-12);
verifyEqual(testCase,string(lateral.runs(1).type),"lateral");
verifyEqual(testCase,string(longitudinal.runs(1).type),"longitudinal");
end

function request = baseRequest(type)
request = struct('rampType',type,'carRole',"auto", ...
    'settings',struct('speeds',5,'mode',"coast"), ...
    'parallelRequested',false,'numWorkers',0,'checkpointPath',"", ...
    'appVersion',"test");
end

function run = captureRun(car,caseInfo,request,~)
run = rampSpeed.makeRun(request.rampType,"coast",request.settings,caseInfo);
run.status = "complete";
run.runMeta.testCarRsf = car.R_sf;
end
