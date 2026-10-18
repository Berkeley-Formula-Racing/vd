function tests = test_vehicleConfigRegistry
%TEST_VEHICLECONFIGREGISTRY Verify shared parameter-name discovery.
tests = functiontests(localfunctions);
end

function testRegistryMatchesCurrentDesign(testCase)
[~,~,X] = carConfig();
names = vehicleConfig.parameterNames();
verifyEqual(testCase,names,X.Properties.VariableNames);
verifyEqual(testCase,numel(names),57);
verifyEqual(testCase,numel(unique(names)),57);
end

function testDoeUsesEveryRegistryName(testCase)
[~,~,~,B] = carConfig();
study = DOEStudyConfig();
resolved = doeResolveStudy(study,B);
verifyTrue(testCase,all(ismember(cellstr(resolved.parameters.name), ...
    vehicleConfig.parameterNames())));
end
