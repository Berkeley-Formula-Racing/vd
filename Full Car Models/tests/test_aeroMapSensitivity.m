function tests = test_aeroMapSensitivity
tests = functiontests(localfunctions);
end

function testRejectsCarsWithoutSolvedGGData(testCase)
modelRoot = fileparts(which('carConfig'));
cd(modelRoot); setup_paths
[baseCell,~] = carConfig();
[carCell,plan] = aeroMapStarCases(baseCell);

verifyError(testCase,@() aeroMapSensitivity(carCell,plan), ...
    'aeroMapSensitivity:unsolvedCars');
end
