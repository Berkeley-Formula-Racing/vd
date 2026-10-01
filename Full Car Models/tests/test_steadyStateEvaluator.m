function tests = test_steadyStateEvaluator
tests = functiontests(localfunctions);
end

function setupOnce(~)
testDir = fileparts(mfilename('fullpath'));
addpath(fileparts(testDir));
setup_paths
end

function testMemoizesRepeatedStateWithoutChangingEquations(testCase)
[cars,~] = carConfig();
car = cars{1,1};
P = [4,0.15,20,0.2,0.45,0,0,0.02,0.02];

evaluator = steadyStateEvaluator(car);
first = evaluator.evaluate(P);
second = evaluator.evaluate(P);
[engineRpm,beta,latAccel,longAccel,yawAccel,wheelAccel,~,~,Fzvirtual] = car.equations(P);

verifyEqual(testCase,evaluator.equationCalls(),1);
verifyEqual(testCase,first,second);
verifyEqual(testCase,first.engineRpm,engineRpm,'AbsTol',1e-12);
verifyEqual(testCase,first.beta,beta,'AbsTol',1e-12);
verifyEqual(testCase,first.latAccel,latAccel,'AbsTol',1e-12);
verifyEqual(testCase,first.longAccel,longAccel,'AbsTol',1e-12);
verifyEqual(testCase,first.yawAccel,yawAccel,'AbsTol',1e-12);
verifyEqual(testCase,first.wheelAccel,wheelAccel,'AbsTol',1e-12);
verifyEqual(testCase,first.Fzvirtual,Fzvirtual,'AbsTol',1e-12);
end

function testCachedConstraintMatchesCarConstraint(testCase)
[cars,~] = carConfig();
car = cars{1,1};
P = [4,0.15,20,0.2,0.45,0,0,0.02,0.02];
targetLatAccel = P(3)*P(5);

evaluator = steadyStateEvaluator(car);
[actualC,actualCeq] = steadyStateConstraint4(evaluator.evaluate(P),P,targetLatAccel);
[expectedC,expectedCeq] = car.constraint4(P,targetLatAccel);

verifyEqual(testCase,actualC,expectedC,'AbsTol',1e-12);
verifyEqual(testCase,actualCeq,expectedCeq,'AbsTol',1e-12);
end

function testEvaluateFullReusesTheCachedEquationEvaluation(testCase)
car = carConfigBaseline();
P = [4,0.15,20,0.2,0.45,0,0,0.02,0.02];

evaluator = steadyStateEvaluator(car);
state = evaluator.evaluate(P);
fullState = evaluator.evaluateFull(P);

verifyEqual(testCase,evaluator.equationCalls(),1);
verifyEqual(testCase,fullState.engineRpm,state.engineRpm,'AbsTol',1e-12);
verifyEqual(testCase,fullState.longAccel,state.longAccel,'AbsTol',1e-12);
verifyEqual(testCase,fullState.Fzvirtual,state.Fzvirtual,'AbsTol',1e-12);
verifyEqual(testCase,fullState.ssInfo,state.ssInfo);
verifyEqual(testCase,fullState.rideHeightContext,state.rideHeightContext);
end
