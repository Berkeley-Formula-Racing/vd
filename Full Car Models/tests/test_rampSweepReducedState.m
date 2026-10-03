function tests = test_rampSweepReducedState
tests = functiontests(localfunctions);
end

function testCoastStateReconstructionKeepsFixedControls(testCase)
z = [4.25,-0.17];
P = rampSweepState(z,20,8.5,false);

verifySize(testCase,P,[1 9]);
verifyEqual(testCase,P,[4.25,0,20,-0.17,8.5/20,0,0,0,0], ...
    'AbsTol',eps);
end

function testBalancedStateReconstructionRetainsIndependentControls(testCase)
z = [5.1,0.23,0.4,0.08];
P = rampSweepState(z,25,9.75,true);

verifyEqual(testCase,P,[5.1,0.4,25,0.23,9.75/25,0,0,0.08,0.08], ...
    'AbsTol',eps);
end

function testReducedStateRejectsWrongDimension(testCase)
verifyError(testCase,@()rampSweepState([1,2,3],20,8.5,false), ...
    'rampSweepState:badState');
end

function testEquationEvaluatorKeepsExactStateAndModeSpecificContext(testCase)
[cars,~] = carConfig();
car = cars{1,1};
P = [2,0.2,10,0.05,0.5,0,0,0.01,0.01];

minimal = rampEquationEvaluator(car,false);
minimalState = minimal.evaluate(P);
[engineRpm,beta,latAccel,longAccel,yawAccel,wheelAccel,~,~,Fzvirtual] = ...
    car.equations(P);
verifyEqual(testCase,minimalState.engineRpm,engineRpm,'AbsTol',0);
verifyEqual(testCase,minimalState.beta,beta,'AbsTol',0);
verifyEqual(testCase,minimalState.latAccel,latAccel,'AbsTol',0);
verifyEqual(testCase,minimalState.longAccel,longAccel,'AbsTol',0);
verifyEqual(testCase,minimalState.yawAccel,yawAccel,'AbsTol',0);
verifyEqual(testCase,minimalState.wheelAccel,wheelAccel,'AbsTol',0);
verifyEqual(testCase,minimalState.Fzvirtual,Fzvirtual,'AbsTol',0);
verifyEqual(testCase,minimal.equationCalls(),1);
verifyEmpty(testCase,minimalState.rideHeightContext);

contextEvaluator = rampEquationEvaluator(car,true);
contextState = contextEvaluator.evaluate(P);
verifyEqual(testCase,contextState.latAccel,latAccel,'AbsTol',0);
verifyNotEmpty(testCase,contextState.rideHeightContext);
end

function testBalancedLowSpeedRampRetainsFeasiblePoints(testCase)
[cars,~] = carConfig();
settings = struct('speeds',5,'nRamp',3,'nBisect',0,'mode','balanced', ...
    'verbose',false);
R = rampSweep(cars{1,1},settings);

verifyGreaterThan(testCase,height(R.points),0);
verifyEqual(testCase,R.status,"completed");
end
