function tests = test_rideHeightAeroContext
tests = functiontests(localfunctions);
end

function testWarmStartMatchesColdResultAtSameOperatingPoint(testCase)
car = carConfigBaseline();
longVel = 20;
wheelTorque = [37 37 61 61];

[cold,coldContext] = car.solveRideHeightAero(longVel,wheelTorque);
[warm,warmContext] = car.solveRideHeightAero(longVel,wheelTorque,coldContext);

verifyTrue(testCase,cold.converged);
verifyTrue(testCase,warm.converged);
verifyTrue(testCase,warmContext.usedWarmStart);
verifyFalse(testCase,warmContext.usedColdFallback);
verifyEqual(testCase,warm.iterations,0);
verifyInfoEquivalent(testCase,warm,cold);
end

function testNonlocalSpeedContextFallsBackToColdStart(testCase)
car = carConfigBaseline();
wheelTorque = [37 37 61 61];
[~,context] = car.solveRideHeightAero(20,wheelTorque);

[nonlocal,nonlocalContext] = car.solveRideHeightAero(30,wheelTorque,context);
cold = car.solveRideHeightAero(30,wheelTorque);

verifyFalse(testCase,nonlocalContext.usedWarmStart);
verifyFalse(testCase,nonlocalContext.usedColdFallback);
verifyInfoEquivalent(testCase,nonlocal,cold);
end

function testContextFromDifferentCarSetupIsRejected(testCase)
car = carConfigBaseline();
wheelTorque = [37 37 61 61];
[~,context] = car.solveRideHeightAero(20,wheelTorque);

otherCar = car;
otherCar.rideHeightAero.wheel_rate_front_Npm = ...
    1.05*otherCar.rideHeightAero.wheel_rate_front_Npm;
[~,otherContext] = otherCar.solveRideHeightAero(20,wheelTorque,context);

verifyFalse(testCase,otherContext.usedWarmStart);
end

function testContextFromDifferentAeroMapIsRejected(testCase)
car = carConfigBaseline();
wheelTorque = [37 37 61 61];
[~,context] = car.solveRideHeightAero(20,wheelTorque);

otherCar = car;
otherCar.aero.map = otherCar.aero.map.withCorrections(1.05,1,0);
[~,otherContext] = otherCar.solveRideHeightAero(20,wheelTorque,context);

verifyFalse(testCase,otherContext.usedWarmStart);
end

function testExplicitInitialGuessAndContextedEquationPath(testCase)
car = carConfigBaseline();
P = [0 0 20 0 0 0 0 0 0];
staticGuess = [car.rideHeightAero.static_front_ride_height_in; ...
    car.rideHeightAero.static_rear_ride_height_in];
seed = car.newRideHeightAeroContext(staticGuess);

[legacyOutputs,legacyContext] = evaluateWithContext(car,P,[]);
[contextOutputs,context] = evaluateWithContext(car,P,legacyContext);

verifyFalse(testCase,legacyContext.usedWarmStart);
verifyTrue(testCase,context.usedWarmStart);
verifyEqual(testCase,contextOutputs{4},legacyOutputs{4},'AbsTol',1e-12);
verifyEqual(testCase,contextOutputs{9},legacyOutputs{9},'AbsTol',1e-12);
verifyEqual(testCase,contextOutputs{16}.downforce, ...
    legacyOutputs{16}.downforce,'AbsTol',1e-12);

[seeded,seedContext] = car.solveRideHeightAero(20,[37 37 61 61],seed);
cold = car.solveRideHeightAero(20,[37 37 61 61]);
verifyTrue(testCase,seedContext.usedWarmStart);
verifyInfoEquivalent(testCase,seeded,cold);
verifyError(testCase,@() car.newRideHeightAeroContext([NaN;0]), ...
    'Car:invalidRideHeightAeroInitialGuess');
end

function [outputs,context] = evaluateWithContext(car,P,context)
outputs = cell(1,16);
[outputs{:},context] = car.equations(P,context);
end

function verifyInfoEquivalent(testCase,actual,expected)
fields = {'frontRideHeightIn','rearRideHeightIn','Fz_front_axle', ...
    'Fz_rear_axle','downforce','drag','downforce_front','downforce_rear', ...
    'cla','cda','D_f','D_r'};
for k = 1:numel(fields)
    field = fields{k};
    verifyEqual(testCase,actual.(field),expected.(field),'AbsTol',1e-7, ...
        ['Ride-height context changed ' field '.']);
end
verifyLessThanOrEqual(testCase,actual.residualIn,1e-7);
end
