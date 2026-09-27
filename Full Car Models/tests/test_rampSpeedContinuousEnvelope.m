function tests = test_rampSpeedContinuousEnvelope
% Contract tests for the ramp-speed continuous-envelope powertrain path.
tests = functiontests(localfunctions);
end

function testEnvelopeUsesContinuousRatioWithinPowertrainBounds(testCase)
[car,~] = carConfigBaseline();
envelope = rampSpeed.buildContinuousEnvelope(car,20);

verifyEqual(testCase,string(envelope.powertrainModel),"continuousEnvelope");
verifyEqual(testCase,envelope.speed_mps,20,"AbsTol",0);
reductions = car.powertrain.gears .* car.powertrain.final_drive .* ...
    car.powertrain.primary_reduction;
verifyGreaterThanOrEqual(testCase,envelope.drivetrainReduction, ...
    min(reductions)-1e-12);
verifyLessThanOrEqual(testCase,envelope.drivetrainReduction, ...
    max(reductions)+1e-12);
verifyLessThanOrEqual(testCase,envelope.engineRpm, ...
    car.powertrain.redline+1e-9);
verifyTrue(testCase,isfinite(envelope.engineRpm));
verifyTrue(testCase,isfinite(envelope.engineTorque_Nm));
verifyTrue(testCase,isfinite(envelope.wheelForce_N));
verifyGreaterThan(testCase,envelope.wheelForce_N,0);
end

function testCarAcceptsContinuousRatioWithoutMutatingTheCar(testCase)
[car,~] = carConfigBaseline();
state = rampSpeed.makePureLongState(car,20,1,0.02);
envelope = rampSpeed.buildContinuousEnvelope(car,20);
original = car;

[engineRpm,~,~,longAccel,~,~,~,currentGear] = car.equations( ...
    state,[],struct("continuousRatio",envelope.drivetrainReduction));

verifyTrue(testCase,isfinite(engineRpm));
verifyTrue(testCase,isfinite(longAccel));
verifyTrue(testCase,isnan(currentGear));
verifyEqual(testCase,car.powertrain.gears,original.powertrain.gears);
verifyEqual(testCase,car.R_sf,original.R_sf,"AbsTol",0);
end

function testEnvelopeRejectsNonpositiveSpeed(testCase)
[car,~] = carConfigBaseline();
verifyError(testCase,@() rampSpeed.buildContinuousEnvelope(car,0), ...
    "rampSpeed:invalidSpeed");
verifyError(testCase,@() rampSpeed.buildContinuousEnvelope(car,-1), ...
    "rampSpeed:invalidSpeed");
end
