function tests = test_camberKinematicsConfig
tests = functiontests(localfunctions);
end

function testRideCamberAddsMirroredCamberUnderCompression(testCase)
cfg = config(-0.4,-0.2,0.68,0.58,-0.592);
base = Camber_Evaluation(10,0,0,0,-1,-1,0,0,cfg,[0;0]);
bump = Camber_Evaluation(10,0,0,0,-1,-1,0,0,cfg,[0.5;1.0]);

verifyEqual(testCase,bump-base,[0.2;-0.2;0.2;-0.2],'AbsTol',1e-12);
end

function testRideCamberUsesEachWheelBumpCompression(testCase)
cfg = config(-0.4,-0.2,0.68,0.58,-0.592);
base = Camber_Evaluation(10,0,0,0,-1,-1,0,0,cfg,zeros(4,1));
bump = Camber_Evaluation(10,0,0,0,-1,-1,0,0,cfg,[0.5;0.25;1.0;0.5]);

verifyEqual(testCase,bump-base,[0.2;-0.1;0.2;-0.1],'AbsTol',1e-12);
end

function testRearRollCamberUsesConfiguredGains(testCase)
off = config(0,0,0,0,0);
on = config(0,0,1.0,2.0,-3.0);
base = Camber_Evaluation(10,0.981,0,0,-1,-1,0,0,off,[0;0]);
roll = Camber_Evaluation(10,0.981,0,0,-1,-1,0,0,on,[0;0]);

% v*yaw/g = 1 g, so the configured roll angle is exactly one degree.
verifyEqual(testCase,roll(3)-base(3),2.0,'AbsTol',1e-12);
verifyEqual(testCase,roll(4)-base(4),-3.0,'AbsTol',1e-12);
end

function testCarConfigStoresCamberKinematics(testCase)
[cars,~] = carConfig();
cfg = cars{1,1}.camberKinematics;

verifyEqual(testCase,cfg.roll_gradient_deg_per_g,0.68);
verifyEqual(testCase,cfg.rear_roll_camber_outer_deg_per_deg,0.58);
verifyEqual(testCase,cfg.rear_roll_camber_inner_deg_per_deg,-0.592);
verifyEqual(testCase,cfg.ride_camber_front_deg_per_in,0);
verifyEqual(testCase,cfg.ride_camber_rear_deg_per_in,0);
end

function cfg = config(frontRide,rearRide,rollGradient,rearOuter,rearInner)
cfg = struct( ...
    'roll_gradient_deg_per_g',rollGradient, ...
    'rear_roll_camber_outer_deg_per_deg',rearOuter, ...
    'rear_roll_camber_inner_deg_per_deg',rearInner, ...
    'ride_camber_front_deg_per_in',frontRide, ...
    'ride_camber_rear_deg_per_in',rearRide);
end
