function tests = test_cuoSteerFromYaw
tests = functiontests(localfunctions);
end

function testReportsSignedSteerResidualAgainstYawKinematics(testCase)
% A positive steer margin is understeer; the sign must also work in a left turn.
[cuoRight,fromYawRight] = cuoSteerFromYaw(7.7106,0.5,10,2);
[cuoLeft,fromYawLeft]   = cuoSteerFromYaw(-5.7106,-0.5,10,2);

verifyEqual(testCase,fromYawRight,5.7106,'AbsTol',1e-4);
verifyEqual(testCase,cuoRight,2.0000,'AbsTol',1e-4);
verifyEqual(testCase,fromYawLeft,-5.7106,'AbsTol',1e-4);
verifyEqual(testCase,cuoLeft,0.0000,'AbsTol',1e-4);
end
