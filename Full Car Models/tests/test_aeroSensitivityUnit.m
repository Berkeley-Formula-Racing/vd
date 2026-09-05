function tests = test_aeroSensitivityUnit
tests = functiontests(localfunctions);
end

function testUsesAeromapParameterUnits(testCase)
S = struct();
S.paramInfo = struct('name',"FrontRideHeightIn",'unit',"in");

verifyEqual(testCase,aeroSensitivityUnit(S,'FrontRideHeightIn','time'),'s per in');
verifyEqual(testCase,aeroSensitivityUnit(S,'FrontRideHeightIn','g'),'g per in');
end
