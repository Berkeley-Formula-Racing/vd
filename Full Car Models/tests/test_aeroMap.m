function tests = test_aeroMap
tests = functiontests(localfunctions);
end

function testReturnsTheExactBaselineCoefficients(testCase)
modelRoot = fileparts(which('carConfig'));
map = AeroMap(fullfile(modelRoot,'aeromap_b26.csv'));

aero = map.evaluate(0,0);

verifyEqual(testCase,aero.cla,2.98,'AbsTol',1e-12);
verifyEqual(testCase,aero.cda,1.27,'AbsTol',1e-12);
verifyEqual(testCase,aero.D_f,0.5364258,'AbsTol',1e-12);
verifyEqual(testCase,aero.D_r,0.4635742,'AbsTol',1e-12);
end

function testReturnsFiniteNearestValueOutsideMeasuredEnvelope(testCase)
modelRoot = fileparts(which('carConfig'));
map = AeroMap(fullfile(modelRoot,'aeromap_b26.csv'));

aero = map.evaluate(-2,2);

verifyTrue(testCase,all(isfinite([aero.cla aero.cda aero.D_f aero.D_r])));
verifyGreaterThan(testCase,aero.cla,0);
verifyGreaterThan(testCase,aero.cda,0);
verifyGreaterThanOrEqual(testCase,aero.D_f,0);
verifyLessThanOrEqual(testCase,aero.D_f,1);
end

function testCorrectionsApplyToEveryInterpolatedMapPoint(testCase)
modelRoot = fileparts(which('carConfig'));
baseMap = AeroMap(fullfile(modelRoot,'aeromap_b26.csv'));
correctedMap = baseMap.withCorrections(1.10,0.90,0.04);

base = baseMap.evaluate(-0.30,0.20);
actual = correctedMap.evaluate(-0.30,0.20);

verifyEqual(testCase,actual.cla,1.10*base.cla,'AbsTol',1e-12);
verifyEqual(testCase,actual.cda,0.90*base.cda,'AbsTol',1e-12);
verifyEqual(testCase,actual.D_f,base.D_f+0.04,'AbsTol',1e-12);
verifyEqual(testCase,actual.D_r,1-actual.D_f,'AbsTol',1e-12);
end
