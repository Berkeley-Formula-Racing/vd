function tests = test_lc0NDUncertaintyScenarios
tests = functiontests(localfunctions);
end
function testBuildsNamedLowNominalHighDonorTransferCases(testCase)
% Break caught: unsupported donor force levels are presented as one exact
% target-tire prediction instead of an explicit uncertainty envelope.
addpath(fileparts(fileparts(mfilename('fullpath'))));
calibration = struct('rhoMu',1,'rhoStiff',1,'couplingExponent',1.3);

scenarios = lc0NDUncertaintyScenarios(calibration);

verifyEqual(testCase,{scenarios.name},{'low','nominal','high'});
verifyLessThan(testCase,scenarios(1).rhoMu,scenarios(2).rhoMu);
verifyLessThan(testCase,scenarios(2).rhoMu,scenarios(3).rhoMu);
verifyLessThan(testCase,scenarios(1).rhoStiff,scenarios(2).rhoStiff);
verifyLessThan(testCase,scenarios(2).rhoStiff,scenarios(3).rhoStiff);
verifyGreaterThanOrEqual(testCase,[scenarios.couplingExponent],1);
end
