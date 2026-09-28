function tests = test_lc0NDValidateDonorChoice
tests = functiontests(localfunctions);
end

function testAcceptsConfiguredLargerSameCompoundDonor(testCase)
addpath(fileparts(fileparts(mfilename('fullpath'))));
cfg = lc0NDConfig();
provenance = lc0NDValidateDonorChoice(cfg);
verifyTrue(testCase,provenance.same_compound);
verifyTrue(testCase,provenance.donor_is_larger);
verifyEqual(testCase,provenance.donor_size_in(1),18);
end

