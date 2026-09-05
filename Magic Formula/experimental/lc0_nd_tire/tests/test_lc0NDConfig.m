function tests = test_lc0NDConfig
tests = functiontests(localfunctions);
end

function testDeclaresProvisionalLC0DonorExplicitly(testCase)
% Break caught: donor selection is implicit or based on an unrelated tyre.
addpath(fileparts(fileparts(mfilename('fullpath'))));
cfg = lc0NDConfig();

verifyTrue(testCase,isfield(cfg,'donor'));
verifyTrue(testCase,isfile(cfg.donor.file));
verifyTrue(testCase,contains(lower(cfg.donor.description),'lc0'));
verifyTrue(testCase,contains(lower(cfg.donor.description),'provisional'));
end
