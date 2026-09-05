function tests = test_lc0NDFitLateral
tests = functiontests(localfunctions);
end

function testBuildsLoadNormalizedPureLateralReference(testCase)
% Break caught: a coupling fit uses raw donor lateral force rather than the
% pure lateral coefficient at its stated reference condition.
addpath(fileparts(fileparts(mfilename('fullpath'))));
donor = struct('slipRatio',zeros(3,1), 'slipAngle_deg',[-1;0;1], ...
    'Fy_N',[-100;0;100], 'Fz_N',100*ones(3,1), ...
    'pressure_psi',12*ones(3,1), 'camber_deg',zeros(3,1));
opts = struct('load_N',100,'loadHalfWidth_N',1, ...
    'pressure_psi',12,'pressureHalfWidth_psi',.1, ...
    'camber_deg',0,'camberHalfWidth_deg',.1, ...
    'maxAbsSlipRatio',.001,'alphaGrid',[-1;0;1], ...
    'alphaHalfWidth_deg',.01,'minSamplesPerBin',1);

fit = lc0NDFitLateral(donor,opts);

verifyEqual(testCase,fit.n_selected,3);
verifyEqual(testCase,fit.curve.mu_y,[-1;0;1],'AbsTol',1e-12);
verifyEqual(testCase,fit.peak_abs_mu_y,1,'AbsTol',1e-12);
end
