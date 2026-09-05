function tests = test_lc0NDFitLongitudinal
tests = functiontests(localfunctions);
end

function testBuildsLoadNormalizedReferenceCurve(testCase)
% Break caught: a donor longitudinal reference is fit in raw force, which
% makes it unusable across target and donor load scales.
addpath(fileparts(fileparts(mfilename('fullpath'))));
donor = struct('slipRatio',[-.1;0;.1], 'slipAngle_deg',zeros(3,1), ...
    'Fx_N',[-100;0;100], 'Fz_N',100*ones(3,1), ...
    'pressure_psi',12*ones(3,1), 'camber_deg',zeros(3,1));
opts = struct('load_N',100,'loadHalfWidth_N',1, ...
    'pressure_psi',12,'pressureHalfWidth_psi',0.1, ...
    'camber_deg',0,'camberHalfWidth_deg',0.1, ...
    'maxAbsSlipAngle_deg',0.1,'kappaGrid',[-.1;0;.1], ...
    'kappaHalfWidth',.001,'minSamplesPerBin',1);

fit = lc0NDFitLongitudinal(donor,opts);

verifyEqual(testCase,fit.n_selected,3);
verifyEqual(testCase,fit.curve.mu_x,[-1;0;1],'AbsTol',1e-12);
verifyEqual(testCase,fit.curve.n_samples,[1;1;1]);
verifyEqual(testCase,fit.peak_drive_mu,1,'AbsTol',1e-12);
verifyEqual(testCase,fit.peak_brake_mu,1,'AbsTol',1e-12);
end
