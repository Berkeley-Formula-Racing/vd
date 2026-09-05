function tests = test_lc0NDLateralSummary
tests = functiontests(localfunctions);
end

function testReportsPeakFrictionAndSmallSlipStiffness(testCase)
% Break caught: summary metrics use TTC sign conventions or fitted points
% outside the requested load/pressure/camber bin.
addpath(fileparts(fileparts(mfilename('fullpath'))));
data = struct( ...
    'slipAngle_deg',[-2;-1;0;1;2], ...
    'Fy_N',[-100;-80;0;80;100], ...
    'Fz_N',100*ones(5,1), ...
    'pressure_psi',12*ones(5,1), ...
    'camber_deg',zeros(5,1));
opts = struct('loadCenters_N',100,'loadHalfWidth_N',1, ...
    'pressureCenters_psi',12,'pressureHalfWidth_psi',0.1, ...
    'camberCenters_deg',0,'camberHalfWidth_deg',0.1, ...
    'smallSlipWindow_deg',1);

summary = lc0NDLateralSummary(data,opts);

verifyEqual(testCase,summary.n_samples,5);
verifyEqual(testCase,summary.peak_abs_Fy_N,100,'AbsTol',1e-12);
verifyEqual(testCase,summary.peak_mu_y,1,'AbsTol',1e-12);
verifyEqual(testCase,summary.stiffness_N_per_deg,80,'AbsTol',1e-12);
verifyEqual(testCase,summary.alpha_at_peak_deg,-2,'AbsTol',1e-12);
end

function testFlagsIncompleteSlipSweep(testCase)
% Break caught: a partial transient is mistaken for a measured peak curve.
addpath(fileparts(fileparts(mfilename('fullpath'))));
data = struct('slipAngle_deg',[0;0.5], 'Fy_N',[0;40], ...
    'Fz_N',100*ones(2,1), 'pressure_psi',12*ones(2,1), ...
    'camber_deg',zeros(2,1));
opts = struct('loadCenters_N',100,'loadHalfWidth_N',1, ...
    'pressureCenters_psi',12,'pressureHalfWidth_psi',0.1, ...
    'camberCenters_deg',0,'camberHalfWidth_deg',0.1, ...
    'smallSlipWindow_deg',1,'minSamples',3,'minSlipAngleSpan_deg',4);

summary = lc0NDLateralSummary(data,opts);

verifyEqual(testCase,summary.slipAngle_span_deg,0.5,'AbsTol',1e-12);
verifyFalse(testCase,summary.is_complete_sweep);
end

function testMarksEmptyBinsIncomplete(testCase)
% Break caught: empty requested bins must remain reportable rather than
% turning a coverage flag into NaN and crashing the analysis runner.
addpath(fileparts(fileparts(mfilename('fullpath'))));
data = struct('slipAngle_deg',0, 'Fy_N',0, 'Fz_N',100, ...
    'pressure_psi',12, 'camber_deg',0);
opts = struct('loadCenters_N',200,'loadHalfWidth_N',1, ...
    'pressureCenters_psi',12,'pressureHalfWidth_psi',0.1, ...
    'camberCenters_deg',0,'camberHalfWidth_deg',0.1, ...
    'smallSlipWindow_deg',1);

summary = lc0NDLateralSummary(data,opts);

verifyEqual(testCase,summary.n_samples,0);
verifyFalse(testCase,summary.is_complete_sweep);
end

function testReportsPositiveCorneringStiffnessMagnitudeForTTCSign(testCase)
% Break caught: TTC's FY-versus-SA sign convention should not make the
% engineering stiffness diagnostic appear negative.
addpath(fileparts(fileparts(mfilename('fullpath'))));
data = struct('slipAngle_deg',[-1;0;1], 'Fy_N',[80;0;-80], ...
    'Fz_N',100*ones(3,1), 'pressure_psi',12*ones(3,1), ...
    'camber_deg',zeros(3,1));
opts = struct('loadCenters_N',100,'loadHalfWidth_N',1, ...
    'pressureCenters_psi',12,'pressureHalfWidth_psi',0.1, ...
    'camberCenters_deg',0,'camberHalfWidth_deg',0.1, ...
    'smallSlipWindow_deg',1);

summary = lc0NDLateralSummary(data,opts);

verifyEqual(testCase,summary.stiffness_N_per_deg,80,'AbsTol',1e-12);
end
