function tests = test_lc0NDFitCouplingExponent
tests = functiontests(localfunctions);
end

function testRecoversEllipseExponentFromSyntheticCombinedPoint(testCase)
% Break caught: the coupling fit does not preserve the stated normalized
% pure-force curves or chooses an exponent unrelated to force error.
addpath(fileparts(fileparts(mfilename('fullpath'))));
longFit = struct('curve',table([-1;1],[-1;1],[true;true], ...
    'VariableNames',{'slip_ratio','mu_x','is_qualified'}));
latFit = struct('curve',table([-1;1],[-1;1],[true;true], ...
    'VariableNames',{'slip_angle_deg','mu_y','is_qualified'}));
combined = struct('slipRatio',1,'slipAngle_deg',1, ...
    'Fx_N',100/sqrt(2),'Fy_N',100/sqrt(2),'Fz_N',100, ...
    'pressure_psi',12,'camber_deg',0);
opts = struct('load_N',100,'loadHalfWidth_N',1, ...
    'pressure_psi',12,'pressureHalfWidth_psi',.1, ...
    'camber_deg',0,'camberHalfWidth_deg',.1, ...
    'minAbsSlipRatio',.01,'minAbsSlipAngle_deg',.1, ...
    'pGrid',[1;2;3],'minPoints',1);

fit = lc0NDFitCouplingExponent(combined,longFit,latFit,opts);

verifyEqual(testCase,fit.n_points,1);
verifyEqual(testCase,fit.best_exponent,2,'AbsTol',1e-12);
verifyEqual(testCase,fit.rmse_mu_at_best,0,'AbsTol',1e-12);
end
