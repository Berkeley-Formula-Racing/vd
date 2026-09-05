function tests = test_lc0NDEvaluate
tests = functiontests(localfunctions);
end

function testUsesMeasuredLateralAndExplicitDonorScaling(testCase)
% Break caught: donor force level leaks directly into the target model, or
% the fitted combined-slip cap alters a feasible pure-axis force.
addpath(fileparts(fileparts(mfilename('fullpath'))));
target = lateralFit([-1;1],[-2;2]);
donorLat = lateralFit([-1;1],[-1;1]);
donorLong = longFit([-1;1],[-1;1]);
cal = struct('rhoMu',1,'rhoStiff',1,'couplingExponent',2);
model = lc0NDBuildModel(target,donorLong,donorLat,cal);

[fxPure,fyPure] = lc0NDEvaluate(model,0,1,100);
verifyEqual(testCase,fxPure,200,'AbsTol',1e-12);
verifyEqual(testCase,fyPure,0,'AbsTol',1e-12);

[fx,fy,info] = lc0NDEvaluate(model,1,1,100);
verifyEqual(testCase,fx,200/sqrt(2),'AbsTol',1e-12);
verifyEqual(testCase,fy,200/sqrt(2),'AbsTol',1e-12);
verifyEqual(testCase,info.longitudinal_mu_scale,2,'AbsTol',1e-12);
end

function fit = lateralFit(alpha,mu)
fit = struct('curve',table(alpha,mu,true(numel(alpha),1), ...
    'VariableNames',{'slip_angle_deg','mu_y','is_qualified'}), ...
    'peak_abs_mu_y',max(abs(mu)));
end

function fit = longFit(kappa,mu)
fit = struct('curve',table(kappa,mu,true(numel(kappa),1), ...
    'VariableNames',{'slip_ratio','mu_x','is_qualified'}), ...
    'peak_drive_mu',max(mu),'peak_brake_mu',-min(mu));
end
