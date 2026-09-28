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

function testAppliesMeasuredTargetLoadScaling(testCase)
addpath(fileparts(fileparts(mfilename('fullpath'))));
target = lateralFit([-2;0;2],[0.5;0;0.5]);
donorLat = lateralFit([-1;1],[-1;1]);
donorLong = longFit([-1;1],[-1;1]);
cal = struct('rhoMu',1,'rhoStiff',1,'couplingExponent',1.3);
model = lc0NDBuildModel(target,donorLong,donorLat,cal);
model.target_scaling = struct( ...
    'reference',struct('pressure_psi',12,'camber_deg',0), ...
    'table',table([100;200],[12;12],[0;0],[1;0.8],[1;0.8], ...
        true(2,1),'VariableNames',{'load_center_N', ...
        'pressure_center_psi','camber_center_deg','mu_scale', ...
        'stiffness_scale','is_qualified'}));

[~,fy,info] = lc0NDEvaluate(model,2,0,200, ...
    struct('outOfRange',"clamp",'pressurePsi',12,'camberDeg',0));
verifyEqual(testCase,fy,64,'AbsTol',1e-12);
verifyEqual(testCase,info.lateral_mu_scale,0.8,'AbsTol',1e-12);
verifyEqual(testCase,info.lateral_stiffness_scale,0.8,'AbsTol',1e-12);
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
