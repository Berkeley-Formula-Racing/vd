function tests = test_lc0NDCompareForceModels
tests = functiontests(localfunctions);
end

function testComparesMatchedOperatingPointWithoutMutatingModel(testCase)
% Break caught: the legacy and experimental forces are evaluated at
% different states, making a visual comparison meaningless.
addpath(fileparts(fileparts(mfilename('fullpath'))));
model = makeModel();
cases = table(1,1,100,0,'VariableNames', ...
    {'alpha_deg','slip_ratio','Fz_N','camber_deg'});
comparison = lc0NDCompareForceModels(model,cases,@legacy);

verifyEqual(testCase,comparison.experimental_Fx_N,200/sqrt(2),'AbsTol',1e-12);
verifyEqual(testCase,comparison.experimental_Fy_N,200/sqrt(2),'AbsTol',1e-12);
verifyEqual(testCase,comparison.legacy_Fx_N,100,'AbsTol',1e-12);
verifyEqual(testCase,comparison.legacy_Fy_N,100,'AbsTol',1e-12);
verifyTrue(testCase,comparison.experimental_supported);
end

function [fx,fy] = legacy(alpha,kappa,fz,~)
fx = kappa.*fz;
fy = alpha.*fz;
end

function model = makeModel()
lat = struct('curve',table([-1;1],[-2;2],true(2,1), ...
    'VariableNames',{'slip_angle_deg','mu_y','is_qualified'}));
donorLat = struct('curve',table([-1;1],[-1;1],true(2,1), ...
    'VariableNames',{'slip_angle_deg','mu_y','is_qualified'}));
long = struct('curve',table([-1;1],[-1;1],true(2,1), ...
    'VariableNames',{'slip_ratio','mu_x','is_qualified'}));
model = lc0NDBuildModel(lat,long,donorLat, ...
    struct('rhoMu',1,'rhoStiff',1,'couplingExponent',2));
end
