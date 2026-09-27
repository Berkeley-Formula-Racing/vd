function tests = test_nondimensionalTire
tests = functiontests(localfunctions);
end

function testZeroLoadIsZeroAndUnsupportedStatesAreReported(testCase)
% Break caught: the vehicle-facing adapter turns contact loss into a force
% singularity or silently extrapolates beyond TTC coverage.
addpath(fileparts(fileparts(mfilename('fullpath'))));
setup_paths;
addpath(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), ...
    'Magic Formula','experimental','lc0_nd_tire'));
model = syntheticModel();
tire = NondimensionalTire(model,12,1,"nominal");

[fx0,fy0,info0] = tire.evaluate(0,0,0,0);
verifyEqual(testCase,fx0,0,'AbsTol',1e-12);
verifyEqual(testCase,fy0,0,'AbsTol',1e-12);
verifyFalse(testCase,info0.is_extrapolated);

[fx,fy,info] = tire.evaluate(8,0.25,100,0);
verifyTrue(testCase,isfinite(fx));
verifyTrue(testCase,isfinite(fy));
verifyTrue(testCase,info.is_extrapolated);
verifyGreaterThan(testCase,info.clamp_count,0);
verifyEqual(testCase,info.scenario,"nominal");
end

function testNegativeLoadIsRejected(testCase)
model = syntheticModel();
tire = NondimensionalTire(model,12,1,"nominal");
verifyError(testCase,@() tire.F_x(0,0,-1,0), ...
    'NondimensionalTire:negativeLoad');
end

function testPureAxesRemainUnmodified(testCase)
model = syntheticModel();
tire = NondimensionalTire(model,12,1,"nominal");

fx = tire.F_x(0,0.1,100,0);
fy = tire.F_y(2,0,100,0);
verifyEqual(testCase,fx,10,'AbsTol',1e-12);
verifyEqual(testCase,fy,50,'AbsTol',1e-12);
end

function model = syntheticModel()
target = struct('curve',table([-2;0;2],[0.5;0;0.5],true(3,1), ...
    'VariableNames',{'slip_angle_deg','mu_y','is_qualified'}), ...
    'peak_abs_mu_y',0.5);
donor = struct('curve',table([-0.1;0;0.1],[-0.1;0;0.1],true(3,1), ...
    'VariableNames',{'slip_ratio','mu_x','is_qualified'}), ...
    'peak_drive_mu',0.1,'peak_brake_mu',0.1);
model = lc0NDBuildModel(target,donor,target, ...
    struct('rhoMu',1,'rhoStiff',1,'couplingExponent',1.3));
end
