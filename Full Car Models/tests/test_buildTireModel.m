function tests = test_buildTireModel
tests = functiontests(localfunctions);
end

function testSelectsNondimensionalModelFromVersionedArtifact(testCase)
% Break caught: the new tire model exists in isolation but the vehicle
% factory silently continues constructing the legacy Tire2 object.
addpath(fileparts(fileparts(mfilename('fullpath'))));
setup_paths;
addpath(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), ...
    'Magic Formula','experimental','lc0_nd_tire'));
model = syntheticModel(); %#ok<NASGU>
root = tempname;
mkdir(root);
cleaner = onCleanup(@() rmdir(root,'s')); %#ok<NASGU>
artifact = fullfile(root,'model.mat');
save(artifact,'model');

tire = buildTireModel(struct('model_type',"lc0_nd", ...
    'model_artifact',artifact,'model_uncertainty',"nominal", ...
    'p_i',12,'friction_scaling_factor',1));

verifyClass(testCase,tire,'NondimensionalTire');
verifyEqual(testCase,tire.uncertainty_mode,"nominal");
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
