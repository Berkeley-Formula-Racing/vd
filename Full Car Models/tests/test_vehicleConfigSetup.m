function tests = test_vehicleConfigSetup
%TEST_VEHICLECONFIGSETUP Verify editable setup defaults and legacy sampling parity.
tests = functiontests(localfunctions);
end

function testSetupEditsAreIndependent(testCase)
[base,eventParams] = vehicleConfig.defaultSetup();
variant = base;
variant.parameters.wheelbase = 1.60;
variant.parameters.final_drive = 3.2;
verifyEqual(testCase,base.parameters.wheelbase,62*0.0254);
verifyEqual(testCase,base.parameters.final_drive,33/11);
verifyEqual(testCase,eventParams.accel_length,75);
verifyEqual(testCase,variant.parameters.wheelbase,1.60);
verifyEqual(testCase,variant.parameters.final_drive,3.2);
end

function testDefaultSetupMatchesLegacyRouteAndCapturedReference(testCase)
[base,eventParams] = vehicleConfig.defaultSetup();
verifyEqual(testCase,fieldnames(base), ...
    {'schemaVersion';'profile';'id';'label';'parameters';'fixed';'aeroMap';'suspension'});
verifyEqual(testCase,base.schemaVersion,1);
verifyEqual(testCase,base.profile,"qss");
verifyEqual(testCase,fieldnames(base.parameters), ...
    vehicleConfig.parameterNames()');

% carConfig remains the independent legacy construction route in this task.
[oldCars,oldEvents,~,oldBaseline] = carConfig();
actualParameters = struct2table(base.parameters,'AsArray',true);
actualParameters = actualParameters(:,vehicleConfig.parameterNames());
verifyEqual(testCase,actualParameters,oldBaseline);
verifyEqual(testCase,eventParams,oldEvents);
verifyEqual(testCase,base.fixed.gears,oldCars{1,1}.powertrain.gears);
verifyEqual(testCase,base.fixed.primary_reduction, ...
    oldCars{1,1}.powertrain.primary_reduction);
verifyEqual(testCase,base.fixed.torque_fn,oldCars{1,1}.powertrain.torque_fn);
verifyEqual(testCase,base.fixed.model_type,"legacy");
verifyEqual(testCase,base.fixed.model_uncertainty,"nominal");
verifyEqual(testCase,base.fixed.model_artifact, ...
    fullfile(fileparts(fileparts(which('carConfig'))), ...
        'Magic Formula','experimental','lc0_nd_tire','results', ...
        'lc0_nd_pacejka_comparison.mat'));
verifyEqual(testCase,base.fixed.Fx_parameters, ...
    oldCars{1,1}.tire.Fx_parameters);
verifyEqual(testCase,base.fixed.Fy_parameters, ...
    oldCars{1,1}.tire.Fy_parameters);
verifyEqual(testCase,base.fixed.friction_scaling_factor, ...
    oldCars{1,1}.tire.friction_scaling_factor);

referencePath = fullfile(fileparts(fileparts(which('carConfig'))),'artifacts', ...
    'car_configuration_refactor','2026-10-04-reference', ...
    'construction-reference.mat');
requireReference = strcmp(getenv('VD_TASK2_REQUIRE_REFERENCE'),'1');
if requireReference && ~isfile(referencePath)
    error('test_vehicleConfigSetup:missingReference', ...
        'VD_TASK2_REQUIRE_REFERENCE is set but the Task 1 snapshot is missing.');
end
if isfile(referencePath)
    loaded = load(referencePath,'reference');
    old = loaded.reference.qss.default;
    current = carConfigurationSnapshot(struct( ...
        'parameters',actualParameters, ...
        'eventParams',eventParams, ...
        'fixed',struct( ...
            'gears',base.fixed.gears, ...
            'primary_reduction',base.fixed.primary_reduction, ...
            'torque_fn',base.fixed.torque_fn, ...
            'Fx_parameters',base.fixed.Fx_parameters, ...
            'Fy_parameters',base.fixed.Fy_parameters, ...
            'friction_scaling_factor',base.fixed.friction_scaling_factor)));
    expected = struct( ...
        'parameters',old.baselineTable, ...
        'eventParams',old.eventParams, ...
        'fixed',struct( ...
            'gears',old.lapCar.values.powertrain.values.gears, ...
            'primary_reduction',old.lapCar.values.powertrain.values.primary_reduction, ...
            'torque_fn',old.lapCar.values.powertrain.values.torque_fn, ...
            'Fx_parameters',old.lapCar.values.tire.values.Fx_parameters, ...
            'Fy_parameters',old.lapCar.values.tire.values.Fy_parameters, ...
            'friction_scaling_factor', ...
                old.lapCar.values.tire.values.friction_scaling_factor));
    verifyEqual(testCase,current,expected);
end
end

function testDefaultTyreAssetsResolveFromModelDirectory(testCase)
modelDir = fileparts(which('carConfig'));
oldFolder = pwd;
restoreFolder = onCleanup(@() cd(oldFolder)); %#ok<NASGU>
cd(tempdir);
[base,~] = vehicleConfig.defaultSetup();
verifyEqual(testCase,base.fixed.Fx_parameters, ...
    cell2mat(load(fullfile(modelDir,'carComponents', ...
        'Fx_combined_parameters_run38_30.mat'),'Xbestcell').Xbestcell));
verifyEqual(testCase,base.fixed.Fy_parameters, ...
    cell2mat(load(fullfile(modelDir,'carComponents', ...
        'Lapsim_Fy_combined_parameters_1965run15.mat'),'Xbestcell').Xbestcell));
verifyEqual(testCase,base.aeroMap.path,fullfile(modelDir,'aeromap_b26.csv'));
verifyTrue(testCase,base.aeroMap.enabled);
verifyEqual(testCase,base.aeroMap.claScale,1);
verifyEqual(testCase,base.aeroMap.cdaScale,1);
verifyEqual(testCase,base.aeroMap.copOffset,0);
verifyEqual(testCase,base.suspension.rollSplitPolicy,"explicit");
end

function testFromLegacyGroupsMapsSampledAndFixedValues(testCase)
base = vehicleConfig.defaultSetup();
[cP,aP,eP,DTp,Bp,tP] = legacyGroupsFromSetup(base);
cP.mass = 177;
aP.cda = 1.23;
eP.redline = 9200;
DTp.final_drive = 3.6;
Bp.brake_distribution = 0.68;
tP.p_i = 11;
eP.gears = [2.1 1.7 1.3 1.0];
eP.primary_reduction = 2.25;
eP.torque_fn = @(rpm) rpm;
tP.Fx_parameters = [1 2 3];
tP.Fy_parameters = [4 5 6];
tP.friction_scaling_factor = 0.91;
tP.model_type = "lc0_nd";
tP.model_artifact = fullfile(fileparts(which('carConfig')), ...
    'Magic Formula','experimental','lc0_nd_tire','results', ...
    'lc0_nd_pacejka_comparison.mat');
tP.model_uncertainty = "high";
aP.use_aeromap = false;
aP.map_path = fullfile(fileparts(which('carConfig')),'aeromap_b26.csv');
aP.claScale = 1.07;
aP.cdaScale = 0.93;
aP.copOffset = 0.02;

setup = vehicleConfig.fromLegacyGroups(cP,aP,eP,DTp,Bp,tP);
verifyEqual(testCase,setup.parameters.mass,177);
verifyEqual(testCase,setup.parameters.cda,1.23);
verifyEqual(testCase,setup.parameters.redline,9200);
verifyEqual(testCase,setup.parameters.final_drive,3.6);
verifyEqual(testCase,setup.parameters.brake_distribution,0.68);
verifyEqual(testCase,setup.parameters.p_i,11);
verifyEqual(testCase,setup.fixed.gears,eP.gears);
verifyEqual(testCase,setup.fixed.primary_reduction,2.25);
verifyEqual(testCase,func2str(setup.fixed.torque_fn),func2str(eP.torque_fn));
verifyEqual(testCase,setup.fixed.Fx_parameters,tP.Fx_parameters);
verifyEqual(testCase,setup.fixed.Fy_parameters,tP.Fy_parameters);
verifyEqual(testCase,setup.fixed.friction_scaling_factor,0.91);
verifyEqual(testCase,setup.fixed.model_type,"lc0_nd");
verifyEqual(testCase,setup.fixed.model_artifact,tP.model_artifact);
verifyEqual(testCase,setup.fixed.model_uncertainty,"high");
verifyFalse(testCase,setup.aeroMap.enabled);
verifyEqual(testCase,setup.aeroMap.path,aP.map_path);
verifyEqual(testCase,setup.aeroMap.claScale,1.07);
verifyEqual(testCase,setup.aeroMap.cdaScale,0.93);
verifyEqual(testCase,setup.aeroMap.copOffset,0.02);
verifyEqual(testCase,fieldnames(setup.parameters), ...
    vehicleConfig.parameterNames()');
end

function testLegacyAeroMapFlagsKeepPrecedenceAndDoNotInventMissingPath(testCase)
base = vehicleConfig.defaultSetup();
[cP,aP,eP,DTp,Bp,tP] = legacyGroupsFromSetup(base);
missingMap = rmfield(aP,{'use_aeromap','map_path'});
setup = vehicleConfig.fromLegacyGroups(cP,missingMap,eP,DTp,Bp,tP);
verifyFalse(testCase,setup.aeroMap.enabled);
verifyEmpty(testCase,setup.aeroMap.path);
missingMap.use_aeromap = true;
setup = vehicleConfig.fromLegacyGroups(cP,missingMap,eP,DTp,Bp,tP);
verifyTrue(testCase,setup.aeroMap.enabled);
verifyEmpty(testCase,setup.aeroMap.path);

legacyMap = rmfield(aP,'use_aeromap');
legacyMap.map_enabled = false;
setup = vehicleConfig.fromLegacyGroups(cP,legacyMap,eP,DTp,Bp,tP);
verifyFalse(testCase,setup.aeroMap.enabled);
legacyMap.map_enabled = true;
setup = vehicleConfig.fromLegacyGroups(cP,legacyMap,eP,DTp,Bp,tP);
verifyTrue(testCase,setup.aeroMap.enabled);

conflicting = aP;
conflicting.map_enabled = true;
conflicting.use_aeromap = false;
setup = vehicleConfig.fromLegacyGroups(cP,conflicting,eP,DTp,Bp,tP);
verifyFalse(testCase,setup.aeroMap.enabled);
conflicting.map_enabled = false;
conflicting.use_aeromap = true;
setup = vehicleConfig.fromLegacyGroups(cP,conflicting,eP,DTp,Bp,tP);
verifyTrue(testCase,setup.aeroMap.enabled);
end

function testExplicitRowsKeepTheRequestedValues(testCase)
base = vehicleConfig.defaultSetup();
Q = table([155;170],[0.30;0.50], ...
    'VariableNames',{'mass','R_sf'});
[setups,X] = vehicleConfig.sampleSetups(base,"Explicit",Q);
verifySize(testCase,setups,[2 1]);
verifyEqual(testCase,X.mass,[155;170]);
verifyEqual(testCase,arrayfun(@(s)s.parameters.mass,setups),[155;170]);
verifyEqual(testCase,arrayfun(@(s)s.parameters.R_sf,setups),[0.30;0.50]);
verifyEqual(testCase,X.Properties.VariableNames, ...
    vehicleConfig.parameterNames());
verifyEqual(testCase,setups(1).id,"case-000001");
verifyEqual(testCase,setups(2).id,"case-000002");
end

function testEmptyExplicitInputReturnsZeroByFiftySevenTable(testCase)
base = vehicleConfig.defaultSetup();
Q = table(zeros(0,1),zeros(0,1), ...
    'VariableNames',{'mass','R_sf'});
[setups,X] = vehicleConfig.sampleSetups(base,"Explicit",Q);
verifySize(testCase,setups,[0 1]);
verifyTrue(testCase,isstruct(setups));
verifySize(testCase,X,[0 57]);
verifyEqual(testCase,X.Properties.VariableNames, ...
    vehicleConfig.parameterNames());
end

function testOneExplicitRowRetainsUnchangedBaselineFields(testCase)
base = vehicleConfig.defaultSetup();
Q = table(155,'VariableNames',{'mass'});
[setups,X] = vehicleConfig.sampleSetups(base,'Explicit',Q);
verifySize(testCase,setups,[1 1]);
verifyEqual(testCase,X.mass,155);
verifyEqual(testCase,X.driver_weight,64);
verifyEqual(testCase,setups.parameters.mass,155);
verifyEqual(testCase,setups.parameters.driver_weight,64);
end

function testUnknownExplicitVariableKeepsLegacyIdentifier(testCase)
base = vehicleConfig.defaultSetup();
Q = table(155,'VariableNames',{'unknown_value'});
verifyError(testCase,@() vehicleConfig.sampleSetups(base,"Explicit",Q), ...
    'parameters_loop:unknownExplicitParameter');
end

function testExplicitRejectsNonScalarBaselineWithLegacyIdentifier(testCase)
base = vehicleConfig.defaultSetup();
base.parameters.mass = [155 170];
Q = table(165,'VariableNames',{'mass'});
verifyError(testCase,@() vehicleConfig.sampleSetups(base,"Explicit",Q), ...
    'parameters_loop:nonScalarBaseline');
end

function testGridPreservesNdgridOrderAndStringCharNames(testCase)
base = vehicleConfig.defaultSetup();
base.parameters.mass = [155 170];
base.parameters.R_sf = [0.30 0.50];
[setups,X] = vehicleConfig.sampleSetups(base,'Grid',[]);
expectedMass = [155;170;155;170];
expectedSplit = [0.30;0.30;0.50;0.50];
verifyEqual(testCase,X.mass,expectedMass);
verifyEqual(testCase,X.R_sf,expectedSplit);
verifyEqual(testCase,arrayfun(@(s)s.parameters.mass,setups),expectedMass);
[~,fullFactorial] = vehicleConfig.sampleSetups(base,"FullFactorial",[]);
verifyEqual(testCase,X,fullFactorial);
end

function testDirectLhsAndRandomMatchLegacySamplerAndAmbientRng(testCase)
base = vehicleConfig.defaultSetup();
base.parameters.mass = [155 170];
base.parameters.R_sf = [0.30 0.50];
[cP,aP,eP,DTp,Bp,tP] = legacyGroupsFromSetup(base);
oldRng = rng;
restoreRng = onCleanup(@() rng(oldRng)); %#ok<NASGU>
for samplingType = ["LHS","Random"]
    rng(73,'twister');
    [~,expected] = parameters_loop(cP,aP,eP,DTp,Bp,tP,samplingType,4);
    expectedState = rng;
    rng(73,'twister');
    [~,actual] = vehicleConfig.sampleSetups(base,samplingType,4);
    actualState = rng;
    verifyEqual(testCase,actual,expected);
    verifyEqual(testCase,actualState,expectedState);
end
end

function [cP,aP,eP,DTp,Bp,tP] = legacyGroupsFromSetup(setup)
p = setup.parameters;
cNames = {'mass','driver_weight','accel_driver_weight','wheelbase', ...
    'weight_dist','track_width','wheel_radius','cg_height', ...
    'roll_center_height_front','roll_center_height_rear','R_sf','I_zz', ...
    'ackermann','camber_compliance_f','camber_compliance_r','static_r_toe', ...
    'roll_gradient_deg_per_g','rear_roll_camber_outer_deg_per_deg', ...
    'rear_roll_camber_inner_deg_per_deg','ride_camber_front_deg_per_in', ...
    'ride_camber_rear_deg_per_in','I_wheel','I_driveline','Crr', ...
    'qss_steering_limit_deg','spring_rate_front_lb_in', ...
    'spring_rate_rear_lb_in','motion_ratio_front','motion_ratio_rear', ...
    'static_front_ride_height_in','static_rear_ride_height_in'};
aNames = {'cda','cla','distribution','cla_p_deg_p','D_p_deg_p', ...
    'accel_cda','accel_cla','acc_cla_p_deg_p','acc_D_p_deg_p', ...
    'map_reference_front_ride_height_in','map_reference_rear_ride_height_in'};
eNames = {'redline','shift_point','shift_time'};
dtNames = {'final_drive','drivetrain_efficiency','G_d1','G_d2_overrun', ...
    'G_d2_driving'};
bNames = {'brake_distribution','max_braking_torque'};
tNames = {'gamma_f','gamma_r','p_i','grip_scaling_front','grip_scaling_rear'};
cP = projectFields(p,cNames);
aP = projectFields(p,aNames);
eP = projectFields(p,eNames);
DTp = projectFields(p,dtNames);
Bp = projectFields(p,bNames);
tP = projectFields(p,tNames);
eP.gears = setup.fixed.gears;
eP.primary_reduction = setup.fixed.primary_reduction;
eP.torque_fn = setup.fixed.torque_fn;
tP.Fx_parameters = setup.fixed.Fx_parameters;
tP.Fy_parameters = setup.fixed.Fy_parameters;
tP.friction_scaling_factor = setup.fixed.friction_scaling_factor;
tP.model_type = setup.fixed.model_type;
tP.model_artifact = setup.fixed.model_artifact;
tP.model_uncertainty = setup.fixed.model_uncertainty;
aP.use_aeromap = setup.aeroMap.enabled;
aP.map_path = setup.aeroMap.path;
aP.claScale = setup.aeroMap.claScale;
aP.cdaScale = setup.aeroMap.cdaScale;
aP.copOffset = setup.aeroMap.copOffset;
end

function output = projectFields(source,names)
output = struct();
for k = 1:numel(names)
    output.(names{k}) = source.(names{k});
end
end

function testFromLegacyGroupsDefaultsExplicitlyEmptyOptionalTireFields(testCase)
base = vehicleConfig.defaultSetup();
[cP,aP,eP,DTp,Bp,tP] = legacyGroupsFromSetup(base);
tP.model_type = [];
tP.model_artifact = [];
tP.model_uncertainty = [];

setup = vehicleConfig.fromLegacyGroups(cP,aP,eP,DTp,Bp,tP);
verifyEqual(testCase,setup.fixed.model_type,"legacy");
verifyEqual(testCase,setup.fixed.model_artifact,"");
verifyEqual(testCase,setup.fixed.model_uncertainty,"nominal");
end
