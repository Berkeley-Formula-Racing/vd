function tests = test_carConfigCompatibility
%TEST_CARCONFIGCOMPATIBILITY Verify named opt-in and legacy wrapper behavior.
tests = functiontests(localfunctions);
end

function testAllLegacyInputAndOutputAritiesKeepCellShape(testCase)
one = carConfig();
[two,eventTwo] = carConfig("Grid");
[three,eventThree,tableThree] = carConfig("FullFactorial",[]);
[four,eventFour,tableFour,baselineFour] = ...
    carConfig("Grid",[],"legacy");
[five,eventFive,tableFive,baselineFive] = ...
    carConfig("FullFactorial",[],"legacy","nominal");

verifySize(testCase,one,[1 2]);
verifySize(testCase,two,[1 2]);
verifySize(testCase,three,[1 2]);
verifySize(testCase,four,[1 2]);
verifySize(testCase,five,[1 2]);
verifyEqual(testCase,eventTwo,eventThree);
verifyEqual(testCase,eventThree,eventFour);
verifyEqual(testCase,eventFour,eventFive);
verifySize(testCase,tableThree,[1 57]);
verifySize(testCase,tableFour,[1 57]);
verifySize(testCase,baselineFour,[1 57]);
verifyEqual(testCase,tableFour,baselineFour);
verifyEqual(testCase,tableFive,baselineFive);
end

function testNamedOutputPacksExactlyLikeLegacyExplicitOutput(testCase)
Q = table([155;170],[0.30;0.50], ...
    'VariableNames',{'mass','R_sf'});
[legacy,eventParams,designTable,baselineTable] = carConfig("Explicit",Q);
[named,eventParams2,designTable2,baselineTable2] = ...
    carConfig("Explicit",Q,"legacy","nominal","cases");
verifySize(testCase,named,[2 1]);
verifyEqual(testCase,string({named.id}).',["case-000001";"case-000002"]);
verifyEqual(testCase,[named.designRow].',[1;2]);
verifyEqual(testCase,eventParams2,eventParams);
verifyEqual(testCase,designTable2,designTable);
verifyEqual(testCase,baselineTable2,baselineTable);
verifyEqual(testCase,legacy{1,1}.M,219);
verifyEqual(testCase,legacy{1,2}.M,214);
packed = vehicleConfig.toLegacyCells(named);
verifySize(testCase,packed,[2 2]);
for i = 1:2
    verifyEqual(testCase,carConfigurationSnapshot(packed{i,1}), ...
        carConfigurationSnapshot(legacy{i,1}));
    verifyEqual(testCase,carConfigurationSnapshot(packed{i,2}), ...
        carConfigurationSnapshot(legacy{i,2}));
    for column = 1:width(designTable)
        name = designTable.Properties.VariableNames{column};
        verifyEqual(testCase,named(i).setup.parameters.(name), ...
            designTable{i,column});
    end
end
end

function testEmptyAndOneRowExplicitTablesHaveStableShapes(testCase)
emptyQ = table(zeros(0,1),zeros(0,1), ...
    'VariableNames',{'mass','R_sf'});
[emptyCells,~,emptyTable,emptyBaseline] = carConfig("Explicit",emptyQ);
[emptyCases,~,emptyTable2,emptyBaseline2] = ...
    carConfig("Explicit",emptyQ,"legacy","nominal","cases");
verifySize(testCase,emptyCells,[0 2]);
verifySize(testCase,emptyCases,[0 1]);
verifySize(testCase,emptyTable,[0 57]);
verifyEqual(testCase,emptyTable2,emptyTable);
verifyEqual(testCase,emptyBaseline2,emptyBaseline);

oneRow = table(155,0.30,'VariableNames',{'mass','R_sf'});
[singleCase,~,singleTable] = ...
    carConfig("Explicit",oneRow,"legacy","nominal","cases");
verifySize(testCase,singleCase,[1 1]);
verifySize(testCase,singleTable,[1 57]);
verifyEqual(testCase,singleCase.id,"case-000001");
verifyEqual(testCase,singleCase.designRow,1);
verifyEqual(testCase,singleCase.setup.parameters.mass,155);
verifyEqual(testCase,singleCase.setup.parameters.R_sf,0.30);
verifySize(testCase,vehicleConfig.toLegacyCells(singleCase),[1 2]);
end

function testDirectParametersLoopPreservesMapValidationAndAmbientLhs(testCase)
base = vehicleConfig.defaultSetup();
[cP,aP,eP,DTp,Bp,tP] = legacyGroupsFromSetup(base);

aP.map_path = fullfile(repositoryRoot(),'Full Car Models','missing-map.csv');
aP.use_aeromap = true;
aP.map_enabled = false;
verifyError(testCase,@() parameters_loop(cP,aP,eP,DTp,Bp,tP, ...
    "FullFactorial",[]),'parameters_loop:aeroMapNotFound');

aP.use_aeromap = false;
aP.map_enabled = true;
[disabledCars,~] = parameters_loop(cP,aP,eP,DTp,Bp,tP,"FullFactorial",[]);
verifyFalse(testCase,disabledCars{1,1}.aero.hasMap());

aP.use_aeromap = true;
aP.map_enabled = false;
aP.map_path = base.aeroMap.path;
[enabledCars,~] = parameters_loop(cP,aP,eP,DTp,Bp,tP,"FullFactorial",[]);
verifyTrue(testCase,enabledCars{1,1}.aero.hasMap());

cP.mass = [160 164];
rng(21,'twister');
U = rand(3,57);
for column = 1:57
    U(:,column) = (randperm(3)' - U(:,column))/3;
end
expectedMass = 160 + U(:,1)*4;
expectedNext = rand;
rng(21,'twister');
[lhsCars,lhsTable] = parameters_loop(cP,aP,eP,DTp,Bp,tP,"LHS",3);
actualNext = rand;
verifySize(testCase,lhsCars,[3 2]);
verifyEqual(testCase,lhsTable.mass,expectedMass,'AbsTol',1e-12);
verifyEqual(testCase,actualNext,expectedNext,'AbsTol',1e-12);
end

function testDeprecatedCarConfigDoeKeepsWarningAndRngState(testCase)
for request = ["LHS","Random"]
    rng(19,'twister');
    stateBefore = rng;
    [cars,~,designTable] = carConfig(request,3);
    stateAfter = rng;
    [~,warningId] = lastwarn;
    verifySize(testCase,cars,[3 2]);
    verifySize(testCase,designTable,[3 57]);
    verifyEqual(testCase,warningId,'carConfig:legacyDOE');
    verifyEqual(testCase,stateAfter,stateBefore);
end
end

function testLegacyErrorIdentifiersAndBadOutputFormat(testCase)
base = vehicleConfig.defaultSetup();
[cP,aP,eP,DTp,Bp,tP] = legacyGroupsFromSetup(base);
verifyError(testCase,@() parameters_loop(cP,aP,eP,DTp,Bp,tP,"bad",[]), ...
    'parameters_loop:badSamplingType');
verifyError(testCase,@() parameters_loop(cP,aP,eP,DTp,Bp,tP,"Explicit",[]), ...
    'parameters_loop:badExplicitInput');
verifyError(testCase,@() carConfig("bad"),'carConfig:badRequest');
verifyError(testCase,@() carConfig("Explicit",table(),"legacy","nominal","other"), ...
    'carConfig:badOutputFormat');
verifyError(testCase,@() carConfig("Explicit",table(),"legacy","nominal",["cells","cases"]), ...
    'carConfig:badOutputFormat');
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

function root = repositoryRoot()
modelRoot = fileparts(fileparts(which('vehicleConfig.defaultSetup')));
root = fileparts(modelRoot);
end

function testCarConfigMissingMapErrorPrecedesDispatch(testCase)
[fixtureModelRoot,cleanup] = installMissingMapFixture(testCase);
verifyEqual(testCase, ...
    normalizeFixturePath(which('carConfig')), ...
    normalizeFixturePath(fullfile(fixtureModelRoot,'carConfig.m')));
verifyEqual(testCase, ...
    normalizeFixturePath(which('vehicleConfig.defaultSetup')), ...
    normalizeFixturePath(fullfile(fixtureModelRoot,'+vehicleConfig','defaultSetup.m')));

missingMapId = 'parameters_loop:aeroMapNotFound';
nonempty = table(155,'VariableNames',{'mass'});
empty = table(zeros(0,1),'VariableNames',{'mass'});
verifyError(testCase,@() carConfig(),missingMapId);
verifyError(testCase,@() carConfig("Explicit",nonempty),missingMapId);
verifyError(testCase,@() carConfig( ...
    "Explicit",empty,"legacy","nominal","cases"),missingMapId);
verifyError(testCase,@() carConfig("not-a-request"),missingMapId);
clear cleanup
end

function [fixtureModelRoot,cleanup] = installMissingMapFixture(testCase)
repositoryDefaultSetup = which('vehicleConfig.defaultSetup');
repositoryModelRoot = fileparts(fileparts(repositoryDefaultSetup));
originalPath = path;
originalFolder = pwd;
fixtureRoot = tempname(tempdir);
fixtureModelRoot = fullfile(fixtureRoot,'Full Car Models');
cleanup = onCleanup(@() restoreMissingMapFixture( ...
    originalPath,originalFolder,fixtureRoot));

[ok,message] = mkdir(fullfile(fixtureModelRoot,'+vehicleConfig'));
assert(ok,message);
[ok,message] = mkdir(fullfile(fixtureModelRoot,'carComponents'));
assert(ok,message);
copyfile(fullfile(repositoryModelRoot,'carConfig.m'),fixtureModelRoot);
copyfile(fullfile(repositoryModelRoot,'+vehicleConfig','defaultSetup.m'), ...
    fullfile(fixtureModelRoot,'+vehicleConfig'));
copyfile(fullfile(repositoryModelRoot,'carComponents', ...
    'Fx_combined_parameters_run38_30.mat'), ...
    fullfile(fixtureModelRoot,'carComponents'));
copyfile(fullfile(repositoryModelRoot,'carComponents', ...
    'Lapsim_Fy_combined_parameters_1965run15.mat'), ...
    fullfile(fixtureModelRoot,'carComponents'));

addpath(fixtureModelRoot,'-begin');
clear carConfig
eval('clear vehicleConfig.defaultSetup');
assert(strcmp(normalizeFixturePath(which('carConfig')), ...
    normalizeFixturePath(fullfile(fixtureModelRoot,'carConfig.m'))));
assert(strcmp(normalizeFixturePath(which('vehicleConfig.defaultSetup')), ...
    normalizeFixturePath(fullfile(fixtureModelRoot,'+vehicleConfig','defaultSetup.m'))));
assert(~isfile(fullfile(fixtureModelRoot,'aeromap_b26.csv')));
end

function restoreMissingMapFixture(originalPath,originalFolder,fixtureRoot)
path(originalPath);
cd(originalFolder);
clear carConfig
eval('clear vehicleConfig.defaultSetup');
[tempParent,~] = fileparts(fixtureRoot);
assert(strcmp(normalizeFixturePath(tempParent),normalizeFixturePath(tempdir)), ...
    'Refusing to remove the fixture outside MATLAB tempdir.');
if isfolder(fixtureRoot)
    rmdir(fixtureRoot,'s');
end
end

function normalized = normalizeFixturePath(pathValue)
normalized = lower(strrep(char(pathValue),'\','/'));
normalized = regexprep(normalized,'/+$','');
end
