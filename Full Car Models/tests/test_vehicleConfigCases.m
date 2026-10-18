function tests = test_vehicleConfigCases
%TEST_VEHICLECONFIGCASES Verify named case identity, shape, and lossless packing.
tests = functiontests(localfunctions);
end

function testNamedCasesPackRowsAndPreserveSolvedState(testCase)
base = vehicleConfig.defaultSetup();
first = base;
first.parameters.mass = 155;
first.label = "light";
second = base;
second.parameters.mass = 170;
second.label = "heavy";

namedCases = vehicleConfig.buildCases([first;second]);
verifySize(testCase,namedCases,[2 1]);
verifyEqual(testCase,string({namedCases.id}).',["case-000001";"case-000002"]);
verifyEqual(testCase,[namedCases.designRow].',[1;2]);
verifyEqual(testCase,string({namedCases.label}).',["light";"heavy"]);
verifyEqual(testCase,[namedCases(1).setup.parameters.mass; ...
    namedCases(2).setup.parameters.mass],[155;170]);
verifyEqual(testCase,[namedCases(1).lapCar.M;namedCases(2).lapCar.M],[219;234]);
verifyEqual(testCase,[namedCases(1).accelCar.M;namedCases(2).accelCar.M],[214;229]);

[eventResult,eventAlias] = savedEvents2(testCase);
verifyTrue(testCase,eventResult == eventAlias);
namedCases(1).lapCar.ggPoints = [1 2;3 4];
namedCases(1).lapCar.ss_info = [5 6 7];
namedCases(1).lapCar.comp = eventResult;
packed = vehicleConfig.toLegacyCells(namedCases);
verifySize(testCase,packed,[2 2]);
verifyEqual(testCase,packed{1,1}.ggPoints,[1 2;3 4]);
verifyEqual(testCase,packed{1,1}.ss_info,[5 6 7]);
verifyClass(testCase,packed{1,1}.comp,'Events2');
verifyTrue(testCase,packed{1,1}.comp == eventResult);
verifyEmpty(testCase,namedCases(2).lapCar.comp);

emptyCases = vehicleConfig.buildCases(repmat(base,0,1));
verifySize(testCase,emptyCases,[0 1]);
verifySize(testCase,vehicleConfig.toLegacyCells(emptyCases),[0 2]);
end

function testDistinctAssetIdentitiesDoNotShareMapOrTireModel(testCase)
legacy = vehicleConfig.defaultSetup("legacy","nominal");
lc0Low = vehicleConfig.defaultSetup("lc0_nd","low");
lc0High = vehicleConfig.defaultSetup("lc0_nd","high");
noMap = legacy;
noMap.aeroMap.enabled = false;
noMap.aeroMap.path = fullfile(repositoryRoot(),'Full Car Models','missing-map.csv');

namedCases = vehicleConfig.buildCases([legacy;lc0Low;lc0High;noMap]);
verifyClass(testCase,namedCases(1).lapCar.tire,'Tire2');
verifyClass(testCase,namedCases(2).lapCar.tire,'NondimensionalTire');
verifyClass(testCase,namedCases(3).lapCar.tire,'NondimensionalTire');
verifyEqual(testCase,namedCases(2).lapCar.tire.uncertainty_mode,"low");
verifyEqual(testCase,namedCases(3).lapCar.tire.uncertainty_mode,"high");
verifyTrue(testCase,namedCases(1).lapCar.aero.hasMap());
verifyTrue(testCase,namedCases(2).lapCar.aero.hasMap());
verifyFalse(testCase,namedCases(4).lapCar.aero.hasMap());
verifyFalse(testCase,namedCases(1).accelCar.aero.hasMap());
verifyFalse(testCase,namedCases(2).accelCar.aero.hasMap());
end

function testEmptySingletonAndTwoCaseMatRoundTrips(testCase)
base = vehicleConfig.defaultSetup();
first = base;
first.parameters.mass = 155;
second = base;
second.parameters.mass = 170;
emptyCases = vehicleConfig.buildCases(repmat(base,0,1));
twoCases = vehicleConfig.buildCases([first;second]);
[eventObject,eventAlias] = savedEvents2(testCase);
twoCases(1).lapCar.comp = eventObject;
singleCase = twoCases(1);
emptyCells = vehicleConfig.toLegacyCells(emptyCases);
singleCells = vehicleConfig.toLegacyCells(singleCase);
twoCells = vehicleConfig.toLegacyCells(twoCases);

folder = tempname;
mkdir(folder);
cleanup = onCleanup(@() removeTempFolder(folder));
matPath = fullfile(folder,'case-roundtrip.mat');
save(matPath,'emptyCases','singleCase','twoCases', ...
    'emptyCells','singleCells','twoCells','eventObject','eventAlias');
loaded = load(matPath);
verifySize(testCase,loaded.emptyCases,[0 1]);
verifySize(testCase,loaded.emptyCells,[0 2]);
verifySize(testCase,loaded.singleCase,[1 1]);
verifySize(testCase,loaded.singleCells,[1 2]);
verifySize(testCase,loaded.twoCases,[2 1]);
verifySize(testCase,loaded.twoCells,[2 2]);
verifyEqual(testCase,loaded.twoCases(1).id,"case-000001");
verifyEqual(testCase,loaded.twoCases(2).setup.parameters.mass,170);
verifyClass(testCase,loaded.eventObject,'Events2');
verifyEqual(testCase,loaded.eventObject.times,eventObject.times);
verifyEqual(testCase,loaded.twoCases(1).lapCar.comp.times,eventObject.times);
verifyTrue(testCase,loaded.eventObject == loaded.eventAlias);
verifyTrue(testCase,loaded.twoCases(1).lapCar.comp == loaded.eventObject);
verifyTrue(testCase,loaded.twoCells{1,1}.comp == loaded.eventObject);
for row = 1:2
    actualLap = carConfigurationSnapshot(loaded.twoCells{row,1});
    expectedLap = carConfigurationSnapshot(twoCells{row,1});
    actualLap.values.comp = [];
    expectedLap.values.comp = [];
    verifyEqual(testCase,actualLap,expectedLap);
    actualAccel = carConfigurationSnapshot(loaded.twoCells{row,2});
    expectedAccel = carConfigurationSnapshot(twoCells{row,2});
    actualAccel.values.comp = [];
    expectedAccel.values.comp = [];
    verifyEqual(testCase,actualAccel,expectedAccel);
end
clear cleanup
end

function testTask1DefaultExplicitAndTireSelectionSnapshots(testCase)
reference = constructionReference(testCase);
base = vehicleConfig.defaultSetup();
defaultCase = vehicleConfig.buildCases(base);
verifyEqual(testCase,carConfigurationSnapshot(defaultCase(1).lapCar), ...
    reference.qss.default.lapCar);
verifyEqual(testCase,carConfigurationSnapshot(defaultCase(1).accelCar), ...
    reference.qss.default.accelCar);

explicit = reference.qss.explicit;
explicitSetups = repmat(base,height(explicit.designTable),1);
for row = 1:height(explicit.designTable)
    explicitSetups(row) = applyDesignRow(base,explicit.designTable,row);
end
explicitCases = vehicleConfig.buildCases(explicitSetups);
for row = 1:height(explicit.designTable)
    verifyEqual(testCase,carConfigurationSnapshot(explicitCases(row).lapCar), ...
        explicit.lapCars{row});
    verifyEqual(testCase,carConfigurationSnapshot(explicitCases(row).accelCar), ...
        explicit.accelCars{row});
end

selectionNames = fieldnames(reference.qss.tireSelections);
selectionSetups = repmat(base,numel(selectionNames),1);
for k = 1:numel(selectionNames)
    item = reference.qss.tireSelections.(selectionNames{k});
    selectionSetups(k) = vehicleConfig.defaultSetup(item.model,item.uncertainty);
end
selectionCases = vehicleConfig.buildCases(selectionSetups);
for k = 1:numel(selectionNames)
    item = reference.qss.tireSelections.(selectionNames{k});
    verifyEqual(testCase,carConfigurationSnapshot(selectionCases(k).lapCar), ...
        item.lapCar);
    verifyEqual(testCase,carConfigurationSnapshot(selectionCases(k).accelCar), ...
        item.accelCar);
end
end

function setup = applyDesignRow(setup,designTable,row)
names = designTable.Properties.VariableNames;
for k = 1:numel(names)
    setup.parameters.(names{k}) = designTable{row,names{k}};
end
end

function reference = constructionReference(testCase)
path = fullfile(repositoryRoot(),'artifacts', ...
    'car_configuration_refactor','2026-10-04-reference', ...
    'construction-reference.mat');
if strcmp(getenv('VD_TASK4_REQUIRE_REFERENCE'),'1') && ~isfile(path)
    error('test_vehicleConfigCases:missingReference', ...
        'VD_TASK4_REQUIRE_REFERENCE is set but the Task 1 snapshot is missing.');
end
assumeTrue(testCase,isfile(path), ...
    'The Task 1 construction reference is not present in this clone.');
loaded = load(path,'reference');
reference = loaded.reference;
end

function root = repositoryRoot()
modelRoot = fileparts(fileparts(which('vehicleConfig.defaultSetup')));
root = fileparts(modelRoot);
end

function removeTempFolder(path)
if isfolder(path)
    rmdir(path,'s');
end
end




function [eventObject,eventAlias] = savedEvents2(testCase)
path = fullfile(repositoryRoot(),'artifacts', ...
    'car_configuration_refactor','2026-10-04-reference-production-28mps', ...
    'numerical-reference.mat');
if strcmp(getenv('VD_TASK3_REQUIRE_NUMERICAL_REFERENCE'),'1') && ~isfile(path)
    error('test_vehicleConfigCases:missingNumericalReference', ...
        'VD_TASK3_REQUIRE_NUMERICAL_REFERENCE is set but the numerical reference is missing.');
end
assumeTrue(testCase,isfile(path), ...
    'The Task 1 numerical reference is not present in this clone.');
loaded = load(path,'numerical');
eventObject = loaded.numerical.qss.eventObject;
eventAlias = loaded.numerical.qss.eventAlias;
verifyClass(testCase,eventObject,'Events2');
verifyClass(testCase,eventAlias,'Events2');
end
