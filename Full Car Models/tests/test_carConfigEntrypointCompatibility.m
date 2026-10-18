%% test_carConfigEntrypointCompatibility
% Construction-only assertion script. Run this file with RUN, not RUNTESTS.
callerPath = path;
callerFolder = pwd;
compatibilityTestsDir = fileparts(mfilename('fullpath'));
compatibilitySourceManifest = fullfile(compatibilityTestsDir, ...
    'car_configuration_entrypoints.csv');
compatibilityOriginalManifestBytes = readManifestBytes( ...
    compatibilitySourceManifest);
runCompatibilityChecks(false);
assert(isequal(readManifestBytes(compatibilitySourceManifest), ...
    compatibilityOriginalManifestBytes), ...
    'Task6:sourceManifestChangedOnSuccess', ...
    'The construction checker changed the tracked entrypoint matrix on success.');
assert(strcmp(path,callerPath),'Task6:callerPathLeaked', ...
    'The compatibility script changed the caller MATLAB path.');
assert(strcmp(pwd,callerFolder),'Task6:callerFolderLeaked', ...
    'The compatibility script changed the caller current directory.');

try
    runCompatibilityChecks(true);
    error('Task6:cleanupProbeMissing', ...
        'The forced-failure cleanup probe did not fail as expected.');
catch ME
    assert(strcmp(ME.identifier,'Task6:cleanupFailureProbe'), ...
        'Task6:cleanupProbeWrongError','Unexpected cleanup probe error: %s',ME.message);
end
assert(isequal(readManifestBytes(compatibilitySourceManifest), ...
    compatibilityOriginalManifestBytes), ...
    'Task6:sourceManifestChangedOnFailure', ...
    'The construction checker changed the tracked entrypoint matrix on failure.');
assert(strcmp(path,callerPath),'Task6:failurePathLeaked', ...
    'Caller MATLAB path was not restored after the forced failure.');
assert(strcmp(pwd,callerFolder),'Task6:failureFolderLeaked', ...
    'Caller directory was not restored after the forced failure.');
fprintf('PASS caller path/cwd unchanged after successful and forced-failure run(filePath).\n');

function runCompatibilityChecks(forceFailure)
oldPath = path;
oldFolder = pwd;
foreignFolder = '';
rampFolder = '';
fixtureFolder = '';
if nargin < 1
    forceFailure = false;
end
try
testFile = mfilename('fullpath');
testsDir = fileparts(testFile);
fullCarModelsDir = fileparts(testsDir);
repoRoot = fileparts(fullCarModelsDir);
addpath(genpath(fullCarModelsDir));
setup_paths;
if forceFailure
    error('Task6:cleanupFailureProbe','Intentional cleanup failure-path probe.');
end
manifestPath = fullfile(testsDir,'car_configuration_entrypoints.csv');
evidenceRoot = fullfile(repoRoot,'artifacts','car_configuration_refactor', ...
    'final-fix-wave');
if ~isfolder(evidenceRoot), mkdir(evidenceRoot); end
evidenceDirectory = tempname(evidenceRoot);
if ~mkdir(evidenceDirectory)
    error('Task6:evidenceDirectory','Could not create isolated checker evidence folder.');
end
manifestCopyPath = fullfile(evidenceDirectory, ...
    'car_configuration_entrypoints.task6-checker-copy.csv');
if ~copyfile(manifestPath,manifestCopyPath)
    error('Task6:manifestCopy','Could not create the isolated matrix copy.');
end
manifestPath = manifestCopyPath;

report = runCarConfigurationChecks(manifestPath,evidenceDirectory);
assert(istable(report),'Task6:badReport','Construction driver must return a table.');
assert(height(report) == 34,'Task6:manifestRows', ...
    'Expected construction evidence for all 34 manifest rows.');
assert(all(ismember(report.status,["construction-pass","baseline-failure","blocked-input"])), ...
    'Task6:badStatus','The driver returned a status outside its construction tier.');
manifest = readtable(manifestPath,'TextType','string');
assert(height(manifest) == height(report));
assert(all(manifest.postStatus == report.status));
relativeEvidenceDirectory = replace(string(evidenceDirectory), ...
    [string(repoRoot) + filesep],"");
relativeEvidenceDirectory = replace(relativeEvidenceDirectory,'\','/');
assert(all(contains(manifest.evidencePath,relativeEvidenceDirectory)));
rampFunctionRow = manifest.path == "Full Car Models/runRampSpeedStudy.m";
rampScriptRow = manifest.path == "Full Car Models/entrypoints/runRampSpeedStudy.m";
assert(nnz(rampFunctionRow) == 1 && nnz(rampScriptRow) == 1);
assert(contains(manifest.blocker(rampFunctionRow), ...
    "Function resolution and load-only behavior only"));
assert(contains(manifest.blocker(rampScriptRow), ...
    "explicit run(scriptPath) invocation is unverified"));
assert(contains(manifest.blocker(rampScriptRow),"Task 8 numerical run"));
assert(~any(manifest.postStatus == "full-pass"), ...
    'Task6:overclaimed','Construction evidence cannot be labeled full-pass.');
fprintf('PASS manifest driver: %d rows, %d construction passes, %d baseline failures, %d blocked inputs.\n', ...
    height(report),nnz(report.status=="construction-pass"), ...
    nnz(report.status=="baseline-failure"),nnz(report.status=="blocked-input"));

% Bootstrap from a foreign working directory and verify canonical class paths.
foreignFolder = tempname;
mkdir(foreignFolder);
cd(foreignFolder);
bootstrapPath = fullfile(fullCarModelsDir,'entrypoints','bootstrap.m');
run(bootstrapPath);
assert(strcmp(normalizePath(pwd),normalizePath(repoRoot)), ...
    'Task6:bootstrapCwd','Bootstrap did not select the repository root.');
for className = ["Car","Tire2","vehicleConfig.buildCar"]
    resolved = which(char(className));
    assert(~isempty(resolved),'Task6:missingClass', ...
        'Could not resolve %s after bootstrap.',className);
    assert(startsWith(normalizePath(resolved), ...
        [normalizePath(fullCarModelsDir) '\']), ...
        'Task6:shadowedClass','%s resolved outside Full Car Models: %s',className,resolved);
end
fprintf('PASS foreign-directory bootstrap and Car/Tire2/vehicleConfig resolution.\n');

% Keep the same-named function and script available through their distinct APIs.
functionPath = which('runRampSpeedStudy');
scriptPath = fullfile(fullCarModelsDir,'entrypoints','runRampSpeedStudy.m');
assert(strcmp(normalizePath(functionPath), ...
    normalizePath(fullfile(fullCarModelsDir,'runRampSpeedStudy.m'))), ...
    'Task6:rampFunctionResolution', ...
    'runRampSpeedStudy must resolve to the options function under Full Car Models.');
assert(isfile(scriptPath),'Task6:rampScriptMissing','Legacy ramp study script is missing.');
allRampPaths = which('runRampSpeedStudy','-all');
if ischar(allRampPaths), allRampPaths = cellstr(allRampPaths); end
if isstring(allRampPaths), allRampPaths = cellstr(allRampPaths(:)); end
normalizedRampPaths = cellfun(@normalizePath,allRampPaths,'UniformOutput',false);
assert(any(strcmp(normalizedRampPaths,normalizePath(scriptPath))), ...
    'Task6:rampScriptResolution','The explicit legacy script path is not available on MATLAB path.');
assert(nargin('runRampSpeedStudy') == 1 && nargout('runRampSpeedStudy') == 4, ...
    'Task6:rampFunctionContract','The options function signature changed.');

% A load-only call proves the options function is callable without starting a solve.
rampFolder = tempname;
mkdir(rampFolder);
study = rampSpeed.makeStudy("task6-entrypoint");
rampCache = fullfile(rampFolder,'study.mat');
save(rampCache,'study','-v7.3');
rampOptions = struct('loadFromPrev',true,'cacheName',rampCache, ...
    'makeFigures',false,'figurePath',"");
[loadedStudy,loadedRuns,figures,events] = runRampSpeedStudy(rampOptions);
assert(loadedStudy.schemaVersion == 1 && isempty(loadedRuns));
assert(isempty(figures) && isempty(events));
% The script path is present, but its invocation has not been tested. Running it
% starts its configured study and remains deferred to bounded Task 8.
fprintf('PASS root ramp function load-only API and legacy script path presence; script invocation unverified, deferred to Task 8.\n');

% Preserve cells by default and lap/acceleration order in named cases.
Q = table([155;170],'VariableNames',{'mass'});
[carCell,eventParams,designTable] = carConfig("Explicit",Q,"legacy","nominal","cells");
[cases,eventParamsCases,caseDesignTable] = carConfig( ...
    "Explicit",Q,"legacy","nominal","cases");
assert(isequal(size(carCell),[2 2]) && isequal(size(cases),[2 1]));
assert(isequal(designTable,caseDesignTable));
assert(isequal(eventParams,eventParamsCases));
for row = 1:2
    assert(cases(row).designRow == row);
    assert(cases(row).lapCar.M == carCell{row,1}.M);
    assert(cases(row).accelCar.M == carCell{row,2}.M);
    assert(cases(row).lapCar.M > cases(row).accelCar.M);
    assert(isa(cases(row).lapCar,'Car') && isa(cases(row).accelCar,'Car'));
end
fprintf('PASS legacy cells and named case role ordering.\n');

% Build independent setup variants and preserve distinct static/map aero policies.
setup = vehicleConfig.defaultSetup();
baseCases = vehicleConfig.buildCases(setup);
baseCell = vehicleConfig.toLegacyCells(baseCases);
variant = setup;
variant.parameters.wheelbase = 1.60;
variant.parameters.final_drive = 3.2;
variantCases = vehicleConfig.buildCases([setup;variant]);
assert(variantCases(1).lapCar.W_b ~= variantCases(2).lapCar.W_b);
assert(variantCases(1).lapCar.powertrain.final_drive ~= ...
    variantCases(2).lapCar.powertrain.final_drive);
assert(setup.parameters.wheelbase ~= variant.parameters.wheelbase);

[staticCars,staticPlan] = aeroStarCases(baseCell, ...
    struct('ClA',[0.90 1.10],'CdA',[],'CoP',[]));
assert(height(staticPlan) == 3 && staticPlan.param(1) == "baseline");
assert(baseCell{1,1}.aero.hasMap());
assert(~staticCars{2,1}.aero.hasMap(), ...
    'Task6:staticAeroPolicy','Static sensitivity should retain its current static Aero replacement.');
assert(baseCell{1,1}.aero.cla == baseCases(1).lapCar.aero.cla);
for row = 1:size(staticCars,1)
    accel = staticCars{row,2};
    assert(accel.M == baseCell{1,2}.M && ...
        accel.aero.cla == baseCell{1,2}.aero.cla && ...
        accel.aero.cda == baseCell{1,2}.aero.cda && ...
        accel.aero.hasMap() == baseCell{1,2}.aero.hasMap());
end

[mapCars,mapPlan] = aeroMapStarCases(baseCell, ...
    struct('FrontRideHeightIn',0.125,'RearRideHeightIn',0.125, ...
    'ClAScale',[0.90 1.10],'CdAScale',1.10,'CoPOffset',0.03));
mapRow = find(mapPlan.swept_parameter == "ClAScale" & ...
    abs(mapPlan.swept_value-1.10) < 1e-12,1);
assert(~isempty(mapRow));
assert(mapCars{mapRow,1}.aero.hasMap());
assert(mapCars{mapRow,1}.aero.map.claScale == 1.10, ...
    'Task6:mapCorrections','Mapped sensitivity lost AeroMap.withCorrections values.');
assert(baseCell{1,1}.aero.map.claScale == 1 && ...
    baseCell{1,1}.aero.map.cdaScale == 1 && baseCell{1,1}.aero.map.copOffset == 0, ...
    'Task6:mutatedBaseline','Mapped variants changed the baseline AeroMap.');
for row = 1:size(mapCars,1)
    accel = mapCars{row,2};
    assert(accel.M == baseCell{1,2}.M && ...
        accel.aero.cla == baseCell{1,2}.aero.cla && ...
        accel.aero.cda == baseCell{1,2}.aero.cda && ...
        accel.aero.hasMap() == baseCell{1,2}.aero.hasMap());
end
fprintf('PASS setup independence and separate static/map sweep policies.\n');

% A synthetic legacy MAT fixture exercises cache, checkpoint and offline analysis.
[carCell,plan] = aeroStarCases(baseCell, ...
    struct('ClA',[0.90 0.95 1.05 1.10],'CdA',[],'CoP',[]));
for row = 1:size(carCell,1)
    car = carCell{row,1};
    car.ggPoints = [row 2;3 4];
    car.ss_info = [row 0.5 0.6];
    carCell{row,1} = car;
end
designTable = table([155;160;165;170;175],'VariableNames',{'mass'});
metricTable = table((1:5)',true(5,1),repmat("",5,1), ...
    [46;47;48;49;50], ...
    'VariableNames',{'case_index','valid','error_message','t_autox'});
times = table((1:5)',[46;47;48;49;50], ...
    'VariableNames',{'case_index','t_autox'});
traces = repmat(struct('time_s',[0;1],'speed_mps',[10;11]),5,1);
fixtureFolder = tempname;
mkdir(fixtureFolder);
legacyPath = fullfile(fixtureFolder,'legacy-cache.mat');
save(legacyPath,'carCell','plan','designTable','metricTable','times','traces','-v7.3');
variables = whos('-file',legacyPath);
savedNames = sort(string({variables.name}));
assert(isequal(savedNames,sort(["carCell","plan","designTable", ...
    "metricTable","times","traces"])), ...
    'Task6:legacyMatSchema','Legacy cache MAT variable names changed.');
[cachedCars,cachedTimes,cachedTraces,haveGG,haveTimes] = ...
    loadAeroCache(legacyPath,plan,baseCell{1,1});
assert(haveGG && haveTimes && isequal(cachedTimes,times) && isequal(cachedTraces,traces));
assert(isequal(cachedCars{1,1}.ggPoints,carCell{1,1}.ggPoints));

resolvedStudy = struct('signature',"task6-legacy-fixture", ...
    'mode',"screening",'maxCases',5,'numWorkers',0,'batchSize',5);
checkpointState = struct();
checkpointState.resolvedStudy = resolvedStudy;
checkpointState.carCell = carCell;
checkpointState.plan = plan;
checkpointState.designTable = designTable;
checkpointPath = fullfile(fixtureFolder,'DOE_checkpoint.mat');
saveDOECheckpoint(checkpointPath,checkpointState);
[loadedCheckpoint,resumeInfo] = loadDOECheckpoint(checkpointPath,resolvedStudy);
assert(isempty(resumeInfo.changed));
assert(isequal(loadedCheckpoint.carCell{1,1}.ss_info,carCell{1,1}.ss_info));
checkpointVariables = whos('-file',checkpointPath);
assert(isequal(string({checkpointVariables.name}),"checkpoint"));

analysis = doeAnalyze(legacyPath,struct('responsesWanted',{{'t_autox'}}, ...
    'predictors',{{'mass'}},'plots',strings(0,1),'visible','off', ...
    'fitGaussianProcesses',false,'sobolSamples',64,'sobolBootstrapSamples',8));
assert(analysis.settings.usedCachedMetrics);
assert(analysis.quality.total == 5 && analysis.quality.valid == 5);

named = struct('id',"legacy-fixture",'label',"legacy fixture", ...
    'designRow',1,'setup',setup,'lapCar',carCell{1,1}, ...
    'accelCar',carCell{1,2});
packed = vehicleConfig.toLegacyCells(named);
assert(isequal(packed{1,1}.ggPoints,carCell{1,1}.ggPoints));
assert(isequal(packed{1,1}.ss_info,carCell{1,1}.ss_info));
fprintf('PASS legacy MAT/cache/checkpoint/analysis and named-case result preservation.\n');

% Run the setup example itself: it must construct named cars and perform no solve.
run(fullfile(fullCarModelsDir,'entrypoints','car_setup_example.m'));
assert(isa(car,'Car') && isa(lapCar,'Car'));
assert(numel(cases) == 2 && cases(2).lapCar.M > cases(1).lapCar.M);

% Keep adaptive-grip test hooks and the one-zero/zero baseline contract.
adaptiveSuite = matlab.unittest.TestSuite.fromFile( ...
    fullfile(testsDir,'test_adaptive_grip_sweep.m'));
runner = matlab.unittest.TestRunner.withTextOutput;
adaptiveResults = runner.run(adaptiveSuite);
assert(all([adaptiveResults.Passed]),'Task6:adaptiveGripFailed', ...
    'The injected adaptive-grip contract test failed.');
fprintf('PASS runnable setup example and adaptive-grip hooks (%d tests).\n', ...
    numel(adaptiveResults));

catch ME
    restoreSessionState(oldPath,oldFolder, ...
        {foreignFolder,rampFolder,fixtureFolder});
    rethrow(ME)
end
restoreSessionState(oldPath,oldFolder, ...
    {foreignFolder,rampFolder,fixtureFolder});
fprintf('Task 6 compatibility assertions completed without solver or Events2 construction.\n');
end
function value = normalizePath(value)
value = lower(strrep(char(string(value)),'/','\'));
while numel(value) > 3 && endsWith(value,'\')
    value(end) = [];
end
end

function bytes = readManifestBytes(filePath)
fid = fopen(filePath,'rb');
if fid < 0
    error('Task6:manifestRead','Could not read entrypoint matrix: %s',filePath);
end
fileCleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
bytes = fread(fid,Inf,'*uint8').';
end


function restoreSessionState(oldPath,oldFolder,temporaryFolders)
path(oldPath);
cd(oldFolder);
for i = 1:numel(temporaryFolders)
    removeFolder(temporaryFolders{i});
end
end
function removeFolder(folder)
if isfolder(folder)
    rmdir(folder,'s');
end
end
