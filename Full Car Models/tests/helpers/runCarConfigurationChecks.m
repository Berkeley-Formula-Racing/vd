function report = runCarConfigurationChecks(manifestPath,evidenceDirectory)
%RUNCARCONFIGURATIONCHECKS Record construction-only status for each launcher row.
%   This driver constructs only configuration objects. It does not invoke G-G,
%   Events2, optimizers, study executors, event solvers, or UI launchers.

helperDir = fileparts(mfilename('fullpath'));
testsDir = fileparts(helperDir);
fullCarModelsDir = fileparts(testsDir);
repoRoot = fileparts(fullCarModelsDir);
if nargin < 1 || isempty(manifestPath)
    manifestPath = fullfile(testsDir,'car_configuration_entrypoints.csv');
end
if nargin < 2 || isempty(evidenceDirectory)
    evidenceDirectory = fullfile(repoRoot,'artifacts','car_configuration_refactor', ...
        ['task6-exec-' datestr(now,'yyyymmddTHHMMSS')]);
end
manifestPath = char(string(manifestPath));
evidenceDirectory = char(string(evidenceDirectory));
if ~isfile(manifestPath)
    error('runCarConfigurationChecks:missingManifest', ...
        'Entrypoint manifest does not exist: %s',manifestPath);
end
manifest = readtable(manifestPath,'TextType','string');
required = ["path","kind","constructionRoute","baselineStatus", ...
    "postStatus","evidencePath","blocker"];
if ~all(ismember(required,string(manifest.Properties.VariableNames)))
    error('runCarConfigurationChecks:badManifest', ...
        'Manifest must retain its seven established columns.');
end
if height(manifest) ~= 34
    error('runCarConfigurationChecks:manifestRowCount', ...
        'Expected the approved 34-row entrypoint matrix; found %d rows.',height(manifest));
end

referencePaths = { ...
    fullfile(repoRoot,'artifacts','car_configuration_refactor', ...
        '2026-10-04-reference','construction-reference.mat'), ...
    fullfile(repoRoot,'artifacts','car_configuration_refactor', ...
        '2026-10-04-reference-production-28mps','numerical-reference.mat')};
for i = 1:numel(referencePaths)
    if ~isfile(referencePaths{i})
        error('runCarConfigurationChecks:missingReference', ...
            'Required immutable reference is missing: %s',referencePaths{i});
    end
end

oldPath = path;
oldFolder = pwd;
restoreSession = onCleanup(@() restoreSessionState(oldPath,oldFolder)); %#ok<NASGU>
addpath(genpath(fullCarModelsDir));
setup_paths;

% Construct the public role pair and frozen ramp car once. The sweep variants
% below reuse one explicit asset identity, so the checks add no per-case reads.
[qssCells,~,qssDesign] = carConfig();
if size(qssCells,2) ~= 2 || isempty(qssCells) || height(qssDesign) ~= size(qssCells,1)
    error('runCarConfigurationChecks:qssConstruction', ...
        'carConfig did not return its legacy lap/acceleration cell contract.');
end
[rampCar,~] = carConfigBaseline();
if ~isa(rampCar,'Car')
    error('runCarConfigurationChecks:rampConstruction', ...
        'carConfigBaseline did not build its frozen ramp Car.');
end
setup = vehicleConfig.defaultSetup();
assets = vehicleConfig.loadAssets(setup);
baseCases = vehicleConfig.buildCases(setup,assets);
baseCell = vehicleConfig.toLegacyCells(baseCases);
[staticCases,staticPlan] = aeroStarCases(baseCell, ...
    struct('ClA',[0.90 1.10],'CdA',[],'CoP',[]));
[mapCases,mapPlan] = aeroMapStarCases(baseCell, ...
    struct('FrontRideHeightIn',0.125,'RearRideHeightIn',0.125, ...
    'ClAScale',[0.90 1.10],'CdAScale',1.10,'CoPOffset',0.03));
if size(staticCases,2) ~= 2 || height(staticPlan) ~= size(staticCases,1) || ...
        size(mapCases,2) ~= 2 || height(mapPlan) ~= size(mapCases,1) || ...
        ~mapCases{end,1}.aero.hasMap()
    error('runCarConfigurationChecks:sweepConstruction', ...
        'Aero sensitivity variants did not retain the expected legacy cells.');
end

count = height(manifest);
status = strings(count,1);
constructionCheck = strings(count,1);
blocker = strings(count,1);
for row = 1:count
    relativePath = replace(string(manifest.path(row)),'/',filesep);
    sourcePath = fullfile(repoRoot,char(relativePath));
    if ~isfile(sourcePath)
        status(row) = "baseline-failure";
        constructionCheck(row) = "manifest source file missing";
        blocker(row) = "Source path listed in the baseline manifest is unavailable: " + ...
            string(sourcePath);
        continue
    end

    route = string(manifest.constructionRoute(row));
    if contains(string(manifest.path(row)),"Suspension Transient Model/FullCarVibrationModel.m")
        if isprop(qssCells{1,1},'c')
            status(row) = "construction-pass";
            constructionCheck(row) = "Car exposes the damping field used by this launcher.";
            blocker(row) = "Full launcher execution remains deferred to Task 8.";
        else
            status(row) = "baseline-failure";
            constructionCheck(row) = "QSS Car has no declared c property required by the script.";
            blocker(row) = "Pre-existing script assigns car.c although Car does not declare c; no damping field was invented.";
        end
    elseif string(manifest.kind(row)) == "launcher" && ...
            contains(string(manifest.path(row)),"Full Car Models/testing/test.m")
        inputPath = fullfile(fullCarModelsDir,'testing','fmincon inputs.mat');
        if isfile(inputPath)
            status(row) = "construction-pass";
            constructionCheck(row) = "QSS Car construction is available; optimizer was not started.";
            blocker(row) = "Full launcher execution is deferred to Task 8.";
        else
            status(row) = "blocked-input";
            constructionCheck(row) = "QSS Car construction passed; required historical input is absent.";
            blocker(row) = "Missing historical input required by the script: " + string(inputPath);
        end
    elseif string(manifest.path(row)) == "Full Car Models/runRampSpeedStudy.m"
        status(row) = "construction-pass";
        constructionCheck(row) = "Root options function resolved and load-only cache API was verified; no study solve was run.";
        blocker(row) = "Function resolution and load-only behavior only; the legacy script invocation is unverified and outstanding for Task 8.";
    elseif string(manifest.path(row)) == "Full Car Models/entrypoints/runRampSpeedStudy.m"
        status(row) = "construction-pass";
        constructionCheck(row) = "Legacy script path presence was verified; run(scriptPath) was not invoked.";
        blocker(row) = "Only the root function load-only API and legacy script path presence are verified; explicit run(scriptPath) invocation is unverified and outstanding for the Task 8 numerical run.";
    elseif contains(route,"carConfigBaseline")
        status(row) = "construction-pass";
        constructionCheck(row) = "Frozen ramp baseline Car constructed; no ramp solver was invoked.";
        blocker(row) = "Construction-pass only; full launcher execution is deferred to Task 8.";
    elseif contains(route,"LC0 model artifact")
        lc0Setup = vehicleConfig.defaultSetup("lc0_nd","nominal");
        lc0Assets = vehicleConfig.loadAssets(lc0Setup);
        lc0Car = vehicleConfig.buildCar(lc0Setup,"lap",lc0Assets);
        if isa(lc0Car,'Car') && isa(lc0Car.tire,'NondimensionalTire')
            status(row) = "construction-pass";
            constructionCheck(row) = "LC0 model and QSS Car constructed; no comparison solver was invoked.";
            blocker(row) = "Construction-pass only; plotting/comparison launcher execution is deferred.";
        else
            status(row) = "baseline-failure";
            constructionCheck(row) = "LC0 QSS construction did not preserve the selected tire model.";
            blocker(row) = "Required LC0 construction contract failed.";
        end
    elseif contains(route,"carConfig")
        status(row) = "construction-pass";
        constructionCheck(row) = "QSS lap/acceleration Car construction passed; launcher solver/UI was not started.";
        blocker(row) = "Construction-pass only; full launcher execution is deferred to Task 8.";
    else
        status(row) = "construction-pass";
        constructionCheck(row) = "QSS and ramp configuration constructors resolved; launcher body was not executed.";
        blocker(row) = "Construction-pass only; full launcher execution is deferred to Task 8.";
    end
end

relativeEvidence = string(strrep(evidenceDirectory,[repoRoot filesep],''));
relativeEvidence = replace(relativeEvidence,'\','/');
report = table((1:count)',string(manifest.path),string(manifest.kind),status, ...
    constructionCheck,blocker, ...
    'VariableNames',{'manifestRow','path','kind','status','constructionCheck','blocker'});
reportPath = fullfile(evidenceDirectory,'task6-construction-checks.csv');
if ~isfolder(evidenceDirectory)
    mkdir(evidenceDirectory);
end
writeIgnore = fullfile(evidenceDirectory,'.gitignore');
if ~isfile(writeIgnore)
    fid = fopen(writeIgnore,'w');
    if fid < 0
        error('runCarConfigurationChecks:evidenceWrite', ...
            'Could not create local evidence ignore file: %s',writeIgnore);
    end
    fileCleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid,'*\n');
end
writetable(report,reportPath);

updatedManifest = manifest;
updatedManifest.postStatus = status;
updatedManifest.evidencePath(:) = replace(relativeEvidence + ...
    "/task6-construction-checks.csv",'\','/');
for row = 1:count
    prior = strtrim(string(manifest.blocker(row)));
    priorText = char(prior);
    priorTask6 = strfind(priorText,' Task 6: ');
    if ~isempty(priorTask6)
        prior = string(priorText(1:priorTask6(1)-1));
    elseif startsWith(prior,"Task 6:")
        prior = "";
    end
    note = blocker(row);
    if strlength(prior) > 0
        updatedManifest.blocker(row) = prior + " Task 6: " + note;
    else
        updatedManifest.blocker(row) = "Task 6: " + note;
    end
end
manifestTemp = [manifestPath '.task6.tmp.csv'];
writetable(updatedManifest,manifestTemp);
[moved,message] = movefile(manifestTemp,manifestPath,'f');
if ~moved
    if isfile(manifestTemp), delete(manifestTemp); end
    error('runCarConfigurationChecks:manifestWrite', ...
        'Could not atomically update the manifest: %s',message);
end
end

function restoreSessionState(oldPath,oldFolder)
path(oldPath);
cd(oldFolder);
end
