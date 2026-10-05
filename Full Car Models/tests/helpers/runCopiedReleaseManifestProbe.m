%RUNCOPIEDRELEASEMANIFESTPROBE Verify declared release inputs in a fresh copy.
% Launch from a fresh MATLAB process with:
%   run(fullfile(modelRoot,'tests','helpers', ...
%       'runCopiedReleaseManifestProbe.m'))

sourceRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
addpath(sourceRoot);
manifest = rampSpeed.releaseManifest(sourceRoot);

runtimeDirectories = ["+rampSpeed"; "+vehicleConfig"; "carComponents"; ...
    "utilities"; "sweeps"];
runtimeFiles = ["setup_paths.m"; "carConfig.m"; ...
    "carConfigBaseline.m"; "runRampSpeedStudy.m"; "aeromap_b26.csv"];
runtimeRoots = [runtimeDirectories; runtimeFiles];
missingRoots = runtimeRoots(~ismember(runtimeRoots, ...
    string(manifest.requiredModelPaths)));
assert(isempty(missingRoots), ...
    "Release manifest omits runtime inputs: %s",strjoin(missingRoots,", "));

inventory = string(manifest.files);
for k = 1:numel(runtimeRoots)
    assert(any(inventory == runtimeRoots(k)), ...
        "Release file inventory omits required input %s.",runtimeRoots(k));
end

runStamp = char(datetime("now","Format","yyyyMMdd'T'HHmmss"));
artifactRoot = char(fullfile(sourceRoot,".superpowers","sdd", ...
    "2026-10-04-car-configuration-refactor","artifacts", ...
    "task7-copy-probe-" + string(runStamp)));
mkdir(artifactRoot);
releaseRoot = char(fullfile(artifactRoot,"release-copy"));
mkdir(releaseRoot);
diary(char(fullfile(artifactRoot,"copy-runtime-probe.log")));
diary on;
fprintf("Source root: %s\n",sourceRoot);
fprintf("Retained release copy: %s\n",releaseRoot);
fprintf("Declared runtime roots: %s\n",strjoin(runtimeRoots,", "));

selectedFiles = strings(0,1);
normalizedInventory = replace(inventory,string(filesep),"/");
normalizedDirectories = replace(runtimeDirectories,string(filesep),"/");
for k = 1:numel(runtimeDirectories)
    prefix = normalizedDirectories(k) + "/";
    inDirectory = startsWith(normalizedInventory,prefix);
    assert(any(inDirectory),"No declared runtime files found under %s.", ...
        runtimeDirectories(k));
    selectedFiles = [selectedFiles; inventory(inDirectory)]; %#ok<AGROW>
end
selectedFiles = unique([selectedFiles; runtimeFiles],"stable");
for k = 1:numel(selectedFiles)
    relative = selectedFiles(k);
    sourceFile = fullfile(sourceRoot,char(relative));
    assert(isfile(sourceFile),"Declared release input is absent: %s",relative);
    targetFile = fullfile(releaseRoot,char(relative));
    targetFolder = fileparts(targetFile);
    if ~isfolder(targetFolder), mkdir(targetFolder); end
    [ok,message] = copyfile(sourceFile,targetFile);
    assert(ok,"Could not copy %s: %s",relative,message);
end
fprintf("Copied %d declared files.\n",numel(selectedFiles));

% Reset the path after inventorying the source. The fresh
% process runtime below sees model code only from the copied inputs. This
% helper is launched in a new MATLAB process, so clearing classes/functions
% would only risk clearing script state and is unnecessary.
restoredefaultpath;
cd(releaseRoot);
run(fullfile(releaseRoot,'setup_paths.m'));
releaseCanonical = char(java.io.File(char(releaseRoot)).getCanonicalPath());
rootPrefix = [releaseCanonical filesep];
names = {"Aero","AeroMap","Car","Powertrain","Tire2","KTM450", ...
    "buildSingleRampCar","buildTireModel","vehicleConfig.defaultSetup", ...
    "vehicleConfig.buildCar","rampSpeed.loadRampAssets", ...
    "rampSpeed.releaseManifest","carConfig","carConfigBaseline", ...
    "runRampSpeedStudy"};
for k = 1:numel(names)
    resolved = which(names{k});
    assert(~isempty(resolved),"Could not resolve %s from the release copy.",names{k});
    resolved = char(java.io.File(resolved).getCanonicalPath());
    assert(strncmpi(resolved,rootPrefix,numel(rootPrefix)), ...
        "%s resolved outside the release copy: %s",names{k},resolved);
end
pathEntries = string(strsplit(path,pathsep));
blockedRoot = lower("C:" + string(filesep) + "VD");
assert(~any(startsWith(lower(pathEntries),blockedRoot)), ...
    "The copied runtime path still contains " + blockedRoot + ".");

assetPaths = [fullfile(releaseRoot,"aeromap_b26.csv"); ...
    fullfile(releaseRoot,"carComponents","Fx_combined_parameters_run38_30.mat"); ...
    fullfile(releaseRoot,"carComponents", ...
        "Lapsim_Fy_combined_parameters_1965run15.mat"); ...
    fullfile(releaseRoot,"carComponents","camber_models_fast.mat"); ...
    fullfile(releaseRoot,"carComponents","camberratiossmoothed.mat")];
for k = 1:numel(assetPaths)
    assert(isfile(assetPaths(k)),"Copied runtime asset is missing: %s",assetPaths(k));
end

qssSetup = vehicleConfig.defaultSetup("legacy","nominal");
assert(startsWith(string(qssSetup.aeroMap.path),string(rootPrefix), ...
    "IgnoreCase",true),"QSS map path is outside the release copy.");
qssAssets = vehicleConfig.loadAssets(qssSetup);
assert(startsWith(string(qssAssets.camberModelPath),string(rootPrefix), ...
    "IgnoreCase",true),"QSS camber data path is outside the release copy.");
[qssCar,qssDerived] = vehicleConfig.buildCar(qssSetup,"lap",qssAssets);
assert(isa(qssCar,"Car") && qssDerived.totalMass_kg > 0, ...
    "QSS baseline construction failed from the release copy.");

[rampCar,rampConfig] = carConfigBaseline();
assert(isa(rampCar,"Car"),"Frozen ramp baseline construction failed.");
rampMapPath = char(java.io.File(char(rampConfig.defaultMapPath)).getCanonicalPath());
assert(strncmpi(rampMapPath,rootPrefix,numel(rootPrefix)), ...
    "Frozen ramp baseline map resolved outside the release copy.");
fprintf("COPY_PROBE_OK: QSS lap baseline and frozen ramp baseline constructed.\n");
fprintf("QSS mass: %.6f kg; ramp baseline version: %s\n", ...
    qssDerived.totalMass_kg,string(rampConfig.baselineVersion));
fprintf("No solver was launched.\n");
diary off;