function tests = test_rampSpeedReleaseManifest
tests = functiontests(localfunctions);
end

function setupOnce(~)
setup_paths;
end

function testManifestRecordsReleaseInventoryAndMetadata(testCase)
rootFolder = fileparts(fileparts(mfilename("fullpath")));
manifest = rampSpeed.releaseManifest(rootFolder);

verifyTrue(testCase,isdatetime(manifest.createdAt));
verifyTrue(testCase,isfield(manifest,"matlab"));
verifyTrue(testCase,strlength(string(manifest.matlab.release)) > 0);
verifyTrue(testCase,strlength(string(manifest.matlab.version)) > 0);
verifyEqual(testCase,manifest.schemaVersion,1);
verifyTrue(testCase,strlength(string(manifest.appVersion)) > 0);

files = string(manifest.files);
verifyTrue(testCase,any(endsWith(files,"RampSpeedApp.mlapp")));
verifyTrue(testCase,any(endsWith(files,"exportStudy.m")));
verifyTrue(testCase,any(endsWith(files,"releaseManifest.m")));
verifyTrue(testCase,any(contains(files,"tests")));
verifyTrue(testCase,any(endsWith(files,"setup_paths.m")));
verifyTrue(testCase,isfield(manifest,"requiredModelPaths"));
verifyTrue(testCase,any(contains(string(manifest.requiredModelPaths), ...
    "carComponents")));

verifyTrue(testCase,isfield(manifest,"git"));
verifyTrue(testCase,isfield(manifest.git,"available"));
verifyTrue(testCase,isfield(manifest.git,"commit"));
verifyTrue(testCase,isfield(manifest.git,"dirty"));
verifyTrue(testCase,islogical(manifest.git.available));
verifyTrue(testCase,islogical(manifest.git.dirty));

verifyTrue(testCase,isfield(manifest,"toolboxAvailability"));
verifyTrue(testCase,isstruct(manifest.toolboxAvailability));
verifyTrue(testCase,isfield(manifest.toolboxAvailability,"matlab"));
end

function testManifestSurvivesMissingGitMetadata(testCase)
[rootFolder,cleanup] = temporaryFolder(); %#ok<ASGLU>
manifest = rampSpeed.releaseManifest(rootFolder);

verifyFalse(testCase,manifest.git.available);
verifyEqual(testCase,string(manifest.git.commit),"");
verifyFalse(testCase,manifest.git.dirty);
verifyTrue(testCase,isdatetime(manifest.createdAt));
verifyTrue(testCase,isstruct(manifest.toolboxAvailability));
verifyTrue(testCase,isfield(manifest,"missingRequiredModelPaths"));
verifyTrue(testCase,any(contains(string(manifest.missingRequiredModelPaths), ...
    "carComponents")));
end

function [folder,cleanup] = temporaryFolder()
folder = tempname;
mkdir(folder);
cleanup = onCleanup(@()removeFolder(folder));
end

function removeFolder(folder)
if isfolder(folder)
    rmdir(folder,"s");
end
end