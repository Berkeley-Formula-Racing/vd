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
requiredArtifacts = ["RampSpeedApp.mlapp";"RampSpeedApp.prj"; ...
    "README_ramp_speed_app.md"];
verifyTrue(testCase,all(ismember(requiredArtifacts, ...
    string(manifest.missingRequiredModelPaths))));
end

function testProjectDescriptorUsesContainingFolderRoot(testCase)
sourceRoot = fileparts(fileparts(mfilename("fullpath")));
sourceProject = fullfile(sourceRoot,"RampSpeedApp.prj");
[temporaryRoot,cleanup] = temporaryFolder(); %#ok<ASGLU>
projectRoot = fullfile(temporaryRoot,"Full Car Models");
mkdir(projectRoot);
nativeProject = matlab.project.createProject(projectRoot);
nativeProject.close;
clear nativeProject;
targetProject = fullfile(projectRoot,"RampSpeedApp.prj");
copyfile(sourceProject,targetProject);

document = xmlread(targetProject);
rootNodes = document.getElementsByTagName("RootFolder");
rootText = string(char(rootNodes.item(0).getTextContent()));
verifyEqual(testCase,rootText,".");
loadedProject = matlab.project.loadProject(projectRoot);
cleanupLoaded = onCleanup(@()closeProject(loadedProject)); %#ok<NASGU>
verifyEqual(testCase,canonicalPath(loadedProject.RootFolder), ...
    canonicalPath(projectRoot));
closeProject(loadedProject);
clear loadedProject;
clear cleanupLoaded;
pause(0.1);
verifyEqual(testCase,canonicalPath(fullfile(projectRoot,char(rootText))), ...
    canonicalPath(projectRoot));
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
function closeProject(project)
if isvalid(project)
    project.close;
    clear project;
end
end

function path = canonicalPath(fileName)
path = string(char(java.io.File(char(fileName)).getCanonicalPath()));
end