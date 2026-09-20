function manifest = releaseManifest(rootDirectory)
%RELEASEMANIFEST Describe the ramp-speed app release inputs and environment.

if nargin < 1 || isempty(rootDirectory)
    rootDirectory = fileparts(fileparts(mfilename("fullpath")));
end
rootDirectory = string(rootDirectory);
if ~isscalar(rootDirectory) || ~isfolder(rootDirectory)
    error("rampSpeed:invalidRootDirectory", ...
        "rootDirectory must be an existing folder.");
end
rootDirectory = absolutePath(rootDirectory);

manifest = struct();
manifest.rootDirectory = rootDirectory;
manifest.createdAt = datetime("now");
manifest.schemaVersion = 1;
manifest.appName = "RampSpeedApp";
manifest.appVersion = "1.0.0";
manifest.matlab = matlabMetadata();
[manifest.files,manifest.requiredModelPaths,manifest.missingRequiredModelPaths] = inventoryFiles(rootDirectory);
manifest.git = gitMetadata(rootDirectory);
manifest.toolboxAvailability = toolboxMetadata();
end

function metadata = matlabMetadata()
metadata = struct();
try
    metadata.release = string(version("-release"));
catch
    metadata.release = "unavailable";
end
try
    metadata.version = string(version);
catch
    metadata.version = "unavailable";
end
try
    metadata.architecture = string(computer("arch"));
catch
    metadata.architecture = "unavailable";
end
end

function [files,requiredPaths,missingPaths] = inventoryFiles(rootDirectory)
files = strings(0,1);
requiredPaths = strings(0,1);
missingPaths = strings(0,1);

required = ["setup_paths.m"; "+rampSpeed"; "tests"; "resources"; ...
    "carComponents"; "events"; "gg"; "CnAy"; "sweeps"; ...
    "utilities"; "testing"; "Phase Plane"; ...
    "Suspension Transient Model"; "Understeer Gradient"; ...
    "CamberCurveGen"; "CamberCurveGen2"];
for i = 1:numel(required)
    relative = required(i);
    candidate = fullfile(rootDirectory,relative);
    requiredPaths(end+1,1) = relative; %#ok<AGROW>
    if isfile(candidate)
        files(end+1,1) = relative; %#ok<AGROW>
    elseif isfolder(candidate)
        files(end+1,1) = relative; %#ok<AGROW>
        files = [files; collectFiles(rootDirectory,candidate)]; %#ok<AGROW>
    else
        missingPaths(end+1,1) = relative; %#ok<AGROW>
    end
end

releaseFiles = ["RampSpeedApp.mlapp";"RampSpeedApp.prj"; ...
    "README_ramp_speed_app.md"];
for i = 1:numel(releaseFiles)
    relative = releaseFiles(i);
    candidate = fullfile(rootDirectory,relative);
    requiredPaths(end+1,1) = relative; %#ok<AGROW>
    if isfile(candidate)
        files(end+1,1) = relative; %#ok<AGROW>
    else
        missingPaths(end+1,1) = relative; %#ok<AGROW>
    end
end
files = unique(sort(files));
requiredPaths = unique(requiredPaths,"stable");
missingPaths = unique(missingPaths,"stable");
end

function files = collectFiles(rootDirectory,folder)
files = strings(0,1);
entries = dir(folder);
for i = 1:numel(entries)
    name = string(entries(i).name);
    if name == "." || name == ".."
        continue
    end
    candidate = fullfile(folder,name);
    if entries(i).isdir
        files = [files; collectFiles(rootDirectory,candidate)]; %#ok<AGROW>
    else
        files(end+1,1) = relativePath(rootDirectory,candidate); %#ok<AGROW>
    end
end
end

function relative = relativePath(rootDirectory,fileName)
rootText = char(rootDirectory);
fileText = char(absolutePath(fileName));
prefixLength = length(rootText) + 1;
if length(fileText) >= prefixLength && ...
        startsWith(fileText,[rootText filesep])
    relative = string(fileText(prefixLength:end));
else
    relative = string(fileText);
end
relative = replace(relative,"/",filesep);
end

function metadata = gitMetadata(rootDirectory)
metadata = struct("available",false,"commit","","branch","", ...
    "dirty",false,"status","unavailable");
[code,commit] = runGit(rootDirectory,"rev-parse HEAD");
if code ~= 0 || strlength(firstLine(commit)) == 0
    return
end
metadata.available = true;
metadata.commit = firstLine(commit);
[branchCode,branch] = runGit(rootDirectory,"rev-parse --abbrev-ref HEAD");
if branchCode == 0
    metadata.branch = firstLine(branch);
end
[statusCode,statusText] = runGit(rootDirectory,"status --porcelain");
if statusCode == 0
    metadata.dirty = strlength(strtrim(string(statusText))) > 0;
    if metadata.dirty
        metadata.status = "dirty";
    else
        metadata.status = "clean";
    end
end
end

function [code,output] = runGit(rootDirectory,arguments)
quotedRoot = strrep(char(rootDirectory),'"','\"');
command = sprintf('git -C "%s" %s',quotedRoot,char(arguments));
[code,output] = system(command);
end

function line = firstLine(value)
lines = splitlines(strtrim(string(value)));
lines = lines(strlength(lines) > 0);
if isempty(lines)
    line = "";
else
    line = lines(1);
end
end

function metadata = toolboxMetadata()
metadata = struct();
metadata.matlab = safeLicenseTest("MATLAB");
metadata.optimization_toolbox = safeLicenseTest("Optimization_Toolbox");
metadata.simulink = safeLicenseTest("Simulink");
metadata.statistics_toolbox = safeLicenseTest("Statistics_Toolbox");
metadata.symbolic_math_toolbox = safeLicenseTest("Symbolic_Toolbox");
metadata.curve_fitting_toolbox = safeLicenseTest("Curve_Fitting_Toolbox");
metadata.parallel_computing_toolbox = safeLicenseTest("Distrib_Computing_Toolbox");
end

function available = safeLicenseTest(feature)

try
    available = logical(license("test",feature));
catch
    available = false;
end
if isempty(available)
    available = false;
end
available = available(1);
end

function path = absolutePath(fileName)
path = string(char(java.io.File(char(fileName)).getAbsolutePath()));
end