function tests = test_rampSpeedExport
tests = functiontests(localfunctions);
end

function setupOnce(~)
setup_paths;
testsFolder = fileparts(mfilename("fullpath"));
addpath(fullfile(testsFolder,"helpers"));
end

function testExportsCanonicalTablesAndTerminalRows(testCase)
fixture = makeRampFixture();
study = rampSpeed.makeStudy("ramp-speed-export-test");
study.cases = fixture.cases;
study.runs = [fixture.lateralRun fixture.longitudinalRun fixture.lateralRun];

study.runs(1).status = "completed";
study.runs(1).perSpeed.valid(1) = false;
study.runs(1).perSpeed.status(1) = "failed";
study.runs(1).perSpeed.reason(1) = "wheel lift";
study.runs(2).status = "failed";
study.runs(2).perSpeed.status(:) = "failed";
study.runs(2).perSpeed.reason(:) = "solver failure";
study.runs(3).caseId = "cancelled";
study.runs(3).status = "cancelled";
study.runs(3).runMeta.errors = "operator cancelled";
study.runs(3).perSpeed.status(:) = "cancelled";
study.runs(3).perSpeed.reason(:) = "operator cancelled";

[outputFolder,cleanup] = temporaryFolder(); %#ok<ASGLU>
result = rampSpeed.exportStudy(study,outputFolder,struct( ...
    "baseName","ramp_speed_fixture"));

verifyTrue(testCase,isfield(result,"perSpeedCsv"));
verifyTrue(testCase,isfield(result,"pointsCsv"));
verifyTrue(testCase,isfield(result,"studyMat"));
verifyTrue(testCase,isfield(result,"metadataJson"));
verifyTrue(testCase,isfield(result,"terminalLogCsv"));
verifyTrue(testCase,isfield(result,"terminalLog"));
verifyTrue(testCase,isAbsolutePath(result.perSpeedCsv));
verifyTrue(testCase,isAbsolutePath(result.pointsCsv));
verifyTrue(testCase,isAbsolutePath(result.studyMat));
verifyTrue(testCase,isAbsolutePath(result.metadataJson));
verifyTrue(testCase,isAbsolutePath(result.terminalLogCsv));
verifyTrue(testCase,isfile(result.perSpeedCsv));
verifyTrue(testCase,isfile(result.pointsCsv));
verifyTrue(testCase,isfile(result.studyMat));
verifyTrue(testCase,isfile(result.metadataJson));
verifyTrue(testCase,isfile(result.terminalLogCsv));

perSpeed = readtable(result.perSpeedCsv,"TextType","string");
variableNames = string(perSpeed.Properties.VariableNames);
verifyTrue(testCase,all(ismember(["case_id","run_status","speed_mps", ...
    "valid","status","reason"],variableNames)));
verifyEqual(testCase,perSpeed.speed_mps(1),5,"AbsTol",0);
verifyFalse(testCase,logical(perSpeed.valid(1)));
verifyEqual(testCase,perSpeed.status(1),"failed");
verifyEqual(testCase,perSpeed.reason(1),"wheel lift");
verifyEqual(testCase,perSpeed.run_status(1),"completed");

points = readtable(result.pointsCsv,"TextType","string");
pointNames = string(points.Properties.VariableNames);
verifyTrue(testCase,all(ismember(["case_id","run_status","point_reason", ...
    "speed_mps","valid","status"],pointNames)));
verifyTrue(testCase,any(points.run_status == "failed"));

verifyEqual(testCase,string(result.terminalLog.status), ...
    ["completed";"failed";"cancelled"]);
verifyEqual(testCase,height(result.terminalLog),3);
verifyTrue(testCase,isfile(result.terminalLogCsv));

metadataText = string(fileread(result.metadataJson));
verifyTrue(testCase,contains(metadataText,"m/s"));
verifyTrue(testCase,contains(metadataText,"N"));

loaded = rampSpeed.loadStudy(result.studyMat,"ignored-app-version");
[ok,issues] = rampSpeed.validateStudy(loaded);
verifyTrue(testCase,ok,strjoin(issues,newline));
verifyEqual(testCase,loaded.schemaVersion,1);
verifyEqual(testCase,loaded.runs(1).perSpeed.speed_mps, ...
    study.runs(1).perSpeed.speed_mps,"AbsTol",0);
verifyEqual(testCase,loaded.runs(1).perSpeed.reason, ...
    study.runs(1).perSpeed.reason);
end

function testRejectsInvalidCurrentSchema(testCase)
study = rampSpeed.makeStudy("ramp-speed-export-test");
invalidStudy = rmfield(study,"runs");
[outputFolder,cleanup] = temporaryFolder(); %#ok<ASGLU>

verifyError(testCase,@()rampSpeed.exportStudy(invalidStudy,outputFolder), ...
    "rampSpeed:invalidStudy");
end

function testBaseNameWithExtensionRemainsTextual(testCase)
study = rampSpeed.makeStudy("ramp-speed-export-test");
[outputFolder,cleanup] = temporaryFolder(); %#ok<ASGLU>
result = rampSpeed.exportStudy(study,outputFolder,struct( ...
    "baseName","release.csv"));

verifyTrue(testCase,endsWith(result.perSpeedCsv, ...
    "release_per_speed.csv"));
verifyTrue(testCase,isfile(result.perSpeedCsv));
end

function testMetadataLabelsCanonicalSIWhenDisplayUnitsDiffer(testCase)
fixture = makeRampFixture();
study = rampSpeed.makeStudy("ramp-speed-export-test");
study.cases = fixture.cases(1);
study.runs = fixture.lateralRun;
study.displayUnits.speed = "mph";
study.displayUnits.force = "lbf";

[outputFolder,cleanup] = temporaryFolder(); %#ok<ASGLU>
result = rampSpeed.exportStudy(study,outputFolder,struct( ...
    "baseName","canonical_units"));

verifyEqual(testCase,string(result.metadata.canonicalUnits.speed),"m/s");
verifyEqual(testCase,string(result.metadata.canonicalUnits.force),"N");
verifyEqual(testCase,result.perSpeed.speed_mps(1),5,"AbsTol",0);
verifyTrue(testCase,contains(string(fileread(result.metadataJson)), ...
    "canonicalUnits"));
end
function testExportsRequestedVisibleFigure(testCase)
assumeTrue(testCase,~isempty(which("exportgraphics")));
fixture = makeRampFixture();
study = rampSpeed.makeStudy("ramp-speed-export-test");
study.cases = fixture.cases(1);
study.runs = fixture.lateralRun;
study.runs.status = "completed";

[outputFolder,cleanupFolder] = temporaryFolder(); %#ok<ASGLU>
fig = figure("Visible","off");
cleanupFigure = onCleanup(@()close(fig)); %#ok<NASGU>
axesHandle = axes(fig); %#ok<LAXES>
plot(axesHandle,study.runs.perSpeed.speed_mps, ...
    study.runs.perSpeed.aLat_sustainable_mps2,"-o");

figureRequest = struct("figure",fig,"name","capability");
options = struct("baseName","ramp_speed_figure", ...
    "visibleFigures",figureRequest);
result = rampSpeed.exportStudy(study,outputFolder,options);

verifyTrue(testCase,isfield(result,"figureFiles"));
verifyEqual(testCase,numel(result.figureFiles),1);
verifyTrue(testCase,isAbsolutePath(result.figureFiles(1)));
verifyTrue(testCase,isfile(result.figureFiles(1)));
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

function tf = isAbsolutePath(fileName)
file = java.io.File(char(fileName));
tf = file.isAbsolute();
end
