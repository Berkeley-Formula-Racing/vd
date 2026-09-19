function tests = test_rampStudyRoundTrip
tests = functiontests(localfunctions);
end

function setupOnce(~)
setup_paths;
end

function testCurrentSchemaRoundTripsWithFinalPathMetadata(testCase)
fixture = makeRampFixture();
study = rampSpeed.makeStudy("roundtrip-app");
study.created = datetime(2026,1,2,0,0,0);
study.cases = fixture.cases(1);
study.runs = fixture.lateralRun;

[target,folder] = temporaryTarget("current.mat");
cleanup = onCleanup(@()deleteIfPresent(target,folder));
rampSpeed.saveStudy(target,study);
loaded = rampSpeed.loadStudy(target,"ignored-app");

verifyEqual(testCase,loaded.schemaVersion,study.schemaVersion);
verifyEqual(testCase,loaded.appVersion,study.appVersion);
verifyEqual(testCase,string(loaded.cases.id),string(study.cases.id));
verifyEqual(testCase,string(loaded.runs.caseId),string(study.runs.caseId));
verifyEqual(testCase,loaded.runs.status,study.runs.status);
verifyEqual(testCase,loaded.runs.perSpeed.speed_mps, ...
    study.runs.perSpeed.speed_mps,"AbsTol",0);
verifyEqual(testCase,loaded.runs.perSpeed.aLat_sustainable_mps2, ...
    study.runs.perSpeed.aLat_sustainable_mps2,"AbsTol",0);
verifyEqual(testCase,loaded.runs.raw,study.runs.raw);
verifyEqual(testCase,loaded.runs.runMeta.fileName,absolutePath(target));
end

function testLoadStudyRejectsUnsupportedSchema(testCase)
study = rampSpeed.makeStudy("unsupported");
study.schemaVersion = 2;
[target,folder] = temporaryTarget("unsupported.mat");
cleanup = onCleanup(@()deleteIfPresent(target,folder));
save(target,"study","-v7.3");

verifyError(testCase,@()rampSpeed.loadStudy(target,"test-app"), ...
    "rampSpeed:unsupportedSchema");
end

function testSaveValidationDoesNotReplaceExistingTarget(testCase)
fixture = makeRampFixture();
study = rampSpeed.makeStudy("atomic");
study.cases = fixture.cases(1);
study.runs = fixture.lateralRun;
[target,folder] = temporaryTarget("atomic.mat");
cleanup = onCleanup(@()deleteIfPresent(target,folder));
rampSpeed.saveStudy(target,study);
bytesBefore = readBytes(target);

invalid = rmfield(study,"displayUnits");
verifyError(testCase,@()rampSpeed.saveStudy(target,invalid), ...
    "rampSpeed:invalidStudy");

verifyEqual(testCase,readBytes(target),bytesBefore);
end

function [target,folder] = temporaryTarget(name)
folder = tempname;
mkdir(folder);
target = fullfile(folder,name);
end

function path = absolutePath(fileName)
path = string(char(java.io.File(fileName).getAbsolutePath()));
end

function bytes = readBytes(fileName)
fid = fopen(fileName,"r");
cleanup = onCleanup(@()fclose(fid));
bytes = fread(fid,Inf,"*uint8");
end

function deleteIfPresent(fileName,folder)
if isfile(fileName)
    delete(fileName);
end
if isfolder(folder)
    rmdir(folder,"s");
end
end
