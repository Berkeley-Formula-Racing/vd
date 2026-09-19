function tests = test_rampStudyMigration
tests = functiontests(localfunctions);
end

function setupOnce(~)
setup_paths;
end

function testMigratesLegacyStudyWithoutMutatingSource(testCase)
fixture = makeRampFixture();
[legacyFile,sourceStudy] = writeLegacyStudy(fixture.legacyRampResult, ...
    "baseline",[5 10]);
cleanup = onCleanup(@()deleteIfPresent(legacyFile));

before = dir(legacyFile);
bytesBefore = readBytes(legacyFile);
study = rampSpeed.migrateLegacyStudy(legacyFile,"test-app");
after = dir(legacyFile);

verifyEqual(testCase,after.bytes,before.bytes);
verifyEqual(testCase,after.datenum,before.datenum);
verifyEqual(testCase,readBytes(legacyFile),bytesBefore);
verifyEqual(testCase,study.schemaVersion,1);
verifyEqual(testCase,study.appVersion,"test-app");
verifyNumElements(testCase,study.cases,1);
verifyNumElements(testCase,study.runs,1);
verifyEqual(testCase,string(study.cases(1).label),"baseline");
verifyNotEmpty(testCase,string(study.cases(1).id));
verifyEqual(testCase,string(study.runs(1).caseId), ...
    string(study.cases(1).id));
verifyEqual(testCase,study.runs(1).type,"lateral");
verifyEqual(testCase,study.runs(1).mode,"coast");
verifyEqual(testCase,study.runs(1).settings.speeds,[5 10]);
verifyEqual(testCase,study.runs(1).settings.nRamp,4);
verifyEqual(testCase,study.runs(1).settings.mode,"coast");
verifyEqual(testCase,study.runs(1).runMeta.legacy.rampOptions, ...
    sourceStudy.rampOptions);
verifyEqual(testCase,study.runs(1).runMeta.legacy.numWorkers,0);
verifyEqual(testCase,study.runs(1).runMeta.legacy.created, ...
    sourceStudy.created);
verifyEqual(testCase,study.runs(1).raw,sourceStudy.results{1});
verifyEqual(testCase,study.runs(1).perSpeed.aLat_sustainable_mps2, ...
    sourceStudy.results{1}.perSpeed.gLat_top*9.80665, ...
    "AbsTol",1e-12);
end

function testMigrationPreservesMissingSpeedAsInvalidGap(testCase)
fixture = makeRampFixture();
raw = fixture.legacyRampResult;
raw.perSpeed = raw.perSpeed(1,:);
raw.points = raw.points(1,:);
[legacyFile,~] = writeLegacyStudy(raw,"baseline",[5 10]);
cleanup = onCleanup(@()deleteIfPresent(legacyFile));

study = rampSpeed.migrateLegacyStudy(legacyFile,"test-app");
perSpeed = study.runs(1).perSpeed;

verifyEqual(testCase,perSpeed.speed_mps,[5;10],"AbsTol",1e-12);
verifyTrue(testCase,perSpeed.valid(1));
verifyFalse(testCase,perSpeed.valid(2));
verifyEqual(testCase,perSpeed.status(2),"missing");
verifyTrue(testCase,isnan(perSpeed.aLat_sustainable_mps2(2)));
verifyThat(testCase,perSpeed.reason(2), ...
    matlab.unittest.constraints.ContainsSubstring("requested speed"));
end

function testLoadStudyRoutesLegacyCacheThroughMigration(testCase)
fixture = makeRampFixture();
[legacyFile,~] = writeLegacyStudy(fixture.legacyRampResult, ...
    "baseline",[5 10]);
cleanup = onCleanup(@()deleteIfPresent(legacyFile));

study = rampSpeed.loadStudy(legacyFile,"loaded-app");

verifyEqual(testCase,study.schemaVersion,1);
verifyEqual(testCase,study.appVersion,"loaded-app");
verifyEqual(testCase,string(study.cases(1).label),"baseline");
verifyEqual(testCase,height(study.runs(1).perSpeed),2);
verifyEqual(testCase,study.runs(1).raw,fixture.legacyRampResult);
end

function [fileName,study] = writeLegacyStudy(raw,label,speeds)
folder = tempname;
mkdir(folder);
fileName = fullfile(folder,"legacy.mat");
created = datetime(2026,1,1,0,0,0);
study = struct();
study.results = {raw};
study.labels = string(label);
study.rampOptions = struct("speeds",speeds,"nRamp",4,"mode","coast");
study.numWorkers = 0;
study.created = created;
save(fileName,"study","-v7.3");
end

function bytes = readBytes(fileName)
fid = fopen(fileName,"r");
cleanup = onCleanup(@()fclose(fid));
bytes = fread(fid,Inf,"*uint8");
end

function deleteIfPresent(fileName)
if isfile(fileName)
    delete(fileName);
end
folder = fileparts(fileName);
if isfolder(folder)
    rmdir(folder,"s");
end
end
