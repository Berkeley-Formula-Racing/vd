function tests = test_rampSpeedAssetResolution
tests = functiontests(localfunctions);
end

function testAssetResolverUsesProjectRelativePaths(testCase)
[~,config] = carConfigBaseline();
assets = rampSpeed.loadRampAssets(config,config.defaultSetup);

verifyClass(testCase,assets.map,'AeroMap');
verifyTrue(testCase,isfile(assets.mapPath));
verifyTrue(testCase,isfile(assets.camberRatioPath));
verifyTrue(testCase,isfile(assets.camberModelPath));
verifyEqual(testCase,strlength(assets.fingerprints.aeroMap),64);
verifyEqual(testCase,strlength(assets.fingerprints.camberRatios),64);
verifyEqual(testCase,strlength(assets.fingerprints.camberModels),64);
end

function testFingerprintChangesWithFileContent(testCase)
path = [tempname '.bin'];
cleanup = onCleanup(@()deleteIfPresent(path));
fid = fopen(path,'w');
fwrite(fid,uint8([1 2 3 4]),'uint8');
fclose(fid);
first = rampSpeed.fingerprintFile(path);

fid = fopen(path,'w');
fwrite(fid,uint8([1 2 3 5]),'uint8');
fclose(fid);
second = rampSpeed.fingerprintFile(path);

verifyNotEqual(testCase,first,second);
verifyEqual(testCase,strlength(first),64);
verifyEqual(testCase,strlength(second),64);
end

function testBuildReturnsDataOnlySetupWithoutStaleDerivedField(testCase)
[~,config] = carConfigBaseline();
setup = config.defaultSetup;
setup.derived = struct('R_sf',-1,'aeroMapId',"stale");
[car,normalized,derived] = rampSpeed.buildCarFromSetup(config,setup);

verifyClass(testCase,car,'Car');
verifyFalse(testCase,isfield(normalized,'derived'));
verifyEqual(testCase,car.R_sf,derived.R_sf,'AbsTol',1e-12);
end

function testBuildFromResolvedAssetsWorksOutsideCurrentDirectory(testCase)
[~,config] = carConfigBaseline();
assets = rampSpeed.loadRampAssets(config,config.defaultSetup);
oldDirectory = pwd;
cleanup = onCleanup(@()cd(oldDirectory));
cd(tempdir);

[car,normalized,derived] = rampSpeed.buildCarFromSetup(config,config.defaultSetup);

verifyClass(testCase,car,'Car');
verifyEqual(testCase,normalized.aeroMapId,"b26");
verifyEqual(testCase,derived.aeroMapFingerprint,assets.fingerprints.aeroMap);
end

function testDuplicateSetupDoesNotRetainDerivedValues(testCase)
[baselineCar,config] = carConfigBaseline();
duplicate = rampSpeed.duplicateSetup(config.defaultSetup,"asset-copy","Asset copy");
duplicate.rearArbStiffness_NmPerRad = config.options.rearArbStiffness_NmPerRad(end);
[~,normalized,derived] = rampSpeed.buildCarFromSetup(config,duplicate);

verifyFalse(testCase,isfield(normalized,'derived'));
verifyEqual(testCase,normalized.id,"asset-copy");
verifyNotEqual(testCase,derived.R_sf,baselineCar.R_sf);
end

function deleteIfPresent(path)
if isfile(path)
    delete(path);
end
end
