function tests = test_rampSpeedFixedSpeedPresets
tests = functiontests(localfunctions);
end

function testCatalogContainsThreeDeterministicFixedPresets(testCase)
presets = rampSpeed.fixedSpeedPresets();

verifyEqual(testCase,string({presets.id}), ...
    ["preview","accurate","highAccuracy"]);
verifyEqual(testCase,[presets.pointCount],[5 11 21]);
verifyTrue(testCase,all([presets.fixed]));
verifyFalse(testCase,any(contains(lower(string({presets.label})),"adaptive")));
end

function testFixedSpeedGridUsesPresetPointCountAndRange(testCase)
[speeds,meta] = rampSpeed.fixedSpeedGrid([5 25],"highAccuracy");

verifyEqual(testCase,numel(speeds),21);
verifyEqual(testCase,speeds(1),5,"AbsTol",0);
verifyEqual(testCase,speeds(end),25,"AbsTol",0);
verifyEqual(testCase,all(diff(speeds) > 0),true);
verifyEqual(testCase,string(meta.mode),"fixed");
verifyEqual(testCase,string(meta.preset),"highAccuracy");
verifyEqual(testCase,meta.pointCount,21);
verifyEqual(testCase,meta.range_mps,[5;25],"AbsTol",0);
end

function testFixedSpeedGridRejectsAdaptiveNames(testCase)
verifyError(testCase,@()rampSpeed.fixedSpeedGrid([5 25],"adaptive"), ...
    "rampSpeed:invalidFixedSpeedPreset");
end

function testSharedResolverAppliesPresetToEitherRampEntryPoint(testCase)
requested = 5:2.5:30;
[speeds,meta] = rampSpeed.resolveFixedSpeedGrid(requested, ...
    struct("mode","highAccuracy"));

verifyEqual(testCase,numel(speeds),21);
verifyEqual(testCase,string(meta.preset),"highAccuracy");
verifyEqual(testCase,speeds([1 end]),[5;30],"AbsTol",0);
end

function testLongitudinalEntryPointEvaluatesEveryHighAccuracyPoint(testCase)
ensureRampSpeedAppPath();
[car,~] = carConfigBaseline();
settings = struct("speeds",[5 25],"speedGrid", ...
    struct("mode","highAccuracy"),"solverProfile","fastPreview", ...
    "verbose",false);
run = rampSpeed.runLongitudinalRamp(car,settings, ...
    struct("id","fixed-preset-test","label","fixed preset test"),struct());

verifyEqual(testCase,height(run.perSpeed),21);
verifyEqual(testCase,string(run.runMeta.speedGrid.preset),"highAccuracy");
verifyEqual(testCase,run.runMeta.speedGrid.pointCount,21);
verifyEqual(testCase,run.perSpeed.speed_mps(1),5,"AbsTol",0);
verifyEqual(testCase,run.perSpeed.speed_mps(end),25,"AbsTol",0);
end

function ensureRampSpeedAppPath()
root = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(root));
end
