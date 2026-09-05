function tests = test_aeroMapEnabled
tests = functiontests(localfunctions);
end

function testUseAeroMapOverridesLegacyMapEnabled(testCase)
% A new explicit switch must win so a user can disable the map without
% deleting any legacy configuration field.
config = struct('use_aeromap',false,'map_enabled',true);

verifyFalse(testCase,aeroMapEnabled(config));
end

function testLegacyMapEnabledRemainsSupported(testCase)
verifyTrue(testCase,aeroMapEnabled(struct('map_enabled',true)));
verifyFalse(testCase,aeroMapEnabled(struct('map_enabled',false)));
end

function testMissingToggleUsesStaticAero(testCase)
verifyFalse(testCase,aeroMapEnabled(struct()));
end
