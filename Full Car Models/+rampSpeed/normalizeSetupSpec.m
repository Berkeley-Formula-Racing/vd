function [spec,metadata] = normalizeSetupSpec(config,setup)
%NORMALIZESETUPSPEC Normalize editable setup data without derived state.
validateattributes(config,{'struct'},{'scalar'},mfilename,'config');
validateattributes(setup,{'struct'},{'scalar'},mfilename,'setup');
spec = setup;
if isfield(spec,'derived')
    spec = rmfield(spec,'derived');
end
if ~isfield(spec,'schemaVersion') || isempty(spec.schemaVersion)
    spec.schemaVersion = config.schemaVersion;
end
if double(spec.schemaVersion) ~= double(config.schemaVersion)
    error('rampSpeed:setupSchemaMismatch', ...
        'Setup schema version %s does not match configuration schema %s.', ...
        string(spec.schemaVersion),string(config.schemaVersion));
end
if ~isfield(spec,'baselineVersion')
    spec.baselineVersion = config.baselineVersion;
end
if string(spec.baselineVersion) ~= string(config.baselineVersion)
    error('rampSpeed:baselineVersionMismatch', ...
        'Setup baseline version %s does not match configuration %s.', ...
        string(spec.baselineVersion),string(config.baselineVersion));
end
if ~isfield(spec,'source') || isempty(spec.source)
    spec.source = "user";
else
    spec.source = string(spec.source);
end
if ~isfield(spec,'id') || ~isfield(spec,'label')
    error('rampSpeed:invalidSetup','Setup must contain id and label.');
end
spec.id = string(spec.id);
spec.label = string(spec.label);
spec.aeroMapId = string(spec.aeroMapId);
spec.isBaseline = spec.id == string(config.id);
metadata = struct('schemaVersion',double(spec.schemaVersion), ...
    'baselineVersion',string(spec.baselineVersion), ...
    'isBaseline',spec.isBaseline);
end
