function map = validateAeroMapFile(csvPath)
%VALIDATEAEROMAPFILE Validate one CSV through the project's AeroMap contract.

csvPath = string(csvPath);
if ~isscalar(csvPath) || strlength(strtrim(csvPath)) == 0
    error('rampSpeed:invalidAeroMap', ...
        'Aero-map path must be a non-empty scalar path.');
end
try
    map = AeroMap(char(csvPath));
catch ME
    error('rampSpeed:invalidAeroMap', ...
        'Aero map %s is invalid: %s',csvPath,ME.message);
end
end
