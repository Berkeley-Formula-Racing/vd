function [ok,issues] = validateAeroMapCatalog(catalog)
%VALIDATEAEROMAPCATALOG Validate every entry in an aero-map catalog.

ok = true;
issues = strings(0,1);
if ~isstruct(catalog)
    ok = false;
    issues(end+1) = "catalog must be a struct array";
    return
end
for i = 1:numel(catalog)
    if ~isfield(catalog(i),'id') || strlength(string(catalog(i).id)) == 0
        ok = false;
        issues(end+1) = sprintf('entry %d is missing an ID',i);
    end
    if ~isfield(catalog(i),'path') || strlength(string(catalog(i).path)) == 0
        ok = false;
        issues(end+1) = sprintf('entry %d is missing a path',i);
        continue
    end
    try
        rampSpeed.validateAeroMapFile(catalog(i).path);
    catch ME
        ok = false;
        issues(end+1) = string(ME.message);
    end
end
issues = unique(issues,'stable');
end
