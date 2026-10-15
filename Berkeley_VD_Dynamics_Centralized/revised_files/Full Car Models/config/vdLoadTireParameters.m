function parameters = vdLoadTireParameters(path,expectedCount)
%VDLOADTIREPARAMETERS Load one fit without importing variables into a workspace.
if ~isfile(path)
    error('vdTireConfig:missingFile','Tire parameter file does not exist: %s',path);
end
data = load(path,'Xbestcell');
if ~isfield(data,'Xbestcell') || ~iscell(data.Xbestcell)
    error('vdTireConfig:badFile','Expected Xbestcell in %s.',path);
end
parameters = cell2mat(data.Xbestcell);
if ~isnumeric(parameters) || ~isreal(parameters) || ...
        ~isvector(parameters) || numel(parameters) ~= expectedCount || ...
        any(~isfinite(parameters(:)))
    error('vdTireConfig:badParameters', ...
        'Expected %d finite real coefficients in %s.',expectedCount,path);
end
end
