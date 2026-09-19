function study = makeStudy(appVersion)
%MAKESTUDY Create an empty, versioned ramp-speed study container.

if nargin < 1 || isempty(appVersion)
    appVersion = "dev";
end
study = struct();
study.schemaVersion = 1;
study.created = datetime('now');
study.appVersion = string(appVersion);
study.cases = emptyCases();
study.runs = emptyRuns();
study.displayUnits = struct( ...
    'speed', "m/s", ...
    'acceleration', "m/s^2", ...
    'force', "N", ...
    'length', "m", ...
    'angle', "rad", ...
    'angularRate', "rad/s");
end

function cases = emptyCases()
cases = struct('id',{},'label',{},'source',{},'designRow',{},'carRole',{});
end

function runs = emptyRuns()
runs = struct('schemaVersion',{},'caseId',{},'type',{},'mode',{}, ...
    'settings',{},'perSpeed',{},'points',{},'runMeta',{},'status',{},'raw',{});
end
