function study = migrateLegacyStudy(fileName,appVersion)
%MIGRATELEGACYSTUDY Convert a legacy ramp_speed_study cache to schema v1.

if nargin < 2 || isempty(appVersion)
    appVersion = "dev";
end
fileName = requireFileName(fileName);
payload = load(fileName,"study");
if ~isfield(payload,"study") || ~isstruct(payload.study) || ...
        ~isscalar(payload.study)
    error("rampSpeed:invalidLegacyStudy", ...
        "The cache must contain a scalar study struct.");
end
legacy = payload.study;
if ~isfield(legacy,"results")
    error("rampSpeed:invalidLegacyStudy", ...
        "The legacy study does not contain study.results.");
end

results = normalizeResults(legacy.results);
nCases = numel(results);
labels = normalizeLabels(legacy,nCases);
rampOptions = legacyStruct(legacy,"rampOptions");
numWorkers = legacyValue(legacy,"numWorkers",[]);
created = legacyValue(legacy,"created",datetime.empty);
sourcePath = absolutePath(fileName);

study = rampSpeed.makeStudy(appVersion);
if isfield(legacy,"created")
    study.created = legacy.created;
end

for i = 1:nCases
    caseId = "legacy-" + string(i);
    caseInfo = struct( ...
        "id",caseId, ...
        "label",labels(i), ...
        "source",sourcePath, ...
        "designRow",i, ...
        "carRole","lap");
    legacyMeta = struct( ...
        "fileName",sourcePath, ...
        "rampOptions",rampOptions, ...
        "numWorkers",numWorkers, ...
        "created",created);
    runMeta = struct( ...
        "source","rampSpeed.migrateLegacyStudy", ...
        "legacy",legacyMeta);
    run = rampSpeed.normalizeRampResult(results{i},"lateral", ...
        rampOptions,caseInfo,runMeta);
    run.runMeta.fileName = sourcePath;
    study.cases(i) = caseInfo;
    study.runs(i) = run;
end
end

function results = normalizeResults(value)
if iscell(value)
    results = value(:);
elseif isstruct(value)
    results = num2cell(value(:));
else
    error("rampSpeed:invalidLegacyStudy", ...
        "study.results must be a cell or struct array.");
end
end

function labels = normalizeLabels(legacy,n)
if n == 0
    labels = strings(0,1);
    return
end
if isfield(legacy,"labels")
    value = legacy.labels;
    if iscell(value)
        labels = string(value(:));
    elseif ischar(value)
        labels = string(cellstr(value));
    else
        labels = string(value(:));
    end
else
    labels = strings(0,1);
end
if isempty(labels)
    labels = "case " + string((1:n).');
elseif isscalar(labels) && n > 1
    labels = repmat(labels,n,1);
elseif numel(labels) < n
    labels(end+1:n,1) = "case " + string((numel(labels)+1:n).');
else
    labels = labels(1:n);
end
end

function value = legacyStruct(legacy,name)
if isfield(legacy,name) && isstruct(legacy.(name)) && ...
        isscalar(legacy.(name))
    value = legacy.(name);
else
    value = struct();
end
end

function value = legacyValue(legacy,name,default)
if isfield(legacy,name)
    value = legacy.(name);
else
    value = default;
end
end

function fileName = requireFileName(fileName)
fileName = string(fileName);
if ~isscalar(fileName) || strlength(strtrim(fileName)) == 0
    error("rampSpeed:invalidFileName", ...
        "fileName must be a non-empty scalar path.");
end
fileName = char(fileName);
if ~isfile(fileName)
    error("rampSpeed:fileNotFound","Cache file not found: %s",fileName);
end
end

function path = absolutePath(fileName)
path = string(char(java.io.File(fileName).getAbsolutePath()));
end
