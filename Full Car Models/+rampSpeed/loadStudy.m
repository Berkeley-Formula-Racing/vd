function study = loadStudy(fileName,appVersion)
%LOADSTUDY Load a schema-v1 study or migrate a legacy cache.

if nargin < 2 || isempty(appVersion)
    appVersion = "dev";
end
fileName = requireFileName(fileName);
payload = load(fileName,"study");
if ~isfield(payload,"study") || ~isstruct(payload.study) || ...
        ~isscalar(payload.study)
    error("rampSpeed:invalidStudy", ...
        "The cache must contain a scalar study struct.");
end
candidate = payload.study;

if isfield(candidate,"schemaVersion")
    version = double(candidate.schemaVersion);
    if ~isscalar(version) || ~isfinite(version)
        error("rampSpeed:invalidStudy", ...
            "study.schemaVersion must be a finite scalar.");
    end
    if version > 1
        error("rampSpeed:unsupportedSchema", ...
            "Unsupported ramp-speed study schema version %g.",version);
    end
    if version == 1
        [ok,issues] = rampSpeed.validateStudy(candidate);
        if ~ok
            error("rampSpeed:invalidStudy", ...
                "Invalid ramp-speed study: %s",joinIssues(issues));
        end
        if ~isfield(candidate,'setupSpecifications') || isempty(candidate.setupSpecifications)
            candidate.readOnly = true;
        elseif ~isfield(candidate,'readOnly')
            candidate.readOnly = false;
        end
        study = candidate;
        return
    end
end

if isfield(candidate,"results")
    study = rampSpeed.migrateLegacyStudy(fileName,appVersion);
    study.readOnly = true;
else
    error("rampSpeed:invalidStudy", ...
        "The cache is neither a schema-v1 study nor a legacy study.");
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

function message = joinIssues(issues)
issues = string(issues(:));
message = char(strjoin(issues,"; "));
end
