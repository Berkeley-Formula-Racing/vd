function saveStudy(fileName,study)
%SAVESTUDY Validate and atomically save a schema-v1 ramp-speed study.

fileName = requireFileName(fileName);
[ok,issues] = rampSpeed.validateStudy(study);
if ~ok
    error("rampSpeed:invalidStudy", ...
        "Invalid ramp-speed study: %s",joinIssues(issues));
end

target = absolutePath(fileName);
folder = fileparts(char(target));
if isempty(folder)
    folder = pwd;
end
if ~isfolder(folder)
    error("rampSpeed:invalidPath","Target folder does not exist: %s",folder);
end

for i = 1:numel(study.runs)
    study.runs(i).runMeta.fileName = target;
end

temporaryFile = [char(tempname(folder)) '.mat'];
cleanup = onCleanup(@()deleteIfPresent(temporaryFile));
save(temporaryFile,'study','-v7.3');
if ~isfile(temporaryFile)
    error("rampSpeed:saveFailed", ...
        "MATLAB did not create the temporary study file.");
end
[moved,message,messageId] = movefile(temporaryFile,char(target),"f");
if ~moved
    if isempty(messageId)
        error("rampSpeed:saveFailed", ...
            "Could not replace %s: %s",target,message);
    end
    error("rampSpeed:saveFailed", ...
        "Could not replace %s (%s): %s",target,messageId,message);
end
end

function fileName = requireFileName(fileName)
fileName = string(fileName);
if ~isscalar(fileName) || strlength(strtrim(fileName)) == 0
    error("rampSpeed:invalidFileName", ...
        "fileName must be a non-empty scalar path.");
end
fileName = char(fileName);
end

function path = absolutePath(fileName)
path = string(char(java.io.File(fileName).getAbsolutePath()));
end

function message = joinIssues(issues)
issues = string(issues(:));
message = char(strjoin(issues,"; "));
end

function deleteIfPresent(fileName)
if isfile(fileName)
    delete(fileName);
end
end
