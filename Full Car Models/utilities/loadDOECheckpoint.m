function [state,resumeInfo] = loadDOECheckpoint(path,resolvedStudy)
%LOADDOECHECKPOINT Load a compatible checkpoint and apply allowed settings.

path = validatePath(path);
saved = load(char(path),'checkpoint');
if ~isfield(saved,'checkpoint') || ~isstruct(saved.checkpoint) || ...
        ~isfield(saved.checkpoint,'resolvedStudy')
    error('loadDOECheckpoint:invalidCheckpoint', ...
        'Checkpoint must contain a state struct with resolvedStudy.')
end

state = saved.checkpoint;
savedStudy = state.resolvedStudy;
if ~isstruct(savedStudy) || ~isscalar(savedStudy) || ...
        ~isstruct(resolvedStudy) || ~isscalar(resolvedStudy) || ...
        ~isfield(savedStudy,'signature') || ~isfield(resolvedStudy,'signature') || ...
        ~isequal(string(savedStudy.signature),string(resolvedStudy.signature))
    error('loadDOECheckpoint:designMismatch', ...
        'Checkpoint design signature does not match the requested study.')
end

allowed = ["mode","maxCases","numWorkers","batchSize","objective","ramps"];
changed = strings(0,1);
for fieldName = allowed
    name = char(fieldName);
    savedHasField = isfield(savedStudy,name);
    requestedHasField = isfield(resolvedStudy,name);
    if savedHasField ~= requestedHasField || ...
            (savedHasField && ~isequaln(savedStudy.(name),resolvedStudy.(name)))
        changed(end+1,1) = fieldName; %#ok<AGROW>
    end
    if requestedHasField
        state.resolvedStudy.(name) = resolvedStudy.(name);
    end
end

resumeInfo = struct('changed',changed);
end

function path = validatePath(path)
if ~(ischar(path) || (isstring(path) && isscalar(path))) || strlength(string(path)) == 0
    error('loadDOECheckpoint:badPath','path must be a nonempty character vector or string scalar.')
end
path = string(path);
end
