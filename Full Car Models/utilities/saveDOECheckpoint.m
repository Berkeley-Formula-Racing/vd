function saveDOECheckpoint(path,state)
%SAVEDOECHECKPOINT Atomically persist adaptive DOE state as a v7.3 MAT-file.

path = validatePath(path);
tmp = path + ".tmp.mat";
cleanup = onCleanup(@() deleteIfPresent(tmp)); %#ok<NASGU>

checkpoint = state;
save(char(tmp),'checkpoint','-v7.3')

temporary = load(char(tmp),'checkpoint');
if ~isfield(temporary,'checkpoint')
    error('saveDOECheckpoint:tempValidation', ...
        'Temporary checkpoint did not contain the checkpoint variable.')
end

[moved,message] = movefile(char(tmp),char(path),'f');
if ~moved
    error('saveDOECheckpoint:replaceFailed', ...
        'Could not replace checkpoint %s: %s',path,message)
end
end

function path = validatePath(path)
if ~(ischar(path) || (isstring(path) && isscalar(path))) || strlength(string(path)) == 0
    error('saveDOECheckpoint:badPath','path must be a nonempty character vector or string scalar.')
end
path = string(path);
end

function deleteIfPresent(path)
if isfile(path)
    delete(path)
end
end
