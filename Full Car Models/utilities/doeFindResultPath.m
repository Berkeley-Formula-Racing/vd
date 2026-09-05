function resultPath = doeFindResultPath(modelRoot,configuredOutput,requestedSource)
%DOEFINDRESULTPATH Select a DOE result/checkpoint without stale legacy picks.
%   Prefers an explicitly requested source, then the configured output
%   directory, then the newest DOE_output* directory. Root-level legacy
%   files are retained only as the final compatibility fallback.

if nargin < 3, requestedSource = ""; end
modelRoot = char(string(modelRoot));
configuredOutput = char(string(configuredOutput));
requestedSource = string(requestedSource);
if strlength(requestedSource) > 0
    resultPath = resolveSource(char(requestedSource));
    return
end

if isfolder(configuredOutput)
    path = sourceInDirectory(configuredOutput);
    if strlength(path) > 0
        resultPath = char(path);
        return
    end
end

listing = dir(fullfile(modelRoot,'DOE_output*'));
candidates = strings(0,1);
timestamps = zeros(0,1);
for i = 1:numel(listing)
    if ~listing(i).isdir, continue, end
    path = sourceInDirectory(fullfile(listing(i).folder,listing(i).name));
    if strlength(path) == 0, continue, end
    info = dir(char(path));
    candidates(end+1,1) = path; %#ok<AGROW>
    timestamps(end+1,1) = info.datenum; %#ok<AGROW>
end
if ~isempty(candidates)
    [~,pick] = max(timestamps);
    resultPath = char(candidates(pick));
    return
end

legacy = [fullfile(modelRoot,'DOE_results.mat'); ...
    fullfile(modelRoot,'DOE_checkpoint.mat'); ...
    fullfile(fileparts(modelRoot),'DOE_results.mat'); ...
    fullfile(fileparts(modelRoot),'DOE_checkpoint.mat')];
for i = 1:numel(legacy)
    if isfile(legacy{i})
        resultPath = legacy{i};
        return
    end
end
error('doeFindResultPath:noResults', ...
    'No DOE results or checkpoint file is available below %s.',modelRoot)
end

function path = resolveSource(source)
if isfolder(source)
    path = sourceInDirectory(source);
    if strlength(path) > 0, path = char(path); return, end
    error('doeFindResultPath:noResults','No DOE result or checkpoint exists in %s.',source)
end
if isfile(source)
    path = source;
    return
end
error('doeFindResultPath:missingSource','DOE source does not exist: %s.',source)
end

function path = sourceInDirectory(directory)
result = fullfile(directory,'DOE_results.mat');
checkpoint = fullfile(directory,'DOE_checkpoint.mat');
if isfile(result)
    path = string(result);
elseif isfile(checkpoint)
    path = string(checkpoint);
else
    path = "";
end
end
