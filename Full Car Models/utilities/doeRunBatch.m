function results = doeRunBatch(carCell,eventParams,study,caseIndices,runMode)
%DOERUNBATCH Run DOE cases serially or with one outer case-level parfor.

if nargin < 4 || isempty(caseIndices)
    caseIndices = (1:size(carCell,1)).';
end
if nargin < 5 || isempty(runMode), runMode = "full"; end

n = size(carCell,1);
if size(carCell,2) < 2
    error('doeRunBatch:badCarCell', ...
        'carCell must have lap and acceleration cars in columns 1 and 2.')
end
if ~isnumeric(caseIndices) || numel(caseIndices) ~= n
    error('doeRunBatch:badCaseIndices', ...
        'caseIndices must contain one numeric index per car case.')
end
if ~isfield(study,'numWorkers') || ~isscalar(study.numWorkers) || ...
        ~isnumeric(study.numWorkers) || study.numWorkers < 0 || ...
        study.numWorkers ~= floor(study.numWorkers)
    error('doeRunBatch:badNumWorkers', ...
        'study.numWorkers must be a nonnegative integer.')
end

caseIndices = caseIndices(:);
results = cell(n,1);
if study.numWorkers == 0
    for j = 1:n
        results{j} = doeRunCase(carCell{j,1},carCell{j,2}, ...
            eventParams,study,caseIndices(j),runMode);
    end
else
    parfor (j = 1:n,study.numWorkers)
        results{j} = doeRunCase(carCell{j,1},carCell{j,2}, ...
            eventParams,study,caseIndices(j),runMode);
    end
end
end
