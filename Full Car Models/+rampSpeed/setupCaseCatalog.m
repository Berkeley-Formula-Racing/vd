function cases = setupCaseCatalog(carCell,designTable,carRole)
%SETUPCASECATALOG Build stable, display-ready setup records.
%   The car configuration convention is one row per design point and two
%   columns per row: lap car first, acceleration car second.

if nargin < 1 || isempty(carCell) || ~iscell(carCell)
    error("rampSpeed:invalidCarCell", ...
        "carCell must be a non-empty cell array of vehicle setups.");
end
if nargin < 2 || isempty(designTable)
    designTable = table();
elseif ~istable(designTable)
    error("rampSpeed:invalidDesignTable", ...
        "designTable must be a table or empty.");
end
if nargin < 3 || isempty(carRole)
    carRole = "auto";
end
role = normalizeRole(carRole);

if ~isempty(designTable)
    nCases = height(designTable);
    if nCases > size(carCell,1) && ~(isvector(carCell) && ...
            numel(carCell) >= nCases)
        error("rampSpeed:designRowCount", ...
            "designTable must have one row per available car setup.");
    end
else
    nCases = size(carCell,1);
    % A 1-by-2 cell array is the conventional single setup with its two
    % role-specific cars, not two independent design points.
    if nCases == 1 && size(carCell,2) > 2
        nCases = numel(carCell);
    end
end

if nCases == 0
    cases = repmat(emptyCase(),0,1);
    return
end

cases = repmat(emptyCase(),nCases,1);
usedIds = strings(0,1);
for i = 1:nCases
    baseId = tableText(designTable,i, ...
        ["id","caseId","case_id","setupId"]);
    if strlength(strtrim(baseId)) == 0
        baseId = "car-" + compose("%03d",i);
    end
    id = uniqueId(baseId,usedIds);
    usedIds(end+1,1) = id; %#ok<AGROW>

    label = tableText(designTable,i, ...
        ["label","name","setup","caseLabel","caseName"]);
    if strlength(strtrim(label)) == 0
        rowName = designRowName(designTable,i);
        if strlength(rowName) > 0
            label = rowName;
        else
            label = "setup " + string(i);
        end
    end

    source = tableText(designTable,i,["source","origin"]);
    if strlength(strtrim(source)) == 0
        source = "carConfig";
    end

    cases(i).id = id;
    cases(i).label = strtrim(label);
    cases(i).source = strtrim(source);
    cases(i).designRow = i;
    cases(i).sourceIndex = i;
    cases(i).carRole = role;
    cases(i).carColumn = roleColumn(role);
end
end

function value = emptyCase()
value = struct('id',"",'label',"",'source',"", ...
    'designRow',NaN,'sourceIndex',NaN,'carRole',"auto",'carColumn',0);
end

function value = normalizeRole(value)
value = lower(strtrim(string(value)));
if ~isscalar(value)
    error("rampSpeed:invalidCarRole", ...
        "carRole must be a scalar text value.");
end
switch value
    case {"","auto","default"}
        value = "auto";
    case {"lap","lateral"}
        value = "lap";
    case {"accel","acceleration","longitudinal","acceleration-car"}
        value = "acceleration";
    otherwise
        error("rampSpeed:invalidCarRole", ...
            "carRole must be auto, lap, or acceleration.");
end
end

function value = roleColumn(role)
switch role
    case "lap"
        value = 1;
    case "acceleration"
        value = 2;
    otherwise
        value = 0;
end
end

function value = tableText(T,index,candidates)
value = "";
if isempty(T)
    return
end
names = string(T.Properties.VariableNames);
for candidate = string(candidates(:).')
    match = find(strcmpi(names,candidate),1);
    if isempty(match)
        continue
    end
    raw = T{index,match};
    if iscell(raw)
        raw = raw{1};
    end
    candidateValue = string(raw);
    if isscalar(candidateValue) && strlength(strtrim(candidateValue)) > 0
        value = candidateValue;
        return
    end
end
end

function value = designRowName(T,index)
value = "";
if isempty(T) || isempty(T.Properties.RowNames)
    return
end
value = string(T.Properties.RowNames{index});
if ~isscalar(value)
    value = "";
end
end

function value = uniqueId(base,used)
value = strtrim(string(base));
if strlength(value) == 0
    value = "case";
end
if ~any(used == value)
    return
end
suffix = 2;
candidate = value + "-" + string(suffix);
while any(used == candidate)
    suffix = suffix + 1;
    candidate = value + "-" + string(suffix);
end
value = candidate;
end
