function data = buildPlotData(runs,metricId,options)
%BUILDPLOTDATA Build solver-independent data for one metric and its setups.

if nargin < 3 || isempty(options)
    options = struct();
end
runList = normalizeRuns(runs);
metric = findMetric(metricId);

data = struct('metricId',metric.id,'metric',metric, ...
    'sourceLevel',metric.sourceLevel,'xLabel',"speed (m/s)", ...
    'yLabel',metric.yLabel,'series',repmat(seriesTemplate(),0,1));

for i = 1:numel(runList)
    run = runList(i);
    if ~isApplicable(metric,run)
        continue
    end
    setupId = setupIdentifier(run,options,i);
    label = setupLabel(run,options,i,setupId);
    if metric.sourceLevel == "points"
        series = buildPointSeries(run,metric,setupId,label);
    else
        series = buildPerSpeedSeries(run,metric,setupId,label);
    end
    data.series(end+1,1) = series;
end
end

function metric = findMetric(metricId)
catalog = rampSpeed.metricCatalog();
id = string(metricId);
matches = strcmp(string({catalog.id}),id);
if ~any(matches)
    error('rampSpeed:unknownMetric','Unknown ramp-speed metric "%s".',id);
end
metric = catalog(find(matches,1,'first'));
end

function runList = normalizeRuns(runs)
if isempty(runs)
    runList = struct([]);
elseif iscell(runs)
    runList = [runs{:}];
else
    runList = runs(:);
end
end

function applicable = isApplicable(metric,run)
if ~isfield(run,'type') || isempty(run.type)
    applicable = true;
    return
end
type = lower(string(run.type));
applicable = any(metric.validTypes == type) || any(metric.validTypes == "any");
end

function id = setupIdentifier(run,options,index)
if isfield(options,'setupIds') && numel(options.setupIds) >= index
    id = string(options.setupIds(index));
elseif isfield(run,'caseId') && ~isempty(run.caseId)
    id = string(run.caseId);
elseif isfield(run,'runMeta') && isstruct(run.runMeta) && ...
        isfield(run.runMeta,'caseInfo') && isstruct(run.runMeta.caseInfo) && ...
        isfield(run.runMeta.caseInfo,'id')
    id = string(run.runMeta.caseInfo.id);
else
    id = "setup" + index;
end
if strlength(id) == 0
    id = "setup" + index;
end
end

function label = setupLabel(run,options,index,setupId)
if isfield(options,'labels') && numel(options.labels) >= index
    label = string(options.labels(index));
elseif isfield(run,'runMeta') && isstruct(run.runMeta) && ...
        isfield(run.runMeta,'caseInfo') && isstruct(run.runMeta.caseInfo) && ...
        isfield(run.runMeta.caseInfo,'label')
    label = string(run.runMeta.caseInfo.label);
else
    label = setupId;
end
if strlength(label) == 0
    label = setupId;
end
end

function series = buildPerSpeedSeries(run,metric,setupId,label)
series = seriesTemplate();
series.id = setupId;
series.setupId = setupId;
series.label = label;
series.sourceLevel = metric.sourceLevel;
if ~isfield(run,'perSpeed') || ~istable(run.perSpeed)
    return
end
T = run.perSpeed;
x = numericColumn(T,"speed_mps",NaN);
[rawValues,hasValues] = metricValues(T,metric);
values = scaledValues(rawValues,metric.scale);
rowValid = logicalColumn(T,"valid",true);
ruleValid = validityMask(T,metric.validityRule);
plotValid = rowValid & ruleValid & isfinite(x) & hasValues & isfinite(values);
values(~plotValid) = NaN;

series.x = x;
series.values = values;
series.valid = rowValid;
series.truncated = logicalColumn(T,"truncated",false);
series.power_limited = logicalColumn(T,"power_limited",false);
series.wheel_lift = logicalColumn(T,"wheel_lift",false);
series.aero_outside_map = logicalColumn(T,"aero_outside_map",false);
series.reason = stringColumn(T,"reason","");
series.status = stringColumn(T,"status","");
end

function series = buildPointSeries(run,metric,setupId,label)
series = seriesTemplate();
series.id = setupId;
series.setupId = setupId;
series.label = label;
series.sourceLevel = metric.sourceLevel;
if ~isfield(run,'points') || ~istable(run.points)
    return
end
T = run.points;
x = numericColumn(T,"speed_mps",NaN);
[rawValues,hasValues] = metricValues(T,metric);
values = scaledValues(rawValues,metric.scale);
rowValid = logicalColumn(T,"valid",true);
ruleValid = validityMask(T,metric.validityRule);
plotValid = rowValid & ruleValid & isfinite(x) & hasValues & isfinite(values);
values(~plotValid) = NaN;
speedIndex = numericColumn(T,"speed_index",NaN);
finiteIndex = isfinite(speedIndex);
ids = unique(speedIndex(finiteIndex),'stable');
groups = repmat(groupTemplate(),0,1);
flatX = zeros(0,1);
flatValues = zeros(0,1);
flatValid = false(0,1);
flatTruncated = false(0,1);
flatPowerLimited = false(0,1);
flatWheelLift = false(0,1);
flatAeroOutsideMap = false(0,1);
flatReason = strings(0,1);
flatStatus = strings(0,1);
for i = 1:numel(ids)
    row = finiteIndex & speedIndex == ids(i);
    group = groupTemplate();
    group.speedIndex = ids(i);
    group.x = x(row);
    group.values = values(row);
    group.valid = rowValid(row);
    group.truncated = logicalColumn(T,"truncated",false);
    group.truncated = group.truncated(row);
    group.power_limited = logicalColumn(T,"power_limited",false);
    group.power_limited = group.power_limited(row);
    group.wheel_lift = logicalColumn(T,"wheel_lift",false);
    group.wheel_lift = group.wheel_lift(row);
    group.aero_outside_map = logicalColumn(T,"aero_outside_map",false);
    group.aero_outside_map = group.aero_outside_map(row);
    group.reason = stringColumn(T,"reason","");
    group.reason = group.reason(row);
    group.status = stringColumn(T,"status","");
    group.status = group.status(row);
    groups(end+1,1) = group; %#ok<AGROW>
    if ~isempty(flatX)
        flatX(end+1,1) = NaN; %#ok<AGROW>
        flatValues(end+1,1) = NaN; %#ok<AGROW>
        flatValid(end+1,1) = false; %#ok<AGROW>
        flatTruncated(end+1,1) = false; %#ok<AGROW>
        flatPowerLimited(end+1,1) = false; %#ok<AGROW>
        flatWheelLift(end+1,1) = false; %#ok<AGROW>
        flatAeroOutsideMap(end+1,1) = false; %#ok<AGROW>
        flatReason(end+1,1) = ""; %#ok<AGROW>
        flatStatus(end+1,1) = ""; %#ok<AGROW>
    end
    flatX = [flatX; group.x]; %#ok<AGROW>
    flatValues = [flatValues; group.values]; %#ok<AGROW>
    flatValid = [flatValid; group.valid]; %#ok<AGROW>
    flatTruncated = [flatTruncated; group.truncated]; %#ok<AGROW>
    flatPowerLimited = [flatPowerLimited; group.power_limited]; %#ok<AGROW>
    flatWheelLift = [flatWheelLift; group.wheel_lift]; %#ok<AGROW>
    flatAeroOutsideMap = [flatAeroOutsideMap; group.aero_outside_map]; %#ok<AGROW>
    flatReason = [flatReason; group.reason]; %#ok<AGROW>
    flatStatus = [flatStatus; group.status]; %#ok<AGROW>
end
series.x = flatX;
series.values = flatValues;
series.valid = flatValid;
series.truncated = flatTruncated;
series.power_limited = flatPowerLimited;
series.wheel_lift = flatWheelLift;
series.aero_outside_map = flatAeroOutsideMap;
series.reason = flatReason;
series.status = flatStatus;
series.groups = groups;
end

function [values,hasValues] = metricValues(T,metric)
hasValues = true(height(T),1);
try
    if isa(metric.derivation,'function_handle')
        values = metric.derivation(T);
    elseif strlength(string(metric.field)) > 0 && ...
            ismember(char(metric.field),T.Properties.VariableNames)
        values = T.(char(metric.field));
    else
        values = NaN(height(T),1);
        hasValues(:) = false;
    end
catch
    values = NaN(height(T),1);
    hasValues(:) = false;
end
if isstring(values) || iscell(values)
    values = str2double(string(values));
end
values = double(values(:));
if numel(values) ~= height(T)
    resized = NaN(height(T),1);
    n = min(numel(values),height(T));
    resized(1:n) = values(1:n);
    values = resized;
    hasValues(:) = false;
    hasValues(1:n) = true;
end
end

function values = scaledValues(values,scale)
values = double(values(:));
if isempty(scale)
    scale = 1;
end
values = values .* double(scale);
end

function mask = validityMask(T,rule)
if isa(rule,'function_handle')
    try
        mask = logical(rule(T));
    catch
        mask = true(height(T),1);
    end
else
    mask = true(height(T),1);
end
mask = mask(:);
if numel(mask) ~= height(T)
    mask = true(height(T),1);
end
end

function values = numericColumn(T,name,default)
if ismember(char(name),T.Properties.VariableNames)
    values = T.(char(name));
else
    values = repmat(double(default),height(T),1);
end
values = double(values(:));
end

function values = logicalColumn(T,name,default)
if ismember(char(name),T.Properties.VariableNames)
    values = logical(T.(char(name)));
else
    values = repmat(logical(default),height(T),1);
end
values = values(:);
end

function values = stringColumn(T,name,default)
if ismember(char(name),T.Properties.VariableNames)
    values = string(T.(char(name)));
else
    values = repmat(string(default),height(T),1);
end
values = values(:);
end

function series = seriesTemplate()
series = struct('id',"",'setupId',"",'label',"",'sourceLevel',"", ...
    'x',zeros(0,1),'values',zeros(0,1),'valid',false(0,1), ...
    'truncated',false(0,1),'power_limited',false(0,1), ...
    'wheel_lift',false(0,1),'aero_outside_map',false(0,1), ...
    'reason',strings(0,1),'status',strings(0,1), ...
    'groups',repmat(groupTemplate(),0,1));
end

function group = groupTemplate()
group = struct('speedIndex',NaN,'x',zeros(0,1),'values',zeros(0,1), ...
    'valid',false(0,1),'truncated',false(0,1), ...
    'power_limited',false(0,1),'wheel_lift',false(0,1), ...
    'aero_outside_map',false(0,1),'reason',strings(0,1), ...
    'status',strings(0,1));
end
