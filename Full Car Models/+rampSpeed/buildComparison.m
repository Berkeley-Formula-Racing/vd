function comparison = buildComparison(runs,metricId,baselineId,options)
%BUILDCOMPARISON Compare each variant with a named baseline on a finite grid.

if nargin < 4 || isempty(options)
    options = struct();
end
plotData = rampSpeed.buildPlotData(runs,metricId,options);
baselineId = string(baselineId);
series = plotData.series;
baselineIndex = find(strcmp(string({series.setupId}),baselineId),1,'first');
if isempty(baselineIndex)
    error('rampSpeed:baselineNotFound', ...
        'Baseline setup "%s" was not found.',baselineId);
end

grid = comparisonGrid(series,baselineIndex,options);
comparison = struct('metricId',plotData.metricId,'metric',plotData.metric, ...
    'sourceLevel',plotData.sourceLevel,'xLabel',plotData.xLabel, ...
    'yLabel',plotData.yLabel,'baselineId',baselineId, ...
    'comparisonGrid_mps',grid,'interpolationMethod',"linear", ...
    'interpolationPolicy',"linear within contiguous valid non-truncated intervals", ...
    'variantMinusBaseline',true, ...
    'series',repmat(comparisonSeriesTemplate(),0,1));

baseline = series(baselineIndex);
[baselineValues,baselineValid] = interpolateSeries(baseline,grid);
baselineTruncated = resampleMask(baseline,"truncated",grid);
baselinePowerLimited = resampleMask(baseline,"power_limited",grid);
baselineWheelLift = resampleMask(baseline,"wheel_lift",grid);
for i = 1:numel(series)
    if i == baselineIndex
        continue
    end
    [variantValues,variantValid] = interpolateSeries(series(i),grid);
    values = variantValues - baselineValues;
    valid = baselineValid & variantValid & isfinite(values);
    values(~valid) = NaN;
    result = comparisonSeriesTemplate();
    result.id = series(i).setupId;
    result.setupId = series(i).setupId;
    result.label = series(i).label;
    result.baselineId = baselineId;
    result.x = grid;
    result.values = values;
    result.variantValues = variantValues;
    result.baselineValues = baselineValues;
    result.valid = valid;
    result.variantValid = variantValid;
    result.baselineValid = baselineValid;
    result.variantTruncated = resampleMask(series(i),"truncated",grid);
    result.variantPowerLimited = resampleMask(series(i),"power_limited",grid);
    result.variantWheelLift = resampleMask(series(i),"wheel_lift",grid);
    result.baselineTruncated = baselineTruncated;
    result.baselinePowerLimited = baselinePowerLimited;
    result.baselineWheelLift = baselineWheelLift;
    result.truncated = result.baselineTruncated | result.variantTruncated;
    result.power_limited = result.baselinePowerLimited | result.variantPowerLimited;
    result.wheel_lift = result.baselineWheelLift | result.variantWheelLift;
    comparison.series(end+1,1) = result;
end
end

function grid = comparisonGrid(series,baselineIndex,options)
if isfield(options,'comparisonGrid_mps') && ~isempty(options.comparisonGrid_mps)
    grid = double(options.comparisonGrid_mps(:));
elseif isfield(options,'grid') && ~isempty(options.grid)
    grid = double(options.grid(:));
else
    grid = commonNativeGrid(series,baselineIndex);
end
grid = grid(isfinite(grid));
grid = unique(sort(grid));
end

function grid = commonNativeGrid(series,baselineIndex)
if isempty(series)
    grid = zeros(0,1);
    return
end
lower = -Inf;
upper = Inf;
hasDomain = true;
for i = 1:numel(series)
    x = finiteValues(series(i).x);
    if isempty(x)
        hasDomain = false;
        break
    end
    lower = max(lower,min(x));
    upper = min(upper,max(x));
end
if ~hasDomain || lower > upper
    grid = zeros(0,1);
    return
end
grid = zeros(0,1);
for i = 1:numel(series)
    x = finiteValues(series(i).x);
    grid = [grid; x(x >= lower & x <= upper)]; %#ok<AGROW>
end
grid = unique(sort(grid));
if isempty(grid) && baselineIndex <= numel(series)
    x = finiteValues(series(baselineIndex).x);
    grid = x(x >= lower & x <= upper);
end
end

function [values,valid] = interpolateSeries(series,grid)
values = NaN(size(grid));
valid = false(size(grid));
x = double(series.x(:));
y = double(series.values(:));
rowValid = logical(series.valid(:));
truncated = logical(series.truncated(:));
n = min([numel(x),numel(y),numel(rowValid),numel(truncated)]);
x = x(1:n);
y = y(1:n);
rowValid = rowValid(1:n);
truncated = truncated(1:n);
good = rowValid & ~truncated & isfinite(x) & isfinite(y);
if isempty(good)
    return
end
starts = find(good & [true; ~good(1:end-1)]);
ends = find(good & [~good(2:end); true]);
for i = 1:numel(starts)
    segment = starts(i):ends(i);
    segmentX = x(segment);
    segmentY = y(segment);
    [segmentX,order] = sort(segmentX);
    segmentY = segmentY(order);
    [segmentX,uniqueIndex] = unique(segmentX,'stable');
    segmentY = segmentY(uniqueIndex);
    if isscalar(segmentX)
        query = grid == segmentX;
        values(query) = segmentY;
        valid(query) = isfinite(segmentY);
    else
        query = grid >= segmentX(1) & grid <= segmentX(end);
        if any(query)
            values(query) = interp1(segmentX,segmentY,grid(query),'linear');
            valid(query) = isfinite(values(query));
        end
    end
end
end

function mask = resampleMask(series,field,grid)
mask = false(size(grid));
if ~isfield(series,field)
    return
end
x = double(series.x(:));
flags = logical(series.(field)(:));
rowValid = logical(series.valid(:));
truncated = logical(series.truncated(:));
n = min([numel(x),numel(flags),numel(rowValid),numel(truncated)]);
x = x(1:n);
flags = flags(1:n);
rowValid = rowValid(1:n);
truncated = truncated(1:n);
exact = flags & isfinite(x);
for i = find(exact).'
    mask(grid == x(i)) = true;
end
good = rowValid & ~truncated & isfinite(x);
if ~any(good)
    return
end
starts = find(good & [true; ~good(1:end-1)]);
ends = find(good & [~good(2:end); true]);
for i = 1:numel(starts)
    segment = starts(i):ends(i);
    segmentX = x(segment);
    segmentFlags = double(flags(segment));
    [segmentX,order] = sort(segmentX);
    segmentFlags = segmentFlags(order);
    [segmentX,uniqueIndex] = unique(segmentX,'stable');
    segmentFlags = segmentFlags(uniqueIndex);
    if isscalar(segmentX)
        query = grid == segmentX;
        mask(query) = mask(query) | logical(segmentFlags);
    else
        query = grid >= segmentX(1) & grid <= segmentX(end);
        if any(query)
            interpolated = interp1(segmentX,segmentFlags,grid(query),'linear');
            mask(query) = mask(query) | interpolated >= 0.5;
        end
    end
end
end

function x = finiteValues(x)
x = double(x(:));
x = x(isfinite(x));
end

function result = comparisonSeriesTemplate()
result = struct('id',"",'setupId',"",'label',"",'baselineId',"", ...
    'x',zeros(0,1),'values',zeros(0,1),'variantValues',zeros(0,1), ...
    'baselineValues',zeros(0,1),'valid',false(0,1), ...
    'variantValid',false(0,1),'baselineValid',false(0,1), ...
    'truncated',false(0,1),'power_limited',false(0,1), ...
    'wheel_lift',false(0,1), ...
    'baselineTruncated',false(0,1),'variantTruncated',false(0,1), ...
    'baselinePowerLimited',false(0,1),'variantPowerLimited',false(0,1), ...
    'baselineWheelLift',false(0,1),'variantWheelLift',false(0,1));
end
