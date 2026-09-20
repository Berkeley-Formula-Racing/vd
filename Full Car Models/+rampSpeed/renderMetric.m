function handles = renderMetric(ax,data,options)
%RENDERMETRIC Render prepared metric data without solver-object dependencies.

if nargin < 3 || isempty(options)
    options = struct();
end
if isempty(ax) || ~isgraphics(ax,'axes')
    error('rampSpeed:invalidAxes','ax must be a valid axes handle.');
end

handles = struct('axes',ax,'lines',gobjects(0,1), ...
    'warningMarkers',gobjects(0,1),'zeroLine',[]);
series = data.series;
colors = colorSet(options,max(numel(series),1));
legendHandles = gobjects(0,1);
legendLabels = strings(0,1);

for i = 1:numel(series)
    color = colors(min(i,size(colors,1)),:);
    if isfield(series(i),'groups') && ~isempty(series(i).groups)
        for j = 1:numel(series(i).groups)
            group = series(i).groups(j);
            if isempty(group.x)
                continue
            end
            lineHandle = drawSeries(ax,group.x,group.values,data.metric.seriesStyle,color);
            handles.lines(end+1,1) = lineHandle;
            if j == 1
                legendHandles(end+1,1) = lineHandle; %#ok<AGROW>
                legendLabels(end+1,1) = series(i).label; %#ok<AGROW>
            else
                set(lineHandle,'HandleVisibility','off');
            end
        end
    elseif ~isempty(series(i).x)
        lineHandle = drawSeries(ax,series(i).x,series(i).values, ...
            data.metric.seriesStyle,color);
        handles.lines(end+1,1) = lineHandle;
        legendHandles(end+1,1) = lineHandle; %#ok<AGROW>
        legendLabels(end+1,1) = series(i).label; %#ok<AGROW>
    end
    if showWarnings(options)
        handles.warningMarkers = [handles.warningMarkers; ... %#ok<AGROW>
            warningMarkers(ax,series(i))];
    end
end

if isfield(data.metric,'zeroLine') && data.metric.zeroLine
    handles.zeroLine = yline(ax,0,'--','HandleVisibility','off');
end
xlabel(ax,data.xLabel);
ylabel(ax,data.yLabel);
title(ax,metricTitle(data.metric));
if showLegend(options) && ~isempty(legendHandles)
    legend(ax,legendHandles,legendLabels,'Location','best');
end
grid(ax,'on');
end

function lineHandle = drawSeries(ax,x,y,style,color)
lineStyle = "-";
marker = "none";
lineWidth = 1.5;
markerSize = 6;
plotStyle = "";
if isstruct(style)
    if isfield(style,'lineStyle') && ~isempty(style.lineStyle)
        lineStyle = string(style.lineStyle);
    end
    if isfield(style,'marker') && ~isempty(style.marker)
        marker = string(style.marker);
    end
    if isfield(style,'lineWidth') && ~isempty(style.lineWidth)
        lineWidth = double(style.lineWidth);
    end
    if isfield(style,'markerSize') && ~isempty(style.markerSize)
        markerSize = double(style.markerSize);
    end
    if isfield(style,'style')
        plotStyle = lower(string(style.style));
    end
end
if plotStyle == "stairs"
    lineHandle = stairs(ax,x,y,'Color',color,'LineStyle',lineStyle, ...
        'LineWidth',lineWidth,'Marker',marker,'MarkerSize',markerSize);
else
    lineHandle = plot(ax,x,y,'Color',color,'LineStyle',lineStyle, ...
        'LineWidth',lineWidth,'Marker',marker,'MarkerSize',markerSize);
end
end

function markers = warningMarkers(ax,series)
markers = gobjects(0,1);
x = double(series.x(:));
y = double(series.values(:));
markerY = y;
markerY(~isfinite(markerY)) = 0;
masks = {logical(series.truncated(:)),logical(series.power_limited(:)), ...
    logical(series.wheel_lift(:)),~logical(series.valid(:))};
styles = {'o',[0.90 0.45 0.05],'^',[0.85 0.65 0.05], ...
    's',[0.75 0.10 0.10],'x',[0.35 0.35 0.35]};
for i = 1:2:numel(masks)
    mask = masks{i} & isfinite(x);
    if ~any(mask)
        continue
    end
    marker = styles{i};
    color = styles{i+1};
    if i == 7
        marker = 'x';
        color = styles{8};
    end
    h = plot(ax,x(mask),markerY(mask),marker,'Color',color, ...
        'LineStyle','none','MarkerSize',7,'LineWidth',1.2, ...
        'HandleVisibility','off');
    markers(end+1,1) = h; %#ok<AGROW>
end
end

function colors = colorSet(options,n)
if isfield(options,'colors') && ~isempty(options.colors)
    colors = double(options.colors);
    if size(colors,2) ~= 3
        colors = lines(n);
    end
else
    colors = lines(n);
end
if size(colors,1) < n
    colors = lines(n);
end
end

function result = showLegend(options)
if isfield(options,'showLegend') && ~isempty(options.showLegend)
    result = logical(options.showLegend);
else
    result = true;
end
end

function result = showWarnings(options)
if isfield(options,'showWarnings') && ~isempty(options.showWarnings)
    result = logical(options.showWarnings);
else
    result = true;
end
end

function titleText = metricTitle(metric)
titleText = string(metric.title);
if isfield(metric,'subtitle') && strlength(string(metric.subtitle)) > 0
    titleText = titleText + newline + string(metric.subtitle);
end
end
