function [fig,data] = plotCurrentAeroMap(options)
%PLOTCURRENTAEROMAP Plot the current CFD aeromap over front/rear ride heights.
%   [FIG,DATA] = PLOTCURRENTAEROMAP() loads the project aeromap and opens
%   three 3-D surfaces for CL, CD, and front CoP fraction.
%
%   Optional fields in OPTIONS are:
%       mapPath  - path to an aeromap CSV (default: the catalog B26 map)
%       visible  - figure visibility, 'on' or 'off' (default: 'on')
%       gridSize - [nFront nRear] query-grid size (default: [31 31])
%
%   The horizontal axes are CSV FFR/RRH ride-height offsets in inches. The
%   CoP column is stored as a percentage in the CSV and is returned/plotted
%   as a front-load fraction, matching AeroMap.evaluate().

if nargin < 1 || isempty(options)
    options = struct();
end
if ~isstruct(options) || ~isscalar(options)
    error('plotCurrentAeroMap:invalidOptions', ...
        'options must be a scalar struct.');
end

modelRoot = fileparts(fileparts(mfilename('fullpath')));
defaultMapPath = fullfile(modelRoot,'aeromap_b26.csv');
mapPath = optionValue(options,'mapPath',defaultMapPath);
visible = optionValue(options,'visible','on');
gridSize = optionValue(options,'gridSize',[31 31]);

mapPath = char(string(mapPath));
visible = lower(char(string(visible)));
if ~ismember(visible,{'on','off'})
    error('plotCurrentAeroMap:invalidVisibility', ...
        'visible must be ''on'' or ''off''.');
end

validateattributes(gridSize,{'numeric'}, ...
    {'vector','numel',2,'integer','>=',2},mfilename,'gridSize');
gridSize = double(reshape(gridSize,1,2));
if ~isfile(mapPath)
    error('plotCurrentAeroMap:mapNotFound','Aeromap file not found: %s',mapPath);
end

T = readtable(mapPath,'VariableNamingRule','preserve');
required = ["FFR offset" "RRH offset" "CL" "CD" "COP"];
if ~all(ismember(required,string(T.Properties.VariableNames)))
    error('plotCurrentAeroMap:missingColumns', ...
        'Aeromap must contain FFR offset, RRH offset, CL, CD, and COP columns.');
end

front = double(T.('FFR offset'));
rear = double(T.('RRH offset'));
cl = double(T.CL);
cd = double(T.CD);
copPercent = double(T.COP);
valid = isfinite(front) & isfinite(rear) & isfinite(cl) & ...
    isfinite(cd) & isfinite(copPercent);
if nnz(valid) < 3
    error('plotCurrentAeroMap:notEnoughPoints', ...
        'Aeromap needs at least three finite CL/CD/COP samples.');
end

front = front(valid);
rear = rear(valid);
cl = cl(valid);
cd = cd(valid);
copPercent = copPercent(valid);

clInterpolator = scatteredInterpolant(front,rear,cl,'linear','nearest');
cdInterpolator = scatteredInterpolant(front,rear,cd,'linear','nearest');

% Keep ride-height bounds and CoP behavior aligned with the production model.
aeroMap = AeroMap(mapPath);
frontOffsetIn = linspace(aeroMap.frontRangeIn(1), ...
    aeroMap.frontRangeIn(2),gridSize(1)).';
rearOffsetIn = linspace(aeroMap.rearRangeIn(1), ...
    aeroMap.rearRangeIn(2),gridSize(2));
[frontQuery,rearQuery] = ndgrid(frontOffsetIn,rearOffsetIn);

clValues = clInterpolator(frontQuery,rearQuery);
cdValues = cdInterpolator(frontQuery,rearQuery);
[~,~,copValues] = aeroMap.evaluateNumeric(frontQuery,rearQuery);

fig = figure('Name','Current Aero Map','Color','w','Visible',visible, ...
    'Position',[100 80 1100 1500]);
layout = tiledlayout(fig,3,1,'TileSpacing','compact','Padding','compact');
metrics = { ...
    'CL',clValues; ...
    'CD',cdValues; ...
    'CoP (front fraction)',copValues};
measuredValues = {cl,cd,copPercent/100};

for metricIndex = 1:size(metrics,1)
    ax = nexttile(layout);
    surf(ax,frontQuery,rearQuery,metrics{metricIndex,2}, ...
        'EdgeColor','none','FaceColor','interp');
    hold(ax,'on');
    plot3(ax,front,rear,measuredValues{metricIndex}, ...
        'k.','MarkerSize',12,'DisplayName','measured samples');
    hold(ax,'off');
    grid(ax,'on');
    box(ax,'on');
    view(ax,3);
    xlabel(ax,'front ride-height offset (in)');
    ylabel(ax,'rear ride-height offset (in)');
    zlabel(ax,metrics{metricIndex,1});
    title(ax,metrics{metricIndex,1});
    set(ax,'FontSize',10);
end
sgtitle(layout,sprintf('Current aeromap: %s',string(mapPath)), ...
    'Interpreter','none');

data = struct();
data.sourcePath = string(mapPath);
data.frontOffsetIn = frontOffsetIn;
data.rearOffsetIn = rearOffsetIn;
data.CL = clValues;
data.CD = cdValues;
data.CoP = copValues;
data.measuredFrontOffsetIn = front;
data.measuredRearOffsetIn = rear;
data.measuredCL = cl;
data.measuredCD = cd;
data.measuredCoP = copPercent/100;
end

function value = optionValue(options,name,defaultValue)
if isfield(options,name) && ~isempty(options.(name))
    value = options.(name);
else
    value = defaultValue;
end
end
