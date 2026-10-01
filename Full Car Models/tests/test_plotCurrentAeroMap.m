function tests = test_plotCurrentAeroMap
tests = functiontests(localfunctions);
end

function testPlotsCurrentMapMetricsOverRideHeightEnvelope(testCase)
modelRoot = fileparts(which('carConfig'));
mapPath = fullfile(modelRoot,'aeromap_b26.csv');

[fig,data] = plotCurrentAeroMap(struct( ...
    'mapPath',mapPath, ...
    'visible','off', ...
    'gridSize',[9 11]));
cleanup = onCleanup(@() close(fig)); %#ok<NASGU>

verifyTrue(testCase,isgraphics(fig));
verifyEqual(testCase,string(data.sourcePath),string(mapPath));
verifySize(testCase,data.frontOffsetIn,[9 1]);
verifySize(testCase,data.rearOffsetIn,[1 11]);
verifySize(testCase,data.CL,[9 11]);
verifySize(testCase,data.CD,[9 11]);
verifySize(testCase,data.CoP,[9 11]);
verifyTrue(testCase,all(isfinite(data.CL(:))));
verifyTrue(testCase,all(isfinite(data.CD(:))));
verifyTrue(testCase,all(isfinite(data.CoP(:))));
verifyGreaterThanOrEqual(testCase,min(data.CoP(:)),0);
verifyLessThanOrEqual(testCase,max(data.CoP(:)),1);

axesHandles = findall(fig,'Type','axes');
verifyEqual(testCase,numel(axesHandles),3);
titles = strings(numel(axesHandles),1);
for i = 1:numel(axesHandles)
    titles(i) = string(axesHandles(i).Title.String);
end
verifyTrue(testCase,any(contains(titles,'CL')));
verifyTrue(testCase,any(contains(titles,'CD')));
verifyTrue(testCase,any(contains(titles,'CoP')));

for ax = reshape(axesHandles,1,[])
    verifyEqual(testCase,string(ax.XLabel.String), ...
        "front ride-height offset (in)");
    verifyEqual(testCase,string(ax.YLabel.String), ...
        "rear ride-height offset (in)");
    verifyEqual(testCase,string(ax.ZLabel.String), ...
        string(ax.Title.String));
    verifyEqual(testCase,numel(findall(ax,'Type','surface')),1);
end
end

function testUsesRawCLAndCDColumnsAndConvertsCopPercentToFraction(testCase)
modelRoot = fileparts(which('carConfig'));
mapPath = fullfile(modelRoot,'aeromap_b26.csv');
T = readtable(mapPath,'VariableNamingRule','preserve');
row = find(double(T.('FFR offset')) == 0 & ...
    double(T.('RRH offset')) == 0,1);
verifyNotEmpty(testCase,row);

[fig,data] = plotCurrentAeroMap(struct( ...
    'mapPath',mapPath, ...
    'visible','off', ...
    'gridSize',[101 101]));
cleanup = onCleanup(@() close(fig)); %#ok<NASGU>

frontIndex = find(data.frontOffsetIn == 0,1);
rearIndex = find(data.rearOffsetIn == 0,1);
verifyNotEmpty(testCase,frontIndex);
verifyNotEmpty(testCase,rearIndex);
verifyEqual(testCase,data.CL(frontIndex,rearIndex),double(T.CL(row)), ...
    'AbsTol',1e-12);
verifyEqual(testCase,data.CD(frontIndex,rearIndex),double(T.CD(row)), ...
    'AbsTol',1e-12);
verifyEqual(testCase,data.CoP(frontIndex,rearIndex),double(T.COP(row))/100, ...
    'AbsTol',1e-12);
end
