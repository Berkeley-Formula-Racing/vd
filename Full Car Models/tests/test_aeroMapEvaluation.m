function tests = test_aeroMapEvaluation
% Focused regression coverage for batched AeroMap interpolation.
tests = functiontests(localfunctions);
end

function testBatchNumericMethodPreservesExistingScatteredResults(testCase)
modelRoot = fileparts(which('carConfig'));
csvPath = fullfile(modelRoot,'aeromap_b26.csv');
map = AeroMap(csvPath).withCorrections(1.10,0.90,0.04);
front = [0 -0.30 0.45 -2 0.22];
rear = [0 0.20 -0.30 2 -0.11];

[cla,cda,D_f,D_r,outsideMap] = map.evaluateNumeric(front,rear);
[expectedCla,expectedCda,expectedDf,expectedDr,expectedOutside] = ...
    scatteredReference(csvPath,front,rear,1.10,0.90,0.04);

verifyEqual(testCase,size(cla),size(front));
verifyEqual(testCase,size(cda),size(front));
verifyEqual(testCase,size(D_f),size(front));
verifyEqual(testCase,size(D_r),size(front));
verifyEqual(testCase,cla,expectedCla,'AbsTol',1e-12);
verifyEqual(testCase,cda,expectedCda,'AbsTol',1e-12);
verifyEqual(testCase,D_f,expectedDf,'AbsTol',1e-12);
verifyEqual(testCase,D_r,expectedDr,'AbsTol',1e-12);
verifyEqual(testCase,outsideMap,expectedOutside);
% The third query is beyond the CSV's measured front-offset range.
verifyEqual(testCase,outsideMap,[false false true true false]);
end

function testBatchNumericMethodSupportsScalarExpansion(testCase)
modelRoot = fileparts(which('carConfig'));
map = AeroMap(fullfile(modelRoot,'aeromap_b26.csv'));

[cla,cda,D_f,D_r,outsideMap] = map.evaluateNumeric(0,[-0.30 0 0.20]);

verifySize(testCase,cla,[1 3]);
verifySize(testCase,cda,[1 3]);
verifySize(testCase,D_f,[1 3]);
verifySize(testCase,D_r,[1 3]);
verifyEqual(testCase,outsideMap,false(1,3));
end

function testRectangularPlanarMapUsesGridAndMatchesScatteredReference(testCase)
csvPath = writeSyntheticMap(testCase,true,true);
map = AeroMap(csvPath);
front = [-0.73 0.22 -2.0 1.4];
rear = [-1.37 1.82 0.3 0.2];

verifyEqual(testCase,map.interpolationMode,"gridded");
[cla,cda,D_f,D_r,outsideMap] = map.evaluateNumeric(front,rear);
[expectedCla,expectedCda,expectedDf,expectedDr,expectedOutside] = ...
    scatteredReference(csvPath,front,rear);

verifyEqual(testCase,cla,expectedCla,'AbsTol',1e-12);
verifyEqual(testCase,cda,expectedCda,'AbsTol',1e-12);
verifyEqual(testCase,D_f,expectedDf,'AbsTol',1e-12);
verifyEqual(testCase,D_r,expectedDr,'AbsTol',1e-12);
verifyEqual(testCase,outsideMap,expectedOutside);
verifyTrue(testCase,outsideMap(3));
end

function testNonPlanarRectangularMapFallsBackToScatteredInterpolation(testCase)
csvPath = writeSyntheticMap(testCase,true,false);
map = AeroMap(csvPath);
front = [-0.7 0.15 0.8];
rear = [-1.6 0.4 1.8];

verifyEqual(testCase,map.interpolationMode,"scattered");
[cla,cda,D_f,D_r,outsideMap] = map.evaluateNumeric(front,rear);
[expectedCla,expectedCda,expectedDf,expectedDr,expectedOutside] = ...
    scatteredReference(csvPath,front,rear);

verifyEqual(testCase,cla,expectedCla,'AbsTol',1e-12);
verifyEqual(testCase,cda,expectedCda,'AbsTol',1e-12);
verifyEqual(testCase,D_f,expectedDf,'AbsTol',1e-12);
verifyEqual(testCase,D_r,expectedDr,'AbsTol',1e-12);
verifyEqual(testCase,outsideMap,expectedOutside);
end

function testIncompleteRectangularMapFallsBackToScatteredInterpolation(testCase)
csvPath = writeSyntheticMap(testCase,false,true);
map = AeroMap(csvPath);
front = [0.2 -0.5];
rear = [0.3 -1.1];

verifyEqual(testCase,map.interpolationMode,"scattered");
[cla,cda,D_f,D_r,outsideMap] = map.evaluateNumeric(front,rear);
[expectedCla,expectedCda,expectedDf,expectedDr,expectedOutside] = ...
    scatteredReference(csvPath,front,rear);

verifyEqual(testCase,cla,expectedCla,'AbsTol',1e-12);
verifyEqual(testCase,cda,expectedCda,'AbsTol',1e-12);
verifyEqual(testCase,D_f,expectedDf,'AbsTol',1e-12);
verifyEqual(testCase,D_r,expectedDr,'AbsTol',1e-12);
verifyEqual(testCase,outsideMap,expectedOutside);
end

function testEvaluateKeepsScalarStructContractAndEnvelopeFlag(testCase)
modelRoot = fileparts(which('carConfig'));
map = AeroMap(fullfile(modelRoot,'aeromap_b26.csv'));

aero = map.evaluate(-2,2);

verifyEqual(testCase,aero.frontOffsetIn,-2);
verifyEqual(testCase,aero.rearOffsetIn,2);
verifyTrue(testCase,aero.outsideMap);
verifyEqual(testCase,aero.D_r,1-aero.D_f,'AbsTol',1e-12);
verifyTrue(testCase,isfinite(aero.cla) && isfinite(aero.cda));
end

function csvPath = writeSyntheticMap(testCase,isComplete,isPlanar)
testFolder = tempname;
mkdir(testFolder);
testCase.addTeardown(@() rmdir(testFolder,'s'));

frontAxis = [-1;0;1];
rearAxis = [-2;0;2];
[frontGrid,rearGrid] = ndgrid(frontAxis,rearAxis);
front = frontGrid(:);
rear = rearGrid(:);
if isPlanar
    cla = 2.5 + 0.3*front - 0.2*rear;
    cda = 1.2 + 0.04*front + 0.03*rear;
    cop = 51 + 2*front - 1.5*rear;
else
    cla = 2.5 + 0.3*front - 0.2*rear + 0.2*front.*rear;
    cda = 1.2 + 0.04*front + 0.03*rear + 0.05*front.*rear;
    cop = 51 + 2*front - 1.5*rear + 1.2*front.*rear;
end

if ~isComplete
    keep = ~(front == 1 & rear == 2);
    front = front(keep);
    rear = rear(keep);
    cla = cla(keep);
    cda = cda(keep);
    cop = cop(keep);
end

rows = [{'Simulation Number','FFR offset','RRH offset','CLA','CDA','COP'}; ...
    num2cell([(1:numel(front))' front rear cla cda cop])];
csvPath = fullfile(testFolder,'synthetic_aeromap.csv');
writecell(rows,csvPath);
end

function [cla,cda,D_f,D_r,outsideMap] = scatteredReference( ...
        csvPath,front,rear,claScale,cdaScale,copOffset)
if nargin < 4, claScale = 1; end
if nargin < 5, cdaScale = 1; end
if nargin < 6, copOffset = 0; end
T = readtable(csvPath,'VariableNamingRule','preserve');
f = double(T.('FFR offset'));
r = double(T.('RRH offset'));
rawCla = double(T.CLA);
rawCda = double(T.CDA);
rawCop = double(T.COP);
valid = isfinite(f) & isfinite(r) & isfinite(rawCla) & ...
    isfinite(rawCda) & isfinite(rawCop);
claMap = scatteredInterpolant(f(valid),r(valid),rawCla(valid),'linear','nearest');
cdaMap = scatteredInterpolant(f(valid),r(valid),rawCda(valid),'linear','nearest');
copMap = scatteredInterpolant(f(valid),r(valid),rawCop(valid),'linear','nearest');
cla = claScale*claMap(front,rear);
cda = cdaScale*cdaMap(front,rear);
D_f = min(max(copMap(front,rear)/100 + copOffset,0),1);
D_r = 1-D_f;
outsideMap = front < min(f(valid)) | front > max(f(valid)) | ...
    rear < min(r(valid)) | rear > max(r(valid));
end
