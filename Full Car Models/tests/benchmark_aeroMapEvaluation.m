function report = benchmark_aeroMapEvaluation(queryCount)
%BENCHMARK_AEROMAPEVALUATION Compare scalar struct and numeric batch queries.
% Run from Full Car Models after adding this folder to the MATLAB path.
if nargin < 1, queryCount = 5000; end
validateattributes(queryCount,{'numeric'},{'scalar','integer','positive'});

modelRoot = fileparts(which('carConfig'));
map = AeroMap(fullfile(modelRoot,'aeromap_b26.csv'));
front = linspace(map.frontRangeIn(1)-0.1,map.frontRangeIn(2)+0.1,queryCount);
rear = linspace(map.rearRangeIn(2)+0.1,map.rearRangeIn(1)-0.1,queryCount);

batchSeconds = timeit(@() batchQuery(map,front,rear));
scalarSeconds = timeit(@() scalarQueries(map,front,rear));
report = struct('queryCount',queryCount, ...
    'interpolationMode',map.interpolationMode, ...
    'batchSeconds',batchSeconds, ...
    'scalarSeconds',scalarSeconds, ...
    'speedup',scalarSeconds/batchSeconds);
fprintf('AeroMap %s: batch %.6g s, scalar %.6g s, speedup %.2fx (%d queries)\n', ...
    report.interpolationMode,report.batchSeconds,report.scalarSeconds, ...
    report.speedup,report.queryCount);
end

function batchQuery(map,front,rear)
[~,~,~,~,~] = map.evaluateNumeric(front,rear);
end

function scalarQueries(map,front,rear)
for k = 1:numel(front)
    map.evaluate(front(k),rear(k));
end
end
