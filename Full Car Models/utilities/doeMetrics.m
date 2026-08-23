function [M,details] = doeMetrics(carCell,opts)
%DOEMETRICS Reduce solved DOE cars to one design-analysis row per case.
%   Uses saved g-g states and event traces by default; no vehicle solves are
%   repeated. Set opts.runRampMetrics=true to add true ramp-steer gradients.

if nargin < 2, opts = struct(); end
runRamp = getOr(opts,'runRampMetrics',false);
rampOpts = getOr(opts,'rampOpts',struct('speeds',[10 25], ...
    'nRamp',8,'nBisect',3,'verbose',false));

n = size(carCell,1);
rows = cell(n,1);
details = struct('ramp',{cell(n,1)});

for i = 1:n
    try
        if size(carCell,2) < 1 || isempty(carCell{i,1})
            [row,caseDetails] = doeCaseMetrics([],i,[]);
        else
            car = carCell{i,1};
            rampResult = [];
            if runRamp
                rampResult = rampSweep(car,rampOpts);
            end
            [row,caseDetails] = doeCaseMetrics(car,i,rampResult);
        end
    catch ME
        [row,caseDetails] = doeCaseMetrics([],i,[]);
        row.valid = false;
        row.error_message = string(ME.message);
    end
    rows{i} = table2struct(row);
    details.ramp{i} = caseDetails.ramp;
end

M = struct2table(vertcat(rows{:}));
end

function value = getOr(s,name,default)
if isfield(s,name), value=s.(name); else, value=default; end
end
