%% SteadyStateLapsim - simple LC0 event runner
clear classes
entrypointDir = fileparts(mfilename('fullpath'));
run(fullfile(entrypointDir, 'bootstrap.m'));

% Edit this list to choose which dynamic events to run for every car.
eventsToRun = ["skidpad","accel","autocross"];

% carConfig is the single source of vehicle and event configuration.  The
% explicit tire-model argument selects the new nondimensional LC0 model.
[carCell,eventParams] = carConfig("FullFactorial",[],"legacy");
ggOptions = ggProductionOptions();
ggWorkers = productionGgWorkers(numel(ggGrid(ggOptions.maxVelocity,ggOptions).velocity));

for caseIndex = 1:size(carCell,1)
    job = simLog.start('single car run');
    car = carCell{caseIndex,1};
    accelCar = carCell{caseIndex,2};

    car = makeGG(gg2(car,ggWorkers,ggOptions),car);
    eventSim = Events2(car,accelCar,eventParams);

    for eventIndex = 1:numel(eventsToRun)
        eventName = char(eventsToRun(eventIndex));
        solverName = eventName;
        solverName(1) = upper(solverName(1));

        fprintf('Case %d/%d: running %s...\n', ...
            caseIndex,size(carCell,1),eventName);
        eventSim.(solverName)();
        fprintf('  %s = %.4f s\n',eventName,eventSim.times.(eventName));
    end

    car.comp = eventSim;
    carCell{caseIndex,1} = car;
    simLog.finish(job, ...
        'events', cellstr(eventsToRun), ...
        'nCases', 1, ...
        'workers', ggWorkers, ...
        'times', eventSim.times, ...
        'car', car);
end

function numWorkers = productionGgWorkers(rowCount)
numWorkers = 0;
if ~license('test','Distrib_Computing_Toolbox')
    return
end

requestedWorkers = min(feature('numcores'),rowCount);
try
    pool = gcp('nocreate');
    if isempty(pool)
        pool = parpool(requestedWorkers);
    end
    numWorkers = min(pool.NumWorkers,rowCount);
catch ME
    warning('SteadyStateLapsim:parallelPoolUnavailable', ...
        'Falling back to serial G-G solve: %s',ME.message)
end
end
