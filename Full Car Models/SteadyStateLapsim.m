%% SteadyStateLapsim - DOE runner
% Samples complete car configurations, solves each g-g diagram serially
% inside an outer parallel loop, and runs skidpad, accel and autocross.

clear classes
setup_paths

%% DOE settings
samplingType = "LHS";       % "LHS", "Random", or "FullFactorial"
numSamples   = 512;         % ignored for FullFactorial
numWorkers   = 16;          % parallel cars; each gg2 remains serial
randomSeed   = 1;           % repeatable sample locations
savePath     = "DOE_results.mat";

rng(randomSeed);
[carCell,eventParams,designTable] = carConfig(samplingType,numSamples);
numCars = size(carCell,1);
numWorkers = min(numWorkers,numCars);

fprintf('Starting %s DOE: %d cars on %d workers.\n', ...
    samplingType,numCars,numWorkers);
job = simLog.start('vehicle DOE');
started = tic;

carOutMain = cell(numCars,1);
carOutAccel = cell(numCars,1);
parfor (i = 1:numCars,numWorkers)
    car = carCell{i,1};
    accelCar = carCell{i,2};

    rawGG = gg2(car,0);
    car = makeGG(rawGG,car);

    comp = Events2(car,accelCar,eventParams);
    comp.Skidpad();
    comp.Accel();
    comp.Autocross();
    car.comp = comp;

    carOutMain{i} = car;
    carOutAccel{i} = accelCar;
    fprintf('Completed DOE car %d of %d.\n',i,numCars);
end

carCell = [carOutMain,carOutAccel];
elapsed = toc(started);
save(savePath,'carCell','designTable','eventParams','samplingType', ...
    'numSamples','numWorkers','randomSeed','elapsed','-v7.3');
simLog.finish(job,'events',{'skidpad','accel','autocross'}, ...
    'workers',numWorkers,'nCases',numCars, ...
    'details',sprintf('%s seed=%d output=%s',samplingType,randomSeed,savePath));

fprintf('DOE complete in %.1f min (%.1f worker-min/car equivalent).\n', ...
    elapsed/60,elapsed*numWorkers/numCars/60);
