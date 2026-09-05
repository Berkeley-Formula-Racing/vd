%% AeroMapSensitivityStudy
% Re-solve QSS performance while perturbing static ride height and the CFD
% aeromap itself. The five star-design inputs are front/rear ride height,
% map ClA scale, map CdA scale, and map CoP offset.

clear classes
setup_paths

runLapTimes = true;
loadFromPrev = false;
numWorkers = 12;
loadName = "aeromap_sensitivity.mat";
saveName = "aeromap_sensitivity.mat";
figPath = "figures/aeromap_sensitivity";
figOpts = struct('formats',{{'png'}},'resolution',200,'stamp',false);

% Ride-height values are inch changes from the map-reference baseline.
% ClA/CdA are multiplicative CFD-map corrections; CoP is an additive front
% downforce-fraction correction applied after interpolating the map.
levels = struct();
levels.FrontRideHeightIn = [-0.25 -0.125 0.125 0.25];
levels.RearRideHeightIn  = [-0.25 -0.125 0.125 0.25];
levels.ClAScale = [0.90 0.95 1.05 1.10];
levels.CdAScale = [0.90 0.95 1.05 1.10];
levels.CoPOffset = [-0.06 -0.03 0.03 0.06];

[configuredCars,eventParams] = carConfig();
[carCell,plan] = aeroMapStarCases(configuredCars,levels);
mapFingerprint = fingerprint(carCell{1,1});
numCases = size(carCell,1);
fprintf('aeromap sensitivity study: %d independent cases\n',numCases);
disp(plan)

tic
job = simLog.start('aeromap sensitivity study');
haveGG = false; haveTimes = false; times = []; traces = []; %#ok<NASGU>
if loadFromPrev
    [cached,times,traces,haveGG,haveTimes] = loadCache(loadName,plan,mapFingerprint);
    if haveGG
        carCell = cached;
        fprintf('  loaded %s: g-g=%d, lap times=%d\n',loadName,haveGG,haveTimes);
    end
end

if ~haveGG
    fprintf('  g-g for %d cases ...\n',numCases);
    solved = cell(numCases,1);
    parfor (i = 1:numCases,numWorkers)
        solved{i} = makeGG(gg2(carCell{i,1},0),carCell{i,1});
    end
    for i = 1:numCases, carCell{i,1} = solved{i}; end
    fprintf('  g-g done (%.0f s elapsed)\n',toc)
end

opts = struct();
if runLapTimes
    if ~haveTimes
        % The acceleration package is static and unchanged by this study, so
        % run it once at baseline rather than repeating the same event 21 times.
        eventList = {'skidpad','autocross','endurance'};
        [times,traces] = aeroLapTimes(carCell,numWorkers,eventList,eventParams);
        [accelTime,~] = aeroLapTimes(carCell(1,:),0,{'accel'},eventParams);
        times.accel(:) = accelTime.accel(1);
        fprintf('  events done (%.0f s elapsed)\n',toc)
    end
    opts.times = times;
    opts.traces = traces;
    opts.fastestBy = 'autocross';
end

if ~haveGG || (runLapTimes && ~haveTimes)
    save(saveName,'carCell','plan','mapFingerprint','-v7.3');
    if runLapTimes, save(saveName,'times','traces','-append'); end
end

S = aeroMapSensitivity(carCell,plan,opts);
aeroMapSensitivityReport(S);
fprintf('done in %.0f s\n',toc)

close all
parameters = fieldnames(S.sens);
for i = 1:numel(parameters)
    parameter = parameters{i};
    plotAeroSensitivity(S,parameter,'both');
    if isfield(S,'traces') && ~isempty(S.traces)
        plotAeroKernel(S,'autocross',parameter);
        plotAeroOverlay(S,parameter,'autocross',12,'seconds');
    end
end
saveFigures(figPath,[],figOpts)

events = 'g-g';
if runLapTimes, events = 'g-g+skidpad+accel+autocross+endurance'; end
simLog.finish(job,'events',events,'workers',numWorkers,'nCases',numCases, ...
    'car',carCell{1,1},'details',sprintf('cache=%s map=%s', ...
    char(saveName),char(mapFingerprint.path)));

function f = fingerprint(car)
path = char(car.aero.map.sourcePath);
file = dir(path);
if isempty(file), error('AeroMapSensitivityStudy:mapMissing','cannot read %s',path); end
f = struct('path',string(path),'bytes',file.bytes,'modified',file.datenum);
end

function [cars,times,traces,haveGG,haveTimes] = loadCache(name,plan,fingerprint)
cars = []; times = []; traces = []; haveGG = false; haveTimes = false;
if ~isfile(name), return, end
cached = load(name);
needed = {'carCell','plan','mapFingerprint'};
if ~all(isfield(cached,needed)) || ~isequaln(cached.plan,plan) || ...
        ~isequaln(cached.mapFingerprint,fingerprint)
    warning('AeroMapSensitivityStudy:staleCache', ...
        'cache does not match the current star design or aeromap CSV; re-solving.');
    return
end
cars = cached.carCell;
haveGG = all(cellfun(@(car) ~isempty(car.ss_info),cars(:,1)));
if isfield(cached,'times') && isfield(cached,'traces')
    times = cached.times; traces = cached.traces; haveTimes = haveGG;
end
end
