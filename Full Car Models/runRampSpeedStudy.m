%% runRampSpeedStudy
% Compare ride-height aero, mechanical balance, and total handling balance
% over vCar. Run this file, then edit the CARS section to compare variants.

clear classes %#ok<CLCLS> % deliberate: reload changed Car/Aero class definitions
setup_paths

%% cars to compare
[carCell,~] = carConfig();

% Default: every lap car returned by carConfig. If carConfig is a single
% baseline this gives one trace; if it contains a sweep it overlays every
% configuration automatically.
cars = carCell(:,1);
labels = "car " + string((1:numel(cars)).');

% Typical explicit comparison:
% baseline = carCell{1,1};
% variant  = baseline; variant.gripScaleF = 0.60;
% cars = {baseline; variant};
% labels = ["baseline"; "front grip 0.60"];

%% solve and save settings
rampOptions = struct('speeds',5:2.5:30,'nRamp',12,'mode','coast', ...
    'nBisect',6,'verbose',false);
numWorkers = 12;                    % 0 = serial; use one worker per car otherwise
loadFromPrev = false;              % true = redraw from cache without re-solving
cacheName = "ramp_speed_study.mat";
figurePath = "figures/ramp_speed_study";
figureOptions = struct('formats',{{'png','fig'}},'resolution',200,'stamp',false);

% Choose and order the outputs shown. Delete entries to make a smaller figure.
% Valid names: aero_front_load, aero_rear_load, aero_balance,
% mechanical_balance, handling_balance, front_ride_height, rear_ride_height,
% front_camber, rear_camber, pitch_angle, front_shock_travel,
% rear_shock_travel.
outputs = ["aero_front_load" "aero_rear_load" "aero_balance" ...
    "mechanical_balance" "handling_balance" "front_ride_height" ...
    "rear_ride_height" "front_camber" "rear_camber" "pitch_angle" ...
    "front_shock_travel" "rear_shock_travel"];

%% run / reload
if loadFromPrev
    loaded = load(cacheName,'study'); %#ok<UNRCH> user toggle defaults false
    if ~isfield(loaded,'study') || ~isfield(loaded.study,'results')
        error('runRampSpeedStudy:badCache','%s does not contain a ramp-speed study.',cacheName)
    end
    study = loaded.study;
    results = study.results;
    labels = string(study.labels);
else
    job = simLog.start('ramp speed study');
    results = cell(numel(cars),1);
    if numWorkers > 0
        parfor (i = 1:numel(cars),numWorkers)
            results{i} = rampSweep(cars{i},rampOptions);
        end
    else
        for i = 1:numel(cars)
            fprintf('ramp study: %s (%d of %d)\n',labels(i),i,numel(cars));
            results{i} = rampSweep(cars{i},rampOptions);
        end
    end
    study = struct('results',{results},'labels',labels,'rampOptions',rampOptions, ...
        'numWorkers',numWorkers,'created',datetime('now'));
    save(cacheName,'study','-v7.3');
    simLog.finish(job,'events',{'ramp-speed'},'workers',numWorkers, ...
        'nCases',numel(cars),'car',cars{1}, ...
        'details',sprintf('cache=%s mode=%s speeds=%d', ...
            char(cacheName),rampOptions.mode,numel(rampOptions.speeds)));
end

%% plot
fig = plotRampSpeedStudy(results,labels,struct('outputs',outputs));
saveFigures(figurePath,fig,figureOptions);
