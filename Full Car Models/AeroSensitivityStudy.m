%% AeroSensitivityStudy
% Discrete re-sweep sensitivity of the QSS model to ClA, CdA and CoP.
%
% Each case is a complete independent g-g re-solve (and optionally a full
% events run), so the sensitivities include re-optimization of the operating
% point rather than a linearization about the baseline.
%
% HOW TO USE
% 1) Leave every aeroParams entry in carConfig.m SCALAR. This script builds
%    the sweep itself, one parameter at a time, via aeroStarCases.
% 2) Set the levels below if the defaults are not the range you want.
% 3) Run this script. Results land in S.
%
%
% ClA, CdA and CoP are all swept in the same run, as a star design: one
% baseline plus four perturbations along each axis, 13 cases.
%
% FIGURES
% Everything drawn is written to figPath, named from the figure titles. Set
% it to "" to keep them on screen only, or figOpts.stamp = true to give each
% run its own timestamped subfolder instead of overwriting.


clear
setup_paths

runLapTimes  = true;   % false = g-g envelope only (much faster)
loadFromPrev = false;   % true = reuse a previous solve, plot only

numWorkers   = 16;
loadName = "aero_sweep.mat";   % read when loadFromPrev
saveName = "aero_sweep.mat";   % written whenever something was solved

% Where the figures land. Relative paths are relative to this folder
figPath = "figures/aero_sensitivity";
figOpts = struct('formats',{{'png'}},'resolution',200,'stamp',false);

% multipliers for ClA/CdA, additive offsets for CoP (already a fraction).

levels = struct();
levels.ClA = [0.90 0.95 1.05 1.10];
levels.CdA = [0.90 0.95 1.05 1.10];
levels.CoP = [-0.06 -0.03 0.03 0.06];

[baseCell,eventParams] = carConfig();
if size(baseCell,1) ~= 1
    error(['carConfig produced %d cars -- set every aeroParams entry back ' ...
           'to a scalar. This script does the sweeping now, so a vector ' ...
           'there multiplies the case count on top of it.'], size(baseCell,1));
end

[carCell,plan] = aeroStarCases(baseCell,levels);
numCars = size(carCell,1);
fprintf('aero sensitivity study: %d cases (ClA, CdA, CoP one at a time)\n',numCars);
disp(plan);
tic
job = simLog.start('aero sensitivity study');

%% reload a previous solve

haveGG = false; haveTimes = false; times = []; traces = []; %#ok
if loadFromPrev
    [cached,times,traces,haveGG,haveTimes] = loadAeroCache(loadName,plan,baseCell{1,1});
    if haveGG
        carCell = cached;
        fprintf('  loaded %s: g-g=%d, lap times=%d\n',loadName,haveGG,haveTimes);
    end
end

%% g-g for every case

if ~haveGG
    fprintf('  g-g for %d cases (parallel over cases) ...\n',numCars);
    gg = cell(numCars,1);
    parfor (i = 1:numCars, numWorkers)
        gg{i} = makeGG(gg2(carCell{i,1},0),carCell{i,1});
    end
    for i = 1:numCars, carCell{i,1} = gg{i}; end
    fprintf('  g-g done (%.0f s elapsed)\n',toc);
end

%% lap times
opts = struct();
if runLapTimes
    if ~haveTimes
        fprintf('  events for %d cases ...\n',numCars);
        [times,traces] = aeroLapTimes(carCell,numWorkers,[],eventParams);
        fprintf('  events done (%.0f s elapsed)\n',toc);
    end
    opts.times     = times;
    opts.traces    = traces;
    opts.fastestBy = 'autocross';   % which event defines "fastest"
end

%% cache

if ~haveGG || (runLapTimes && ~haveTimes)
    fprintf('  saving %s ...\n',saveName);
    save(saveName,'carCell','plan','-v7.3');
    if runLapTimes
        save(saveName,'times','traces','-append');
    end
end


%% sensitivities
S = aeroSensitivity(carCell,opts);
fprintf('done in %.0f s\n',toc);

% one row for the whole study; the g-g and events flags say how much was
% reused from loadName versus solved fresh
ev = 'g-g';
if runLapTimes, ev = 'g-g+skidpad+accel+autocross+endurance'; end
simLog.finish(job,'events',ev,'workers',numWorkers,'nCases',numCars, ...
    'car',baseCell{1,1}, ...
    'details',sprintf('loadedGG=%d loadedTimes=%d cache=%s',haveGG,haveTimes,char(loadName)));

%% report
aeroSensitivityReport(S);

%% plots

swept = fieldnames(S.sens);

if isfield(S,'traces') && ~isempty(S.traces)

    close all
    for k = 1:numel(swept)
        plotAeroKernel(S,'autocross',swept{k});
        plotAeroOverlay(S,swept{k},'autocross',12,'seconds');
    end
    saveFigures(figPath,[],figOpts);
else
    warning('AeroSensitivityStudy:noPlots', ...
        'both figures need lap traces -- set runLapTimes = true');
end
