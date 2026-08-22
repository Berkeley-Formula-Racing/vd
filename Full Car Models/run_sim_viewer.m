%% run_sim_viewer
% Opens the channel viewer on a solved lap.
%

clear
setup_paths

%% ---------------- what to look at ----------------
carSource = "aero_sweep.mat";
caseIdx   = 1;
event     = "autocross";

cacheName    = "sim_viewer.mat";
loadFromPrev = true;

%%

if ~isfile(carSource)
    error('run_sim_viewer:noCar', ...
        ['%s not found. Run AeroSensitivityStudy or run_grip_sweep first, ' ...
         'or point carSource at a .mat holding a solved carCell.'],carSource);
end
L = load(carSource,'carCell');
if ~isfield(L,'carCell') || size(L.carCell,1) < caseIdx
    error('run_sim_viewer:badCache','%s has no carCell case %d',carSource,caseIdx);
end
car = L.carCell{caseIdx,1};
acc = L.carCell{caseIdx,2};
if isempty(car.ss_info)
    error('run_sim_viewer:unsolved', ...
        'case %d in %s never went through gg2/makeGG',caseIdx,carSource);
end

%% ---------------- build or reload ----------------

ev = [];
if loadFromPrev && isfile(cacheName)
    P = load(cacheName);
    if isfield(P,'ev') && isa(P.ev,'Events2') && ...
            isequal(P.ev.car,car) && isequal(P.ev.accelCar,acc)
        ev = P.ev;
        fprintf('loaded %s\n',cacheName);
    else
        warning('run_sim_viewer:staleCache', ...
            ['%s was built from a different car than %s case %d holds now -- ' ...
             'rebuilding.'],cacheName,carSource,caseIdx);
    end
end

if isempty(ev)
    [~,eventParams] = carConfig();
    fprintf('building Events2 (~35 s) ...\n'); t = tic;
    ev = Events2(car,acc,eventParams);
    fprintf('  %.0f s\n',toc(t));
end

%%
if isempty(ev.(event))
    solver = char(event); solver(1) = upper(solver(1));
    fprintf('solving %s ...\n',event); t = tic;
    ev.(solver)();
    fprintf('  %.0f s, %s = %.4f s\n',toc(t),event,ev.times.(char(event)));
end

% provenance only now -- the reload gate above is isequal on the car itself,
% so a stale key can no longer let an old lap through
key = {char(carSource),caseIdx};
save(cacheName,'ev','car','key');

%% ---------------- open ----------------
% third argument is the car, which unlocks the envelope-derived channels
% (grip used, longitudinal power) on top of the solved and algebraic ones
simViewer(ev,char(event),car);

%% ---------------- using it ----------------
% channels    ctrl-click in the list to add or remove. A leading * marks a
%             channel the solver produced directly; everything else is exact
%             algebra on those, except "grip used" which also reads the g-g.
% lap position  drag the slider, or click anywhere on a trace. The dot on the
%             track map follows, and the readout shows every selected channel
%             at that point.
% colour map by   recolours the track by any channel -- speed to see the fast
%             sections, "grip used" to see where the lap is on the limit.
% x: distance / time   switches the axis; the slider range follows.
% Export CSV  writes the selected channels against distance and time.
