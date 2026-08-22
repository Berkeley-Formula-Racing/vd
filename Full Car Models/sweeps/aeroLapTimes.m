function [times,traces] = aeroLapTimes(carCell,numWorkers,events,eventParams)
% Runs the dynamic events for every case and returns lap times plus the
% solved velocity/acceleration traces.
%
% Deliberately skips Events2.Understeer (a 41-point constant-radius sweep)
% and computePoints -- neither feeds a lap time, and Understeer roughly
% doubles the cost of an events run.
%
% inputs  carCell:    cell array from carConfig, cars already through makeGG
%         numWorkers: workers for the outer parfor over cases (0 = serial)
%         eventParams: second output of carConfig. Required -- Events2 needs
%                     it and will say so.
%         events:     which events to run (optional). Default all four.
%                     Endurance is ~11 laps and costs more than the other
%                     three together, so a sweep that only needs autocross
%                     should say so. Events not run come back NaN.
% outputs times:      table with one row per case, seconds
%         traces:     n-by-1 cell of structs with .autocross and .endurance,
%                     each holding t (s), v (m/s), gLat (g), gLong (g).
%                     These are what lapKernel turns into time-at-condition
%                     densities.

if nargin < 2, numWorkers = 0; end
% Checked here rather than left to Events2: this runs inside a parfor, where a
% missing variable surfaces as MATLAB:parfor:MissingBroadcastVariable and says
% nothing about what is actually wrong.
if nargin < 4 || isempty(eventParams)
    error('aeroLapTimes:noEventParams', ...
        ['eventParams is required -- Events2 needs it. Take it from ' ...
         'carConfig:\n    [carCell,eventParams] = carConfig();']);
end
if nargin < 3 || isempty(events)
    events = {'skidpad','accel','autocross','endurance'};
end
events = lower(cellstr(events));
bad = events(~ismember(events,{'skidpad','accel','autocross','endurance'}));
if ~isempty(bad)
    error('aeroLapTimes:badEvent','unknown event(s): %s',strjoin(bad,', '));
end
% Endurance is ~11 laps and dominates the run; a sweep that does not need it
% should not pay for it. Booleans rather than repeated ismember so the parfor
% body stays transparent to the slicing analysis.
doSkid = ismember('skidpad',events);   doAcc  = ismember('accel',events);
doAuto = ismember('autocross',events); doEnd  = ismember('endurance',events);

n = size(carCell,1);
skidpad = nan(n,1); accel = nan(n,1);
autocross = nan(n,1); endurance = nan(n,1);
traces = cell(n,1);

% Split into single-column cells so parfor can SLICE them. Indexed as
% carCell{i,1} the two-subscript form defeats the slicing analysis, and the
% whole array -- every case's g-g lookup tables, tens of MB -- is broadcast
% to every worker instead of one case going to the worker that needs it.
cars   = carCell(:,1);
accels = carCell(:,2);

parfor (i = 1:n, numWorkers)
    car      = cars{i};
    accelCar = accels{i};
    tr = struct('autocross',[],'endurance',[]);
    try
        comp = Events2(car,accelCar,eventParams);
        % Order is preserved: Events2 accumulates state on the object across
        % calls, so running a subset is not always identical to picking rows
        % out of a full run.
        if doSkid, comp.Skidpad();   skidpad(i)   = comp.times.skidpad;   end
        if doAcc,  comp.Accel();     accel(i)     = comp.times.accel;     end
        if doAuto, comp.Autocross(); autocross(i) = comp.times.autocross; end
        if doEnd,  comp.Endurance(); endurance(i) = comp.times.endurance; end

        if doAuto, tr.autocross = grabTrace(comp.autocross,car.g); end
        if doEnd,  tr.endurance = grabTrace(comp.endurance,car.g); end
    catch ME
        % a diverged case must not take the whole sweep down; it shows up as
        % NaN in the sensitivity and is reported by aeroSensitivityReport
        warning('aeroLapTimes:caseFailed','case %d failed: %s',i,ME.message);
    end
    traces{i} = tr;
end

times = table((1:n)',skidpad,accel,autocross,endurance, ...
    'VariableNames',{'case','skidpad','accel','autocross','endurance'});
end

function s = grabTrace(ev,g)
% pull the solved profile out of an Events2 result into a flat struct
s = struct('t',[],'v',[],'gLat',[],'gLong',[]);
if isempty(ev) || ~isfield(ev,'long_vel'), return, end
v = ev.long_vel(:);
t = ev.time_vec(:);
% time_vec and the profiles can differ by one sample depending on how the
% solver terminates; trim to the shorter of the two rather than guessing
m = min(numel(v),numel(t));
s.t = t(1:m);
s.v = v(1:m);
if isfield(ev,'lat_accel') && numel(ev.lat_accel) >= m
    s.gLat = ev.lat_accel(1:m)/g;  s.gLat = s.gLat(:);
end
if isfield(ev,'long_accel') && numel(ev.long_accel) >= m
    s.gLong = ev.long_accel(1:m)/g; s.gLong = s.gLong(:);
end
end
