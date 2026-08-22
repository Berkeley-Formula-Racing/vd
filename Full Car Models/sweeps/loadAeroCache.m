function [carCell,times,traces,haveGG,haveTimes] = loadAeroCache(loadName,plan,baseCar)
% Reads a cached aero sweep and checks it against the case plan the caller is
% about to run.
%
% The g-g and the lap times are validated independently, so a run that solved
% the g-g and then died in the events reloads the g-g and re-runs only the
% events -- which is the expensive half.
%
% Every rejection is a warning naming the reason, never a silent miss. The
% failure this exists to prevent is plotting a stale sweep against a changed
% configuration and believing the numbers.
%
% inputs  loadName: .mat path written by AeroSensitivityStudy
%         plan:     table from aeroStarCases for the run about to happen
%         baseCar:  optional. carConfig's current baseline car; the cache is
%                   rejected if it was solved for a different one.
% outputs carCell:  the cached cases, or {} if unusable
%         times:    lap time table, or []
%         traces:   lap traces, or []
%         haveGG:   the cached cars carry usable g-g results
%         haveTimes: ...and usable lap times

carCell = {}; times = []; traces = [];
haveGG = false; haveTimes = false;

if ~isfile(loadName)
    warning('loadAeroCache:noCache', ...
        '%s does not exist -- solving from scratch and writing it.',loadName);
    return
end

% into a struct, never straight into the caller: a bare load() would
% overwrite carCell, plan and the run's own flags with whatever the file
% happened to hold, and the run would use the old case set without saying so
prev = load(loadName);

if ~all(isfield(prev,{'carCell','plan'}))
    warning('loadAeroCache:staleCache', ...
        '%s has no carCell/plan -- it predates this format. Re-solving.',loadName);
    return
end

if ~planMatches(prev.plan,plan)
    warning('loadAeroCache:planChanged', ...
        ['%s holds a different case set than the current levels ask for -- ' ...
         're-solving. Point loadName elsewhere if you meant to keep the old ' ...
         'sweep.'],loadName);
    return
end

% A file written before the g-g finished passes every check above while
% carrying unsolved cars, and aeroSensitivity would then fail somewhere far
% from the cause. Confirm there are actually results attached.
if isempty(prev.carCell) || ~isprop(prev.carCell{1,1},'ss_info') || ...
        isempty(prev.carCell{1,1}.ss_info)
    warning('loadAeroCache:noResults', ...
        '%s holds cases that were never solved -- re-solving.',loadName);
    return
end

% The plan pins ClA, CdA and CoP and nothing else. Everything else carConfig
% sets -- mass, grip scaling, final drive -- can change under a cache that
% still passes the plan check, and the study would plot cars solved for a
% different configuration without saying so. Adding grip_scaling_front/rear to
% carConfig moved autocross by 0.84 s, eight times the largest aero effect in
% this sweep, so this is not a quiet failure to leave in place.
if nargin >= 3 && ~isempty(baseCar) && ~sameSetup(prev.carCell{1,1},baseCar)
    warning('loadAeroCache:configChanged', ...
        ['%s was solved for a different carConfig than the one loaded now -- ' ...
         're-solving.'],loadName);
    return
end

carCell = prev.carCell;
haveGG  = true;

if all(isfield(prev,{'times','traces'})) && ~isempty(prev.times) && ...
        height(prev.times) == size(carCell,1)
    times     = prev.times;
    traces    = prev.traces;
    haveTimes = true;
end
end

function ok = planMatches(a,b)
ok = false;
if ~istable(a) || ~istable(b) || ~isequal(size(a),size(b)), return, end
if ~all(ismember({'ClA','CdA','CoP'},a.Properties.VariableNames)), return, end
ok = all(abs(a.ClA-b.ClA) < 1e-9) && ...
     all(abs(a.CdA-b.CdA) < 1e-9) && ...
     all(abs(a.CoP-b.CoP) < 1e-9);
end

function ok = sameSetup(a,b)
% Everything carConfig controls except the aero package -- planMatches already
% covers that -- and the g-g results, which only the cached car has.
ok = false;
skip = {'aero','ggPoints','ss_info','accel_info','decel_info', ...
        'longAccelLookup','longDecelLookup'};
p = setdiff(properties(a),skip);
for k = 1:numel(p)
    if ~isequal(a.(p{k}),b.(p{k})), return, end
end
ok = true;
end
