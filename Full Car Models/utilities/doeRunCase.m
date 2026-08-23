function result = doeRunCase(car,accelCar,eventParams,study,caseIndex,runMode)
%DOERUNCASE Run one DOE vehicle case, including optional ramp metrics.

if nargin < 6 || isempty(runMode), runMode = "full"; end
runMode = lower(string(runMode));
startTime = tic;
result = emptyResult(car,accelCar,caseIndex);

try
    if ~isscalar(runMode) || ~any(runMode == ["full","ramponly"])
        error('doeRunCase:badRunMode', ...
            'runMode must be "full" or "rampOnly".')
    end

    if runMode == "full"
        rawGG = gg2(car,0);
        car = makeGG(rawGG,car);
        comp = Events2(car,accelCar,eventParams);
        events = string(study.events);
        if any(events == "skidpad"), comp.Skidpad(); end
        if any(events == "accel"), comp.Accel(); end
        if any(events == "autocross"), comp.Autocross(); end
        if any(events == "endurance"), comp.Endurance(); end
        if all(ismember(["skidpad","accel","autocross","endurance"],events))
            points = comp.computePoints();
        else
            points = emptyPoints();
        end
        car.comp = comp;
    else
        if isempty(car) || isempty(car.ggPoints) || isempty(car.comp)
            error('doeRunCase:unsolvedCar', ...
                'rampOnly requires solved g-g data and car.comp.')
        end
        points = cachedPoints(car.comp);
    end

    rampResult = [];
    if study.ramps.enabled
        rampOptions = study.ramps;
        rampOptions.verbose = false;
        rampResult = rampSweep(car,rampOptions);
        result.rampSummary = rampResult.perSpeed;
        if getOr(study.ramps,'saveFullPoints',false)
            result.rampPoints = rampResult.points;
        end
    end

    [metricRow,~] = doeCaseMetrics(car,caseIndex,rampResult);
    [~,scoreBreakdown] = doeScoreCase(points,metricRow,study.objective);
    result.car = car;
    result.accelCar = accelCar;
    result.metricRow = appendScore(metricRow,scoreBreakdown);
    result.points = points;
    result.scoreBreakdown = scoreBreakdown;
    result.status = "complete";
catch ME
    [metricRow,~] = doeCaseMetrics([],caseIndex,[]);
    metricRow.valid = false;
    metricRow.error_message = string(ME.message);
    result.metricRow = appendScore(metricRow,result.scoreBreakdown);
    result.status = "failed";
    result.errorIdentifier = string(ME.identifier);
    result.errorMessage = string(ME.message);
end

result.elapsed = toc(startTime);
end

function result = emptyResult(car,accelCar,caseIndex)
[metricRow,~] = doeCaseMetrics([],caseIndex,[]);
scoreBreakdown = emptyScoreBreakdown();
result = struct( ...
    'car',car, ...
    'accelCar',accelCar, ...
    'metricRow',appendScore(metricRow,scoreBreakdown), ...
    'rampSummary',[], ...
    'rampPoints',[], ...
    'points',emptyPoints(), ...
    'scoreBreakdown',scoreBreakdown, ...
    'status',"failed", ...
    'errorIdentifier',"", ...
    'errorMessage',"", ...
    'elapsed',NaN);
end

function points = cachedPoints(comp)
if isempty(comp.points) || ~isstruct(comp.points) || ...
        ~all(isfield(comp.points,{'skidpad','accel','autocross','endurance','total'}))
    error('doeRunCase:missingCachedPoints', ...
        'rampOnly requires cached event points in car.comp.points.')
end
points = comp.points;
end

function points = emptyPoints()
points = struct('skidpad',NaN,'accel',NaN,'autocross',NaN, ...
    'endurance',NaN,'total',NaN);
end

function breakdown = emptyScoreBreakdown()
breakdown = struct( ...
    'points_skidpad',NaN, ...
    'points_accel',NaN, ...
    'points_autocross',NaN, ...
    'points_endurance',NaN, ...
    'modeled_dynamic_points',NaN, ...
    'penalty_invalid_case',NaN, ...
    'penalty_solve_failure',NaN, ...
    'penalty_wheel_lift',NaN, ...
    'penalty_energy',NaN, ...
    'penalty_understeer',NaN, ...
    'penalty_rebalance',NaN, ...
    'total_penalty_points',NaN, ...
    'objective_score',NaN);
end

function row = appendScore(row,scoreBreakdown)
row = [row struct2table(scoreBreakdown)];
end

function value = getOr(s,name,default)
if isfield(s,name), value = s.(name); else, value = default; end
end
