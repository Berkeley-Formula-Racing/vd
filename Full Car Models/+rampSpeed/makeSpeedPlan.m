function plan = makeSpeedPlan(request)
%MAKESPEEDPLAN Create one stable planned task for each normalized speed.
request = rampSpeed.normalizeRequest(request);
speeds = request.settings.speeds_mps;
speedCount = numel(speeds);

if request.settings.speedPolicy == "adaptive"
    origin = "seed";
else
    origin = "requested";
end

tasks = table((1:speedCount)', speeds, repmat(origin, speedCount, 1), ...
    ones(speedCount, 1), repmat("planned", speedCount, 1), ...
    'VariableNames', {'speedIndex', 'speed_mps', 'origin', 'passIndex', 'status'});
plan = struct('tasks', tasks);
end
