function candidates = candidateGears(car,speed_mps)
%CANDIDATEGEARS Enumerate explicit gearbox branches for one vehicle speed.
%
% Every physical gear is returned, including gears that are already above
% redline. The longitudinal solver records those as rejected attempts rather
% than silently dropping a branch at a gear transition.

if ~isobject(car) || ~isprop(car,'powertrain') || isempty(car.powertrain)
    error('rampSpeed:invalidCar', ...
        'candidateGears requires a solver-ready Car with a powertrain.');
end
validateattributes(speed_mps,{'numeric'},{'real','finite','scalar','positive'}, ...
    mfilename,'speed_mps');

powertrain = car.powertrain;
gearCount = numel(powertrain.gears);
if gearCount < 1
    error('rampSpeed:invalidPowertrain', ...
        'The powertrain must define at least one gear.');
end

template = struct('gear',0,'predictedEngineRpm',NaN, ...
    'withinRedline',false,'reason',"");
candidates = repmat(template,gearCount,1);
for gear = 1:gearCount
    reduction = powertrain.drivetrain_reduction(gear);
    rpm = speed_mps/car.R*reduction*30/pi;
    candidates(gear).gear = gear;
    candidates(gear).predictedEngineRpm = rpm;
    candidates(gear).withinRedline = isfinite(rpm) && ...
        rpm <= powertrain.redline + 1e-9;
    if candidates(gear).withinRedline
        candidates(gear).reason = "eligible at zero rear slip";
    else
        candidates(gear).reason = "engine rpm exceeds redline";
    end
end
end
