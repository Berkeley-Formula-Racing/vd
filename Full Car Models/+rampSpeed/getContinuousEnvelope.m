function [envelope,found] = getContinuousEnvelope(model,speed_mps,car)
%GETCONTINUOUSENVELOPE Look up one compiled envelope in a ramp model.

envelope = struct();
found = false;
if nargin < 2 || ~isstruct(model) || ~isscalar(model) || ...
        ~isfield(model,'continuousEnvelope') || ...
        isempty(model.continuousEnvelope)
    return
end
if ~isnumeric(speed_mps) || ~isscalar(speed_mps) || ...
        ~isfinite(speed_mps) || speed_mps <= 0
    return
end
if nargin >= 3 && ~isempty(car) && ~matchesCar(model,car)
    return
end

cached = model.continuousEnvelope;
if ~isstruct(cached) || isempty(cached) || ...
        ~isfield(cached,'speed_mps')
    return
end
speeds = double([cached.speed_mps].');
scale = max([1;abs(speeds);abs(double(speed_mps))]);
index = find(abs(speeds-double(speed_mps)) <= 32*eps(scale),1);
if isempty(index)
    return
end
envelope = cached(index);
found = true;
end

function tf = matchesCar(model,car)
tf = false;
if ~isobject(car) || ~isscalar(car) || ~isprop(car,'R') || ...
        ~isprop(car,'powertrain') || ~isstruct(model) || ...
        ~isfield(model,'powertrain') || ~isstruct(model.powertrain) || ...
        ~isfield(model,'vehicle') || ~isstruct(model.vehicle) || ...
        ~isfield(model.vehicle,'wheelRadius_m')
    return
end
powertrain = model.powertrain;
required = {'redline','gears','primary_reduction','final_drive', ...
    'wheel_radius','drivetrain_efficiency','torque_fn'};
if ~all(isfield(powertrain,required))
    return
end
actual = car.powertrain;
for name = required
    field = name{1};
    if ~isprop(actual,field) || ~isequaln(double(powertrain.(field)), ...
            double(actual.(field)))
        return
    end
end
tf = isequaln(double(model.vehicle.wheelRadius_m),double(car.R));
end
