function envelope = buildContinuousEnvelope(car,speed_mps)
%BUILDCONTINUOUSENVELOPE Compile the best bounded continuous powertrain ratio.
%
% The envelope is deliberately a small, immutable numeric record. It is
% compiled once per requested vehicle speed and can then be passed to
% Car.equations through the continuousRatio evaluation option.

if ~isobject(car) || ~isscalar(car) || ~isprop(car,'powertrain') || ...
        ~isprop(car,'R')
    error('rampSpeed:invalidCar', ...
        'A scalar Car with powertrain and wheel-radius data is required.');
end
if ~isnumeric(speed_mps) || ~isreal(speed_mps) || ~isscalar(speed_mps) || ...
        ~isfinite(speed_mps) || speed_mps <= 0
    error('rampSpeed:invalidSpeed', ...
        'speed_mps must be a finite positive scalar.');
end

powertrain = car.powertrain;
if isempty(powertrain) || ~isscalar(powertrain)
    error('rampSpeed:invalidPowertrain', ...
        'The Car must contain one scalar powertrain.');
end

gearRatios = double(powertrain.gears(:));
if isempty(gearRatios) || any(~isfinite(gearRatios)) || any(gearRatios <= 0)
    error('rampSpeed:invalidPowertrain', ...
        'Powertrain gear ratios must be finite and positive.');
end
if ~isscalar(powertrain.final_drive) || ~isfinite(powertrain.final_drive) || ...
        powertrain.final_drive <= 0 || ~isscalar(powertrain.primary_reduction) || ...
        ~isfinite(powertrain.primary_reduction) || powertrain.primary_reduction <= 0
    error('rampSpeed:invalidPowertrain', ...
        'Powertrain reductions must be finite and positive scalars.');
end
if ~isscalar(powertrain.wheel_radius) || ...
        ~isfinite(powertrain.wheel_radius) || powertrain.wheel_radius <= 0
    error('rampSpeed:invalidPowertrain', ...
        'Powertrain wheel_radius must be a finite positive scalar.');
end
if ~isscalar(powertrain.drivetrain_efficiency) || ...
        ~isfinite(powertrain.drivetrain_efficiency) || ...
        powertrain.drivetrain_efficiency <= 0
    error('rampSpeed:invalidPowertrain', ...
        'Powertrain drivetrain_efficiency must be a finite positive scalar.');
end

torqueMap = double(powertrain.torque_fn);
if size(torqueMap,1) ~= 2 || size(torqueMap,2) < 2 || ...
        any(~isfinite(torqueMap(:)))
    error('rampSpeed:invalidTorqueMap', ...
        'Powertrain torque_fn must be a finite 2-by-N map with N >= 2.');
end
[rpm,order] = sort(torqueMap(1,:));
torqueFtLb = torqueMap(2,order);
[rpm,uniqueIndex] = unique(rpm,'stable');
torqueFtLb = torqueFtLb(uniqueIndex);
rpm = rpm(:);
torqueFtLb = torqueFtLb(:);
if numel(rpm) < 2 || any(diff(rpm) <= 0)
    error('rampSpeed:invalidTorqueMap', ...
        'Powertrain torque_fn RPM values must contain two unique points.');
end

redline = double(powertrain.redline);
if ~isscalar(redline) || ~isfinite(redline) || redline <= 0
    error('rampSpeed:invalidPowertrain', ...
        'Powertrain redline must be a finite positive scalar.');
end

totalReductions = gearRatios*double(powertrain.final_drive)* ...
    double(powertrain.primary_reduction);
rpmPerReduction = double(speed_mps)/double(car.R)*30/pi;
rpmLower = max(rpm(1),0);
rpmUpper = min(rpm(end),redline);
ratioLower = max(min(totalReductions),rpmLower/rpmPerReduction);
ratioUpper = min(max(totalReductions),rpmUpper/rpmPerReduction);
if ~isfinite(ratioLower) || ~isfinite(ratioUpper) || ratioLower > ratioUpper
    error('rampSpeed:invalidRatioRange', ...
        'No continuous powertrain ratio is feasible at %.6g m/s.',speed_mps);
end

% On each torque-map segment, full-throttle wheel force is a quadratic in
% engine RPM. Checking segment endpoints and the concave vertex gives the
% exact maximum for the piecewise-linear torque map without a nested solver.
rpmCandidates = [rpmLower; rpmUpper];
internal = rpm(rpm > rpmLower & rpm < rpmUpper);
segmentEdges = [rpmLower; internal; rpmUpper];
rpmCandidates = [rpmCandidates; internal];
for i = 1:numel(segmentEdges)-1
    a = segmentEdges(i);
    b = segmentEdges(i+1);
    if b <= a
        continue
    end
    mid = (a+b)/2;
    segment = find(rpm <= mid,1,'last');
    segment = min(segment,numel(rpm)-1);
    slope = (torqueFtLb(segment+1)-torqueFtLb(segment))/ ...
        (rpm(segment+1)-rpm(segment));
    intercept = torqueFtLb(segment)-slope*rpm(segment);
    if slope < 0
        vertex = -intercept/(2*slope);
        if vertex > a && vertex < b
            rpmCandidates(end+1,1) = vertex; %#ok<AGROW>
        end
    end
end

% Convert the RPM candidates back into the bounded continuous-ratio domain.
ratioCandidates = min(max(rpmCandidates/rpmPerReduction,ratioLower),ratioUpper);
candidateRpm = ratioCandidates*rpmPerReduction;
candidateTorqueFtLb = interp1(rpm,torqueFtLb,candidateRpm,'linear');
candidateTorqueNm = candidateTorqueFtLb*1.35581795;
candidateWheelForce = candidateTorqueNm*double(powertrain.drivetrain_efficiency).* ...
    ratioCandidates/double(powertrain.wheel_radius);
[wheelForce,best] = max(candidateWheelForce);
if isempty(best) || ~isfinite(wheelForce) || wheelForce <= 0
    error('rampSpeed:invalidTorqueMap', ...
        'The torque map has no positive continuous-envelope force at %.6g m/s.', ...
        speed_mps);
end

envelope = struct( ...
    'powertrainModel',"continuousEnvelope", ...
    'speed_mps',double(speed_mps), ...
    'drivetrainReduction',double(ratioCandidates(best)), ...
    'ratioLowerBound',double(ratioLower), ...
    'ratioUpperBound',double(ratioUpper), ...
    'rpmLowerBound',double(rpmLower), ...
    'rpmUpperBound',double(rpmUpper), ...
    'engineRpm',double(candidateRpm(best)), ...
    'engineTorque_Nm',double(candidateTorqueNm(best)), ...
    'wheelForce_N',double(wheelForce), ...
    'exitflag',1, ...
    'candidateCount',numel(ratioCandidates));
end
