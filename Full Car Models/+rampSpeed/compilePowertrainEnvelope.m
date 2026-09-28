function [envelopes,failures] = compilePowertrainEnvelope(car,speeds, ...
        allowFailures)
%COMPILEPOWERTRAINENVELOPE Compile the continuous powertrain at each speed.
%
% ENVELOPES is a column struct array of rampSpeed.buildContinuousEnvelope
% records. Duplicate speeds are compiled once, in their first-occurrence
% order. The returned records contain data only and do not retain CAR.
% When ALLOWFAILURES is true, speeds with no feasible ratio are returned in
% FAILURES instead of aborting the complete fixed-speed plan.

if nargin < 3 || isempty(allowFailures)
    allowFailures = false;
end
if ~islogical(allowFailures) && ...
        ~(isnumeric(allowFailures) && isreal(allowFailures))
    error('rampSpeed:invalidFailurePolicy', ...
        'allowFailures must be a scalar logical.');
end
if ~isscalar(allowFailures) || ~isfinite(double(allowFailures))
    error('rampSpeed:invalidFailurePolicy', ...
        'allowFailures must be a scalar logical.');
end
allowFailures = logical(allowFailures);

if ~isnumeric(speeds) || ~isreal(speeds) || ~isvector(speeds) || ...
        isempty(speeds) || any(~isfinite(speeds(:))) || any(speeds(:) <= 0)
    error('rampSpeed:invalidSpeed', ...
        'speeds must be a nonempty vector of finite positive numeric values.');
end

if ~isobject(car) || ~isscalar(car) || ~isprop(car,'powertrain') || ...
        ~isprop(car,'R')
    error('rampSpeed:invalidCar', ...
        'A scalar solver-ready Car with powertrain and wheel-radius data is required.');
end

uniqueSpeeds = unique(double(speeds(:)),'stable');
envelopes = struct([]);
failures = repmat(struct('speed_mps',NaN,'identifier',"",'message',""), ...
    0,1);
for i = 1:numel(uniqueSpeeds)
    speed = uniqueSpeeds(i);
    try
        one = rampSpeed.buildContinuousEnvelope(car,speed);
        if isempty(envelopes)
            envelopes = one;
        else
            envelopes(end+1,1) = one; %#ok<AGROW>
        end
    catch ME
        if ~allowFailures
            rethrow(ME);
        end
        failures(end+1,1) = struct('speed_mps',speed, ...
            'identifier',string(ME.identifier),'message',string(ME.message)); %#ok<AGROW>
    end
end
end
