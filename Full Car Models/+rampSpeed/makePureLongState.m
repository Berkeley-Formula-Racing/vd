function state = makePureLongState(car,speed_mps,throttle,rearSlip)
%MAKEPURELONGSTATE Assemble the legacy nine-element pure-longitudinal state.
%
% The model's state/control order is
% [steer, throttle, speed, lateral velocity, yaw rate, kappa_FL,
%  kappa_FR, kappa_RL, kappa_RR]. Only throttle and the common rear slip
% ratio are independent for this ramp mode.

if nargin < 4
    error('rampSpeed:invalidPureLongitudinalState', ...
        'car, speed_mps, throttle, and rearSlip are required.');
end
if ~isobject(car) || ~isprop(car,'R')
    error('rampSpeed:invalidCar', ...
        'makePureLongState requires a solver-ready Car.');
end
validateattributes(speed_mps,{'numeric'},{'real','finite','scalar','positive'}, ...
    mfilename,'speed_mps');
validateattributes(throttle,{'numeric'},{'real','finite','scalar'}, ...
    mfilename,'throttle');
validateattributes(rearSlip,{'numeric'},{'real','finite','scalar'}, ...
    mfilename,'rearSlip');

state = [0,double(throttle),double(speed_mps),0,0,0,0, ...
    double(rearSlip),double(rearSlip)];
end
