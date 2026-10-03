function P = rampSweepState(z,v,ayTarget,holdSpeed)
%RAMPSWEEPSTATE Reconstruct the full steady-state vector from ramp variables.
%   The reduced state supports coast [steer, latVel] and balanced
%   [steer, latVel, throttle, rearSlip] reconstruction.  The runtime keeps
%   balanced optimization in the historical full solver because its two rear
%   wheel residuals remain independent after the exact identities are removed.
%   The remaining entries are prescribed by the ramp kinematics and the
%   original rampSweep bounds.

if nargin < 4 || ~isscalar(holdSpeed)
    error('rampSweepState:badMode','holdSpeed must be a scalar logical flag.');
end
holdSpeed = logical(holdSpeed);
z = z(:).';
expected = 2 + 2*holdSpeed;
if numel(z) ~= expected
    error('rampSweepState:badState', ...
        'Expected %d reduced variables, got %d.',expected,numel(z));
end
if ~isscalar(v) || ~isscalar(ayTarget) || ~isfinite(v) || v == 0 || ...
        ~isfinite(ayTarget)
    error('rampSweepState:badKinematics', ...
        'v and ayTarget must be finite scalar values with nonzero v.');
end

if holdSpeed
    throttle = z(3);
    rearSlip = z(4);
else
    throttle = 0;
    rearSlip = 0;
end

P = [z(1),throttle,v,z(2),ayTarget/v,0,0,rearSlip,rearSlip];
end
