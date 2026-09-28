function run = runLongitudinalRamp(car,settings,caseInfo,callbacks)
%RUNLONGITUDINALRAMP Run pure-longitudinal acceleration on a fixed speed grid.

if nargin < 2, settings = struct(); end
if nargin < 3, caseInfo = struct(); end
if nargin < 4, callbacks = struct(); end
run = rampSpeed.runCanonicalLongitudinalRamp(car,settings,caseInfo,callbacks);
end
