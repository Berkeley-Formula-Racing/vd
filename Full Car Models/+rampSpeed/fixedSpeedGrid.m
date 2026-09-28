function [speeds,meta] = fixedSpeedGrid(range_mps,preset)
%FIXEDSPEEDGRID Create a deterministic fixed grid for a named preset.

if nargin < 2 || isempty(preset)
    preset = "accurate";
end
range_mps = double(range_mps(:));
if numel(range_mps) ~= 2 || any(~isfinite(range_mps)) || ...
        any(range_mps <= 0) || range_mps(2) < range_mps(1)
    error("rampSpeed:invalidFixedSpeedRange", ...
        "Fixed speed range must contain two positive finite endpoints.");
end

requested = lower(strtrim(string(preset)));
if requested == "highaccuracy"
    requested = "highaccuracy";
end
presets = rampSpeed.fixedSpeedPresets();
ids = lower(string({presets.id}));
index = find(ids == requested,1);
if isempty(index)
    error("rampSpeed:invalidFixedSpeedPreset", ...
        "Unknown fixed speed preset ""%s"".",string(preset));
end

pointCount = presets(index).pointCount;
if range_mps(2) == range_mps(1)
    speeds = range_mps(1);
else
    speeds = linspace(range_mps(1),range_mps(2),pointCount).';
end
speeds = double(speeds(:));
meta = struct("mode","fixed", ...
    "preset",string(presets(index).id), ...
    "label",string(presets(index).label), ...
    "pointCount",pointCount, ...
    "range_mps",range_mps);
end
