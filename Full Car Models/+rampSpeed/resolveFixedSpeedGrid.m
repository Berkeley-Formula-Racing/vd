function [speeds,meta] = resolveFixedSpeedGrid(requested,grid)
%RESOLVEFIXEDSPEEDGRID Resolve a custom grid or named fixed preset.

if nargin < 1 || isempty(requested)
    requested = zeros(0,1);
end
if nargin < 2 || isempty(grid)
    grid = struct("mode","fixed");
elseif ischar(grid) || (isstring(grid) && isscalar(grid))
    grid = struct("mode",grid);
end
if ~isstruct(grid) || ~isscalar(grid)
    error("rampSpeed:invalidFixedSpeedGrid", ...
        "speedGrid must be a scalar struct or mode.");
end
if ~isfield(grid,"mode") || isempty(grid.mode)
    grid.mode = "fixed";
end

requested = double(requested(:));
range = zeros(0,1);
if isfield(grid,"range_mps") && ~isempty(grid.range_mps)
    range = double(grid.range_mps(:));
end
mode = lower(strtrim(string(grid.mode)));
if mode == "highaccuracy"
    preset = "highAccuracy";
elseif any(mode == ["preview","accurate"])
    preset = char(mode);
elseif mode == "fixed"
    preset = "";
else
    error("rampSpeed:invalidFixedSpeedPreset", ...
        "Speed grid mode ""%s"" is not one of the three fixed presets.", ...
        string(grid.mode));
end

if strlength(string(preset)) > 0
    if isempty(range)
        if isempty(requested)
            range = [5;30];
        elseif numel(requested) == 1
            range = [requested(1);requested(1)];
        else
            range = [min(requested);max(requested)];
        end
    end
    [speeds,baseMeta] = rampSpeed.fixedSpeedGrid(range,preset);
    meta = baseMeta;
else
    if isempty(range)
        if isempty(requested)
            requested = (5:2.5:30).';
        end
        speeds = requested;
    else
        speeds = range;
    end
    speeds = unique(double(speeds(:)),'stable');
    if isempty(speeds)
        error("rampSpeed:invalidSpeedDomain", ...
            "At least one fixed ramp speed is required.");
    end
    meta = struct("mode","fixed","preset","","label","Custom fixed speeds", ...
        "pointCount",numel(speeds), ...
        "range_mps",[speeds(1);speeds(end)]);
end

speeds = double(speeds(:));
if isempty(speeds)
    error("rampSpeed:invalidSpeedDomain", ...
        "At least one fixed ramp speed is required.");
end
if strlength(string(meta.preset)) > 0 && ...
        (any(~isfinite(speeds)) || any(speeds <= 0) || ...
        any(diff(speeds) <= 0))
    error("rampSpeed:invalidSpeedDomain", ...
        "Preset ramp speeds must be finite and strictly increasing.");
end
if strlength(string(meta.preset)) > 0
    source = "fixed-preset";
    reason = "preset_" + lower(string(meta.preset));
else
    source = "requested";
    reason = "user_requested";
end
meta.mode = "fixed";
meta.requestedSpeeds_mps = speeds;
meta.seedSpeeds_mps = speeds;
meta.finalSpeeds_mps = speeds;
meta.passes = 0;
meta.provenance = table(speeds,repmat(source,numel(speeds),1), ...
    zeros(numel(speeds),1),repmat(reason,numel(speeds),1), ...
    'VariableNames',{'speed_mps','source','pass','reason'});
meta.refinementHistory = struct([]);
meta.stopReason = "";
meta.exactRetrySpeeds_mps = zeros(0,1);
end
