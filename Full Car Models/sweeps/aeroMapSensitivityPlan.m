function info = aeroMapSensitivityPlan(plan)
%AEROMAPSENSITIVITYPLAN Validate and describe an aeromap star-design plan.

required = {'case','swept_parameter','swept_value', ...
    'front_static_ride_height_in','rear_static_ride_height_in', ...
    'cla_scale','cda_scale','cop_offset'};
if ~istable(plan) || ~all(ismember(required,plan.Properties.VariableNames))
    error('aeroMapSensitivityPlan:badPlan', ...
        'plan must be the table returned by aeroMapStarCases.');
end

baseline = find(string(plan.swept_parameter) == "baseline");
if numel(baseline) ~= 1
    error('aeroMapSensitivityPlan:badBaseline', ...
        'plan must contain exactly one baseline row.');
end

spec = { ...
    "FrontRideHeightIn", "front_static_ride_height_in", "front static ride height", "in"; ...
    "RearRideHeightIn",  "rear_static_ride_height_in",  "rear static ride height",  "in"; ...
    "ClAScale",          "cla_scale",                    "map ClA scale",            "scale"; ...
    "CdAScale",          "cda_scale",                    "map CdA scale",            "scale"; ...
    "CoPOffset",         "cop_offset",                   "map CoP offset",           "front-balance fraction"};

swept = unique(string(plan.swept_parameter),"stable");
swept(swept == "baseline") = [];
expected = string(spec(:,1));
if ~all(ismember(swept,expected)) || ~all(ismember(expected,swept))
    error('aeroMapSensitivityPlan:unknownParameter', ...
        'plan must contain one star for each aeromap sensitivity input.');
end

parameters = repmat(struct('name',"",'column',"",'label',"",'unit',""),1,size(spec,1));
for i = 1:size(spec,1)
    parameters(i).name = spec{i,1};
    parameters(i).column = spec{i,2};
    parameters(i).label = spec{i,3};
    parameters(i).unit = spec{i,4};
end

info = struct('baseIdx',baseline,'parameters',parameters);
end
