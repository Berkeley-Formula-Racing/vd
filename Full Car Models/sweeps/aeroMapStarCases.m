function [carCell,plan] = aeroMapStarCases(baseCell,levels)
%AEROMAPSTARCASES Build one-at-a-time ride-height and aeromap corrections.
%   The baseline is the car whose static ride heights match the aeromap
%   reference. Every other case changes exactly one input from that baseline.

if nargin < 2, levels = struct(); end
if size(baseCell,2) < 2 || isempty(baseCell)
    error('aeroMapStarCases:badCell','baseCell must be carConfig''s two-column cell array.');
end

[base,baseAcc] = referenceCase(baseCell);
cfg = base.rideHeightAero;
if ~base.hasRideHeightAero()
    error('aeroMapStarCases:noMap', ...
        'the selected car has no enabled ride-height aeromap. Set use_aeromap = true.');
end

lv = defaultLevels(levels);
names = ["FrontRideHeightIn" "RearRideHeightIn" "ClAScale" "CdAScale" "CoPOffset"];
values = {lv.FrontRideHeightIn,lv.RearRideHeightIn,lv.ClAScale, ...
    lv.CdAScale,lv.CoPOffset};

front = cfg.static_front_ride_height_in;
rear = cfg.static_rear_ride_height_in;
claScale = 1;
cdaScale = 1;
copOffset = 0;
caseId = 1;
sweptParameter = "baseline";
sweptValue = 0;

for p = 1:numel(names)
    for k = 1:numel(values{p})
        value = values{p}(k);
        caseId(end+1,1) = numel(caseId)+1; %#ok<AGROW>
        sweptParameter(end+1,1) = names(p); %#ok<AGROW>
        sweptValue(end+1,1) = value; %#ok<AGROW>
        front(end+1,1) = front(1); %#ok<AGROW>
        rear(end+1,1) = rear(1); %#ok<AGROW>
        claScale(end+1,1) = 1; %#ok<AGROW>
        cdaScale(end+1,1) = 1; %#ok<AGROW>
        copOffset(end+1,1) = 0; %#ok<AGROW>
        switch names(p)
            case "FrontRideHeightIn", front(end) = front(1) + value;
            case "RearRideHeightIn",  rear(end) = rear(1) + value;
            case "ClAScale",          claScale(end) = value;
            case "CdAScale",          cdaScale(end) = value;
            case "CoPOffset",         copOffset(end) = value;
        end
    end
end

plan = table(caseId,sweptParameter,sweptValue,front,rear,claScale,cdaScale,copOffset, ...
    'VariableNames',{'case','swept_parameter','swept_value', ...
    'front_static_ride_height_in','rear_static_ride_height_in', ...
    'cla_scale','cda_scale','cop_offset'});

carCell = cell(height(plan),2);
for i = 1:height(plan)
    car = base;
    car.rideHeightAero.static_front_ride_height_in = plan.front_static_ride_height_in(i);
    car.rideHeightAero.static_rear_ride_height_in = plan.rear_static_ride_height_in(i);
    aero = car.aero;
    aero.map = aero.map.withCorrections(plan.cla_scale(i), ...
        plan.cda_scale(i),plan.cop_offset(i));
    car.aero = aero;
    carCell{i,1} = car;
    carCell{i,2} = baseAcc;
end
end

function [base,baseAcc] = referenceCase(baseCell)
cars = baseCell(:,1);
distance = inf(numel(cars),1);
for i = 1:numel(cars)
    car = cars{i};
    if car.hasRideHeightAero()
        cfg = car.rideHeightAero;
        distance(i) = hypot(cfg.static_front_ride_height_in - ...
            cfg.map_reference_front_ride_height_in, ...
            cfg.static_rear_ride_height_in - cfg.map_reference_rear_ride_height_in);
    end
end
[minimum,index] = min(distance);
if ~isfinite(minimum) || minimum > 1e-10
    error('aeroMapStarCases:noReferenceCase', ...
        ['carConfig must include one car at the aeromap reference ride heights ' ...
         'before this study can construct perturbations.']);
end
base = baseCell{index,1};
baseAcc = baseCell{index,2};
end

function lv = defaultLevels(levels)
lv.FrontRideHeightIn = getOr(levels,'FrontRideHeightIn',[-0.25 -0.125 0.125 0.25]);
lv.RearRideHeightIn = getOr(levels,'RearRideHeightIn',[-0.25 -0.125 0.125 0.25]);
lv.ClAScale = getOr(levels,'ClAScale',[0.90 0.95 1.05 1.10]);
lv.CdAScale = getOr(levels,'CdAScale',[0.90 0.95 1.05 1.10]);
lv.CoPOffset = getOr(levels,'CoPOffset',[-0.06 -0.03 0.03 0.06]);
fields = fieldnames(lv);
for i = 1:numel(fields)
    validateattributes(lv.(fields{i}),{'numeric'},{'vector','real','finite'}, ...
        mfilename,fields{i});
end
if any(lv.ClAScale <= 0) || any(lv.CdAScale <= 0)
    error('aeroMapStarCases:badScale','ClA and CdA scale factors must be positive.');
end
end

function value = getOr(s,name,default)
if isfield(s,name), value = s.(name); else, value = default; end
end
