function T = buildInspectorTable(run,selection,units)
%BUILDINSPECTORTABLE Return one SI/display row for a selected ramp point.

if nargin < 2 || isempty(selection)
    selection = struct();
end
if nargin < 3 || isempty(units)
    units = struct();
end
if ~isstruct(run) || ~isscalar(run)
    error('rampSpeed:invalidRun','run must be a scalar normalized run.');
end

perSpeed = getTable(run,'perSpeed');
points = getTable(run,'points');
speedRowIndex = selectSpeedRow(perSpeed,selection);
pointRowIndex = selectPointRow(points,selection,perSpeed,speedRowIndex);
if isempty(speedRowIndex)
    speedRowIndex = 1;
end

if height(perSpeed) >= speedRowIndex
    speedRow = perSpeed(speedRowIndex,:);
else
    speedRow = perSpeed([],:);
end
if ~isempty(pointRowIndex) && height(points) >= pointRowIndex
    pointRow = points(pointRowIndex,:);
else
    pointRow = points([],:);
end

speedUnit = requestedUnit(units,'speed',"m/s");
forceUnit = requestedUnit(units,'force',"N");
lengthUnit = requestedUnit(units,'length',"m");
angleUnit = requestedUnit(units,'angle',"rad");
accelerationUnit = requestedUnit(units,'acceleration',"m/s^2");
angularRateUnit = requestedUnit(units,'angularRate',"rad/s");

row = struct();
row.setup_id = setupIdentifier(run);
row.case_id = row.setup_id;
row.speed_index = sourceValue(pointRow,speedRow,"speed_index",speedRowIndex);
row.point_index = sourceValue(pointRow,speedRow,"point_index",NaN);
row.speed_mps = sourceValue(pointRow,speedRow,"speed_mps",NaN);
row.speed_display = convertValue(row.speed_mps,"m/s",speedUnit,1);
row.valid = sourceLogical(pointRow,speedRow,"valid",false);
row.status = sourceString(pointRow,speedRow,"status","");
row.reason = sourceString(speedRow,tableEmpty(),"reason","");
row.truncated = sourceLogical(speedRow,tableEmpty(),"truncated",false);
row.power_limited = sourceLogical(speedRow,tableEmpty(),"power_limited",false);
row.wheel_lift = sourceLogical(pointRow,speedRow,"wheel_lift",false);
row.aero_outside_map = sourceLogical(pointRow,speedRow,"aero_outside_map",false);
row.max_constraint_residual = sourceValue(pointRow,speedRow, ...
    "max_constraint_residual",sourceValue(speedRow,tableEmpty(), ...
    "max_constraint_residual",NaN));
row.max_equality_residual = sourceValue(pointRow,speedRow, ...
    "max_equality_residual",NaN);
row.max_inequality_violation = sourceValue(pointRow,speedRow, ...
    "max_inequality_violation",NaN);
row.aero_residual_m = sourceValue(pointRow,speedRow,"aero_residual_m",NaN);
row.aero_status = aeroStatus(row.aero_outside_map,row.aero_residual_m);

forceFields = ["Fz_FL_N","Fz_FR_N","Fz_RL_N","Fz_RR_N", ...
    "Fx_FL_N","Fx_FR_N","Fx_RL_N","Fx_RR_N", ...
    "Fy_FL_N","Fy_FR_N","Fy_RL_N","Fy_RR_N", ...
    "downforce_N","drag_N","front_downforce_N","rear_downforce_N", ...
    "Fz_front_axle_N","Fz_rear_axle_N","min_Fz_N"];
for i = 1:numel(forceFields)
    name = char(forceFields(i));
    value = sourceValue(pointRow,speedRow,name,NaN);
    row.(name) = value;
    row.(displayFieldName(name)) = convertValue(value,"N",forceUnit,1);
end

angleFields = ["alpha_FL_rad","alpha_FR_rad","alpha_RL_rad","alpha_RR_rad", ...
    "gamma_FL_rad","gamma_FR_rad","gamma_RL_rad","gamma_RR_rad", ...
    "steer_rad","pitch_rad"];
for i = 1:numel(angleFields)
    name = char(angleFields(i));
    value = sourceValue(pointRow,speedRow,name,NaN);
    row.(name) = value;
    row.(displayFieldName(name)) = convertValue(value,"rad",angleUnit,1);
end

slipFields = ["kappa_FL","kappa_FR","kappa_RL","kappa_RR"];
for i = 1:numel(slipFields)
    name = char(slipFields(i));
    value = sourceValue(pointRow,speedRow,name,NaN);
    row.(name) = value;
    row.(displayFieldName(name)) = value;
end

rateFields = ["yaw_rate_rps","omega_FL_rps","omega_FR_rps", ...
    "omega_RL_rps","omega_RR_rps"];
for i = 1:numel(rateFields)
    name = char(rateFields(i));
    value = sourceValue(pointRow,speedRow,name,NaN);
    row.(name) = value;
    row.(displayFieldName(name)) = convertValue(value,"rad/s",angularRateUnit,1);
end

accelerationFields = ["aLat_mps2","aLong_mps2","aLat_achieved_mps2", ...
    "aLong_max_mps2","aLat_force_residual_mps2"];
for i = 1:numel(accelerationFields)
    name = char(accelerationFields(i));
    value = sourceValue(pointRow,speedRow,name,NaN);
    row.(name) = value;
    row.(displayFieldName(name)) = convertValue(value,"m/s^2",accelerationUnit,1);
end

lengthFields = ["aero_residual_m","front_ride_height_m","rear_ride_height_m", ...
    "front_shock_travel_m","rear_shock_travel_m"];
for i = 1:numel(lengthFields)
    name = char(lengthFields(i));
    value = sourceValue(pointRow,speedRow,name,NaN);
    row.(name) = value;
    row.(displayFieldName(name)) = convertValue(value,"m",lengthUnit,1);
end

balanceFields = ["aero_balance_front","mechanical_balance_front", ...
    "LLTD","ramp_complete_fraction"];
for i = 1:numel(balanceFields)
    name = char(balanceFields(i));
    value = sourceValue(pointRow,speedRow,name,NaN);
    row.(name) = value;
    row.(displayFieldName(name)) = convertValue(value,"fraction","%",100);
end

T = struct2table(row,'AsArray',true);
end

function index = selectSpeedRow(T,selection)
if isempty(T)
    index = [];
    return
end
n = height(T);
if isfield(selection,'speedIndex')
    requested = double(selection.speedIndex);
    if isfinite(requested) && requested >= 1 && requested <= n && requested == round(requested)
        index = requested;
        return
    end
end
if isfield(selection,'speed_index')
    requested = double(selection.speed_index);
    if isfinite(requested) && requested >= 1 && requested <= n && requested == round(requested)
        index = requested;
        return
    end
end
if isfield(selection,'speed_mps') && ismember('speed_mps',T.Properties.VariableNames)
    target = double(selection.speed_mps);
    [~,index] = min(abs(double(T.speed_mps) - target));
else
    index = min(1,n);
end
end

function index = selectPointRow(T,selection,~,speedRowIndex)
if isempty(T)
    index = [];
    return
end
mask = true(height(T),1);
if isfield(selection,'speedIndex') && ismember('speed_index',T.Properties.VariableNames)
    mask = mask & double(T.speed_index) == double(selection.speedIndex);
elseif isfield(selection,'speed_index') && ismember('speed_index',T.Properties.VariableNames)
    mask = mask & double(T.speed_index) == double(selection.speed_index);
elseif ~isempty(speedRowIndex) && ismember('speed_index',T.Properties.VariableNames)
    mask = mask & double(T.speed_index) == speedRowIndex;
elseif isfield(selection,'speed_mps') && ismember('speed_mps',T.Properties.VariableNames)
    mask = mask & double(T.speed_mps) == double(selection.speed_mps);
end
if isfield(selection,'pointIndex') && ismember('point_index',T.Properties.VariableNames)
    mask = mask & double(T.point_index) == double(selection.pointIndex);
elseif isfield(selection,'point_index') && ismember('point_index',T.Properties.VariableNames)
    mask = mask & double(T.point_index) == double(selection.point_index);
end
candidates = find(mask);
if isempty(candidates)
    index = 1;
else
    index = candidates(1);
end
end

function T = getTable(run,name)
if isfield(run,name) && istable(run.(name))
    T = run.(name);
else
    T = table();
end
end

function T = tableEmpty()
T = table();
end

function value = sourceValue(primary,secondary,name,default)
if ~isempty(primary) && ismember(name,primary.Properties.VariableNames)
    column = primary.(name);
elseif ~isempty(secondary) && ismember(name,secondary.Properties.VariableNames)
    column = secondary.(name);
else
    value = default;
    return
end
if isempty(column)
    value = default;
else
    value = column(1);
end
if isstring(value)
    value = double(str2double(value));
end
value = double(value);
end

function value = sourceLogical(primary,secondary,name,default)
value = sourceValue(primary,secondary,name,double(default));
value = logical(value);
end

function value = sourceString(primary,secondary,name,default)
if ~isempty(primary) && ismember(name,primary.Properties.VariableNames)
    column = primary.(name);
elseif ~isempty(secondary) && ismember(name,secondary.Properties.VariableNames)
    column = secondary.(name);
else
    value = string(default);
    return
end
if isempty(column)
    value = string(default);
else
    value = string(column(1));
end
end

function value = requestedUnit(units,name,default)
if isstruct(units) && isfield(units,name) && ~isempty(units.(name))
    value = string(units.(name));
elseif isa(units,'containers.Map') && isKey(units,name)
    value = string(units(name));
else
    value = string(default);
end
end

function value = convertValue(value,fromUnits,toUnits,scale)
if isempty(value) || ~isfinite(double(value))
    value = double(value);
    return
end
fromUnits = string(fromUnits);
toUnits = string(toUnits);
if fromUnits == "fraction" && toUnits == "%"
    value = double(value) * scale;
elseif fromUnits == "flag" || toUnits == "flag" || fromUnits == toUnits
    value = double(value) * scale;
else
    value = rampSpeed.displayUnits(double(value),fromUnits,toUnits);
end
end

function status = aeroStatus(outsideMap,residual)
if outsideMap
    status = "outside_map";
elseif isfinite(residual)
    status = "within_map";
else
    status = "unknown";
end
end

function id = setupIdentifier(run)
if isfield(run,'caseId') && ~isempty(run.caseId)
    id = string(run.caseId);
elseif isfield(run,'runMeta') && isstruct(run.runMeta) && ...
        isfield(run.runMeta,'caseInfo') && isstruct(run.runMeta.caseInfo) && ...
        isfield(run.runMeta.caseInfo,'id')
    id = string(run.runMeta.caseInfo.id);
else
    id = "";
end
end

function output = displayFieldName(name)
output = regexprep(char(name),'_(N|rad|mps2|rps|m)$','');
output = [output '_display'];
end
