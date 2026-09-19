function run = makeRun(type,mode,settings,caseInfo)
%MAKERUN Create a typed normalized run with canonical SI columns.

if nargin < 1 || isempty(type)
    type = "lateral";
end
if nargin < 2 || isempty(mode)
    mode = "";
end
if nargin < 3 || isempty(settings)
    settings = struct();
end
if nargin < 4 || isempty(caseInfo)
    caseInfo = struct();
end
type = lower(string(type));
if ~isscalar(type) || ~any(type == ["lateral","longitudinal"])
    error('rampSpeed:unsupportedType', ...
        'type must be "lateral" or "longitudinal".');
end
if ~isstruct(settings) || ~isscalar(settings)
    error('rampSpeed:invalidSettings','settings must be a scalar struct.');
end
if ~isstruct(caseInfo) || ~isscalar(caseInfo)
    error('rampSpeed:invalidCaseInfo','caseInfo must be a scalar struct.');
end

if type == "lateral" && strlength(string(mode)) == 0
    mode = "coast";
end
mode = lower(string(mode));

speeds = requestedSpeeds(settings);
run = struct();
run.schemaVersion = 1;
run.caseId = getString(caseInfo,'id',getString(caseInfo,'caseId',""));
run.type = type;
run.mode = mode;
run.settings = settings;
run.perSpeed = emptyPerSpeedTable(speeds,type);
run.points = emptyPointTable();
run.runMeta = defaultRunMeta(caseInfo,type,speeds);
run.status = "pending";
run.raw = struct();
end

function speeds = requestedSpeeds(settings)
if isfield(settings,'speeds') && ~isempty(settings.speeds)
    speeds = double(settings.speeds(:));
else
    speeds = zeros(0,1);
end
end

function value = getString(s,name,default)
if isfield(s,name) && ~isempty(s.(name))
    value = string(s.(name));
else
    value = string(default);
end
value = value(1);
end

function meta = defaultRunMeta(caseInfo,type,speeds)
meta = struct();
meta.caseInfo = caseInfo;
meta.lateralMetricsApplicable = type == "lateral";
meta.requestedSpeeds_mps = speeds;
meta.created = datetime('now');
meta.started = datetime.empty;
meta.completed = datetime.empty;
meta.source = "rampSpeed.makeRun";
meta.warnings = strings(0,1);
meta.errors = strings(0,1);
meta.provenance = struct();
meta.solver = struct();
end

function T = emptyPerSpeedTable(speeds,type)
names = perSpeedNames();
types = perSpeedTypes();
T = typedTable(numel(speeds),names,types);
T.speed_mps = speeds;
T.status(:) = "pending";
T.lateral_metrics_applicable(:) = type == "lateral";
end

function names = perSpeedNames()
names = { ...
    'speed_mps','valid','status','reason', ...
    'aLat_mps2','aLong_mps2','engine_rpm','current_gear','throttle', ...
    'downforce_N','drag_N','ClA_m2','CdA_m2','LoD', ...
    'aero_balance_front','aero_outside_map','aero_residual_m', ...
    'Fz_front_axle_N','Fz_rear_axle_N','min_Fz_N','wheel_lift', ...
    'LLTD','LLT_front_N','LLT_rear_N','long_load_transfer_N', ...
    'aLat_free_mps2','aLat_sustainable_mps2','ramp_complete_fraction', ...
    'truncated','power_limited','K_linear_rad_per_mps2','K_linear_r2', ...
    'K_at_limit_rad_per_mps2','cuo_steer_linear_rad','cuo_steer_limit_rad', ...
    'mechanical_balance_front','grip_balance_mid','grip_balance_limit', ...
    'alpha_balance_mid_rad','alpha_balance_limit_rad', ...
    'LLT_norm_balance_mid','LLT_norm_balance_limit', ...
    'front_Fz_fraction_mid','front_Fz_fraction_limit', ...
    'front_downforce_N','rear_downforce_N','front_ride_height_m', ...
    'rear_ride_height_m','pitch_rad','front_shock_travel_m', ...
    'rear_shock_travel_m','front_camber_rad','rear_camber_rad', ...
    'min_Fz_limit_N','max_constraint_residual','n_exitflag1','n_exitflag2', ...
    'aLong_max_mps2','aLat_achieved_mps2','aLat_force_residual_mps2', ...
    'pure_ay0','steer_zero','lat_velocity_zero','yaw_rate_zero', ...
    'rear_slip_ratio','throttle_upper_active','rear_slip_upper_active', ...
    'traction_limited','lateral_metrics_applicable'};
end

function types = perSpeedTypes()
names = perSpeedNames();
types = repmat({'double'},1,numel(names));
logicalNames = {'valid','aero_outside_map','wheel_lift','truncated', ...
    'power_limited','pure_ay0','steer_zero','lat_velocity_zero', ...
    'yaw_rate_zero','throttle_upper_active','rear_slip_upper_active', ...
    'traction_limited','lateral_metrics_applicable'};
for i = 1:numel(logicalNames)
    types{strcmp(names,logicalNames{i})} = 'logical';
end
types{strcmp(names,'status')} = 'string';
types{strcmp(names,'reason')} = 'string';
end

function names = pointNames()
names = { ...
    'speed_mps','speed_index','point_index','valid','status', ...
    'exitflag','max_constraint_residual','max_equality_residual', ...
    'max_inequality_violation','aLat_mps2','aLong_mps2','steer_rad', ...
    'lat_velocity_mps','yaw_rate_rps','engine_rpm','current_gear','throttle', ...
    'downforce_N','drag_N','ClA_m2','CdA_m2','LoD','aero_balance_front', ...
    'aero_outside_map','aero_residual_m','Fz_front_axle_N','Fz_rear_axle_N', ...
    'min_Fz_N','wheel_lift','LLTD','LLT_front_N','LLT_rear_N', ...
    'long_load_transfer_N','aLong_max_mps2','aLat_achieved_mps2', ...
    'aLat_force_residual_mps2','pure_ay0','steer_zero','lat_velocity_zero', ...
    'yaw_rate_zero','lateral_metrics_applicable', ...
    'Fz_FL_N','Fz_FR_N','Fz_RL_N','Fz_RR_N', ...
    'Fx_FL_N','Fx_FR_N','Fx_RL_N','Fx_RR_N', ...
    'Fy_FL_N','Fy_FR_N','Fy_RL_N','Fy_RR_N', ...
    'alpha_FL_rad','alpha_FR_rad','alpha_RL_rad','alpha_RR_rad', ...
    'gamma_FL_rad','gamma_FR_rad','gamma_RL_rad','gamma_RR_rad', ...
    'kappa_FL','kappa_FR','kappa_RL','kappa_RR', ...
    'T_FL_Nm','T_FR_Nm','T_RL_Nm','T_RR_Nm', ...
    'omega_FL_rps','omega_FR_rps','omega_RL_rps','omega_RR_rps'};
end

function types = pointTypes()
names = pointNames();
types = repmat({'double'},1,numel(names));
logicalNames = {'valid','aero_outside_map','wheel_lift','pure_ay0', ...
    'steer_zero','lat_velocity_zero','yaw_rate_zero', ...
    'lateral_metrics_applicable'};
for i = 1:numel(logicalNames)
    types{strcmp(names,logicalNames{i})} = 'logical';
end
types{strcmp(names,'status')} = 'string';
end

function T = emptyPointTable()
T = typedTable(0,pointNames(),pointTypes());
end

function T = typedTable(n,names,types)
columns = cell(1,numel(names));
for i = 1:numel(names)
    switch types{i}
        case 'logical'
            columns{i} = false(n,1);
        case 'string'
            columns{i} = strings(n,1);
        otherwise
            columns{i} = NaN(n,1);
    end
end
T = table(columns{:},'VariableNames',names);
end
