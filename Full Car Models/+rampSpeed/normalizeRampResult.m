function run = normalizeRampResult(raw,type,settings,caseInfo,runMeta)
%NORMALIZERAMPRESULT Map legacy solver output into the versioned SI schema.

if nargin < 2 || isempty(type)
    type = "lateral";
end
if nargin < 3 || isempty(settings)
    settings = struct();
end
if nargin < 4 || isempty(caseInfo)
    caseInfo = struct();
end
if nargin < 5 || isempty(runMeta)
    runMeta = struct();
end

type = lower(string(type));
if ~isscalar(type) || ~any(type == ["lateral","longitudinal"])
    error('rampSpeed:unsupportedType', ...
        'type must be "lateral" or "longitudinal".');
end
if ~isstruct(settings) || ~isscalar(settings)
    error('rampSpeed:invalidSettings','settings must be a scalar struct.');
end

rawPerSpeed = block(raw,'perSpeed');
rawSpeedCount = blockHeight(rawPerSpeed);
rawSpeed = firstNumeric(rawPerSpeed,rawSpeedCount, ...
    {'speed_mps','vCar','speed','long_vel'},1);
requested = requestedSpeeds(settings);
if isempty(requested)
    requested = rawSpeed;
end
if isempty(requested)
    requested = zeros(0,1);
end

[rawPerSpeed,sourcePresent,matchReasons,unusedRawSpeeds,tolerance] = ...
    alignPerSpeed(rawPerSpeed,requested,settings);

settingsForRun = settings;
settingsForRun.speeds = requested(:).';
if type == "lateral"
    defaultMode = "coast";
else
    defaultMode = "";
end
mode = getString(settings,'mode',defaultMode);
if type == "longitudinal"
    mode = "";
end
run = rampSpeed.makeRun(type,mode,settingsForRun,caseInfo);
run.raw = raw;
run.runMeta = mergeStruct(run.runMeta,runMeta);
run.runMeta.source = getString(run.runMeta,'source',"rampSpeed.normalizeRampResult");
run.runMeta.requestedSpeeds_mps = requested(:);
run.runMeta.lateralMetricsApplicable = type == "lateral";
if ~isfield(runMeta,'created')
    run.runMeta.created = datetime.empty;
end
if ~isfield(runMeta,'completed')
    run.runMeta.completed = datetime.empty;
end
run.runMeta.speedMatchTolerance_mps = tolerance;
run.runMeta.speedMatchPolicy = "first unused raw row within absolute tolerance";
run.runMeta.pointGroupingPolicy = ...
    "explicit speed_index when present; otherwise stable repeated speed values";
if ~isempty(unusedRawSpeeds)
    run.runMeta.warnings(end+1,1) = "unused raw speed row(s): " + ...
        join(string(unusedRawSpeeds),", ");
end

T = run.perSpeed;
n = height(T);
rawSpeedCount = n;
T.speed_mps = requested;

T.valid = firstLogical(rawPerSpeed,rawSpeedCount,{'valid'},sourcePresent);
T.valid = fitLogical(T.valid,n) & sourcePresent;
T.status = firstString(rawPerSpeed,rawSpeedCount,{'status'},strings(rawSpeedCount,1));
T.status = fitString(T.status,n);
missingStatus = strlength(strtrim(T.status)) == 0;
T.status(missingStatus & T.valid) = "complete";
T.status(missingStatus & ~T.valid) = "missing";
T.reason = firstString(rawPerSpeed,rawSpeedCount,{'reason'},strings(rawSpeedCount,1));
T.reason = fitString(T.reason,n);
missingReason = strlength(strtrim(T.reason)) == 0 & ~T.valid;
T.reason(missingReason) = "missing source result";
hasMatchReason = strlength(strtrim(matchReasons)) > 0;
T.reason(hasMatchReason) = matchReasons(hasMatchReason);

g = 9.80665;
T.aLat_mps2 = firstNumeric(rawPerSpeed,rawSpeedCount, ...
    {'aLat_mps2','lat_accel'},1);
T.aLat_mps2 = fillMissing(T.aLat_mps2,firstNumeric(rawPerSpeed,rawSpeedCount, ...
    {'gLat','gLat_top'},g));
T.aLong_mps2 = firstNumeric(rawPerSpeed,rawSpeedCount, ...
    {'aLong_mps2','long_accel'},1);
T.aLong_mps2 = fillMissing(T.aLong_mps2,firstNumeric(rawPerSpeed,rawSpeedCount, ...
    {'gLong','gLong_max'},g));
T.engine_rpm = firstNumeric(rawPerSpeed,rawSpeedCount,{'engine_rpm'},1);
T.current_gear = firstNumeric(rawPerSpeed,rawSpeedCount,{'current_gear'},1);
T.throttle = firstNumeric(rawPerSpeed,rawSpeedCount,{'throttle','throttle_top'},1);
T.downforce_N = firstNumeric(rawPerSpeed,rawSpeedCount,{'downforce_N','downforce'},1);
T.drag_N = firstNumeric(rawPerSpeed,rawSpeedCount,{'drag_N','drag'},1);
T.ClA_m2 = firstNumeric(rawPerSpeed,rawSpeedCount,{'ClA_m2','ClA'},1);
T.CdA_m2 = firstNumeric(rawPerSpeed,rawSpeedCount,{'CdA_m2','CdA'},1);
T.LoD = firstNumeric(rawPerSpeed,rawSpeedCount,{'LoD'},1);
T.aero_balance_front = firstNumeric(rawPerSpeed,rawSpeedCount, ...
    {'aero_balance_front','aero_balance','CoP'},1);
T.aero_outside_map = firstLogical(rawPerSpeed,rawSpeedCount, ...
    {'aero_outside_map'},false(rawSpeedCount,1));
T.aero_residual_m = firstNumeric(rawPerSpeed,rawSpeedCount,{'aero_residual_m'},1);
T.aero_residual_m = fillMissing(T.aero_residual_m, ...
    firstNumeric(rawPerSpeed,rawSpeedCount,{'aero_residual_in'},0.0254));
T.Fz_front_axle_N = firstNumeric(rawPerSpeed,rawSpeedCount, ...
    {'Fz_front_axle_N','Fz_front_axle'},1);
T.Fz_rear_axle_N = firstNumeric(rawPerSpeed,rawSpeedCount, ...
    {'Fz_rear_axle_N','Fz_rear_axle'},1);
T.min_Fz_N = firstNumeric(rawPerSpeed,rawSpeedCount,{'min_Fz_N','min_Fz'},1);
T.wheel_lift = firstLogical(rawPerSpeed,rawSpeedCount,{'wheel_lift'}, ...
    finite(T.min_Fz_N) & T.min_Fz_N <= 0);
T.LLTD = firstNumeric(rawPerSpeed,rawSpeedCount,{'LLTD','mech_balance'},1);
T.LLT_front_N = firstNumeric(rawPerSpeed,rawSpeedCount,{'LLT_front_N','LLT_front'},1);
T.LLT_rear_N = firstNumeric(rawPerSpeed,rawSpeedCount,{'LLT_rear_N','LLT_rear'},1);
T.long_load_transfer_N = firstNumeric(rawPerSpeed,rawSpeedCount, ...
    {'long_load_transfer_N','long_load_transfer'},1);

% Lateral metrics are mapped only for a lateral run. The longitudinal schema
% retains the columns for stable table shapes, but marks numeric values NaN.
T.aLat_free_mps2 = firstNumeric(rawPerSpeed,rawSpeedCount,{'aLat_free_mps2'},1);
T.aLat_free_mps2 = fillMissing(T.aLat_free_mps2, ...
    firstNumeric(rawPerSpeed,rawSpeedCount,{'gLat_max'},g));
T.aLat_sustainable_mps2 = firstNumeric(rawPerSpeed,rawSpeedCount, ...
    {'aLat_sustainable_mps2'},1);
T.aLat_sustainable_mps2 = fillMissing(T.aLat_sustainable_mps2, ...
    firstNumeric(rawPerSpeed,rawSpeedCount,{'gLat_top'},g));
T.ramp_complete_fraction = firstNumeric(rawPerSpeed,rawSpeedCount, ...
    {'ramp_complete_fraction','ramp_complete'},1);
T.truncated = firstLogical(rawPerSpeed,rawSpeedCount,{'truncated'}, ...
    finite(T.ramp_complete_fraction) & T.ramp_complete_fraction < 1-1e-12);
T.power_limited = firstLogical(rawPerSpeed,rawSpeedCount,{'power_limited'}, ...
    false(rawSpeedCount,1));
T.K_linear_rad_per_mps2 = firstNumeric(rawPerSpeed,rawSpeedCount, ...
    {'K_linear_rad_per_mps2'},1);
T.K_linear_rad_per_mps2 = fillMissing(T.K_linear_rad_per_mps2, ...
    firstNumeric(rawPerSpeed,rawSpeedCount,{'K_linear'},pi/180/g));
T.K_linear_r2 = firstNumeric(rawPerSpeed,rawSpeedCount,{'K_linear_r2','K_r2'},1);
T.K_at_limit_rad_per_mps2 = firstNumeric(rawPerSpeed,rawSpeedCount, ...
    {'K_at_limit_rad_per_mps2'},1);
T.K_at_limit_rad_per_mps2 = fillMissing(T.K_at_limit_rad_per_mps2, ...
    firstNumeric(rawPerSpeed,rawSpeedCount,{'K_at_limit'},pi/180/g));
T.cuo_steer_linear_rad = firstNumeric(rawPerSpeed,rawSpeedCount, ...
    {'cuo_steer_linear_rad'},1);
T.cuo_steer_linear_rad = fillMissing(T.cuo_steer_linear_rad, ...
    firstNumeric(rawPerSpeed,rawSpeedCount, ...
    {'CUOsteerFromYaw_linear_deg'},pi/180));
T.cuo_steer_limit_rad = firstNumeric(rawPerSpeed,rawSpeedCount, ...
    {'cuo_steer_limit_rad'},1);
T.cuo_steer_limit_rad = fillMissing(T.cuo_steer_limit_rad, ...
    firstNumeric(rawPerSpeed,rawSpeedCount, ...
    {'CUOsteerFromYaw_limit_deg'},pi/180));
T.mechanical_balance_front = firstNumeric(rawPerSpeed,rawSpeedCount, ...
    {'mechanical_balance_front','mech_balance'},1);
T.grip_balance_mid = firstNumeric(rawPerSpeed,rawSpeedCount,{'grip_balance_mid'},1);
T.grip_balance_limit = firstNumeric(rawPerSpeed,rawSpeedCount,{'grip_balance_limit'},1);
T.alpha_balance_mid_rad = firstNumeric(rawPerSpeed,rawSpeedCount, ...
    {'alpha_balance_mid_rad'},1);
T.alpha_balance_mid_rad = fillMissing(T.alpha_balance_mid_rad, ...
    firstNumeric(rawPerSpeed,rawSpeedCount,{'alpha_balance_mid'},pi/180));
T.alpha_balance_limit_rad = firstNumeric(rawPerSpeed,rawSpeedCount, ...
    {'alpha_balance_limit_rad'},1);
T.alpha_balance_limit_rad = fillMissing(T.alpha_balance_limit_rad, ...
    firstNumeric(rawPerSpeed,rawSpeedCount,{'alpha_balance_limit'},pi/180));
T.LLT_norm_balance_mid = firstNumeric(rawPerSpeed,rawSpeedCount, ...
    {'LLT_norm_balance_mid'},1);
T.LLT_norm_balance_limit = firstNumeric(rawPerSpeed,rawSpeedCount, ...
    {'LLT_norm_balance_limit'},1);
T.front_Fz_fraction_mid = firstNumeric(rawPerSpeed,rawSpeedCount, ...
    {'front_Fz_fraction_mid','front_Fz_frac_mid'},1);
T.front_Fz_fraction_limit = firstNumeric(rawPerSpeed,rawSpeedCount, ...
    {'front_Fz_fraction_limit','front_Fz_frac_limit'},1);
T.front_downforce_N = firstNumeric(rawPerSpeed,rawSpeedCount, ...
    {'front_downforce_N','aero_downforce_front_N'},1);
T.rear_downforce_N = firstNumeric(rawPerSpeed,rawSpeedCount, ...
    {'rear_downforce_N','aero_downforce_rear_N'},1);
T.front_ride_height_m = firstNumeric(rawPerSpeed,rawSpeedCount, ...
    {'front_ride_height_m'},1);
T.front_ride_height_m = fillMissing(T.front_ride_height_m, ...
    firstNumeric(rawPerSpeed,rawSpeedCount,{'front_ride_height_in'},0.0254));
T.rear_ride_height_m = firstNumeric(rawPerSpeed,rawSpeedCount, ...
    {'rear_ride_height_m'},1);
T.rear_ride_height_m = fillMissing(T.rear_ride_height_m, ...
    firstNumeric(rawPerSpeed,rawSpeedCount,{'rear_ride_height_in'},0.0254));
T.pitch_rad = firstNumeric(rawPerSpeed,rawSpeedCount,{'pitch_rad'},1);
T.pitch_rad = fillMissing(T.pitch_rad, ...
    firstNumeric(rawPerSpeed,rawSpeedCount,{'pitch_angle_deg'},pi/180));
T.front_shock_travel_m = firstNumeric(rawPerSpeed,rawSpeedCount, ...
    {'front_shock_travel_m'},1);
T.front_shock_travel_m = fillMissing(T.front_shock_travel_m, ...
    firstNumeric(rawPerSpeed,rawSpeedCount,{'front_shock_travel_in'},0.0254));
T.rear_shock_travel_m = firstNumeric(rawPerSpeed,rawSpeedCount, ...
    {'rear_shock_travel_m'},1);
T.rear_shock_travel_m = fillMissing(T.rear_shock_travel_m, ...
    firstNumeric(rawPerSpeed,rawSpeedCount,{'rear_shock_travel_in'},0.0254));
T.front_camber_rad = firstNumeric(rawPerSpeed,rawSpeedCount, ...
    {'front_camber_rad'},1);
T.front_camber_rad = fillMissing(T.front_camber_rad, ...
    firstNumeric(rawPerSpeed,rawSpeedCount,{'front_camber_deg'},pi/180));
T.rear_camber_rad = firstNumeric(rawPerSpeed,rawSpeedCount, ...
    {'rear_camber_rad'},1);
T.rear_camber_rad = fillMissing(T.rear_camber_rad, ...
    firstNumeric(rawPerSpeed,rawSpeedCount,{'rear_camber_deg'},pi/180));
T.min_Fz_limit_N = firstNumeric(rawPerSpeed,rawSpeedCount, ...
    {'min_Fz_limit_N','min_Fz_limit'},1);
T.max_constraint_residual = firstNumeric(rawPerSpeed,rawSpeedCount, ...
    {'max_constraint_residual','max_equality_residual','max_ceq'},1);
T.n_exitflag1 = firstNumeric(rawPerSpeed,rawSpeedCount,{'n_exitflag1','n_exit1'},1);
T.n_exitflag2 = firstNumeric(rawPerSpeed,rawSpeedCount,{'n_exitflag2','n_exit2'},1);

% Longitudinal diagnostics and pure-Ay enforcement flags.
T.aLong_max_mps2 = firstNumeric(rawPerSpeed,rawSpeedCount, ...
    {'aLong_max_mps2','aLong_max','long_accel'},1);
T.aLong_max_mps2 = fillMissing(T.aLong_max_mps2, ...
    firstNumeric(rawPerSpeed,rawSpeedCount,{'gLong_max'},g));
T.aLat_achieved_mps2 = firstNumeric(rawPerSpeed,rawSpeedCount, ...
    {'aLat_achieved_mps2','aLat_achieved','lat_accel'},1);
T.aLat_achieved_mps2 = fillMissing(T.aLat_achieved_mps2, ...
    firstNumeric(rawPerSpeed,rawSpeedCount,{'aLat_mps2'},1));
T.aLat_force_residual_mps2 = firstNumeric(rawPerSpeed,rawSpeedCount, ...
    {'aLat_force_residual_mps2','aLat_force_residual', ...
    'lat_accel_residual_mps2','lat_accel_residual'},1);
zeroTolerance = getNumericSetting(settings, ...
    {'zeroStateTolerance','stateTolerance'},1e-12);
T.pure_ay0 = firstLogical(rawPerSpeed,rawSpeedCount,{'pure_ay0'}, ...
    finite(T.aLat_achieved_mps2) & abs(T.aLat_achieved_mps2) <= zeroTolerance);
steerRad = firstNumeric(rawPerSpeed,rawSpeedCount,{'steer_rad'},1);
steerRad = fillMissing(steerRad,firstNumeric(rawPerSpeed,rawSpeedCount, ...
    {'steer_angle','steer_avg'},pi/180));
latVelocity = firstNumeric(rawPerSpeed,rawSpeedCount, ...
    {'lat_velocity_mps','lat_vel'},1);
yawRate = firstNumeric(rawPerSpeed,rawSpeedCount,{'yaw_rate_rps','yaw_rate'},1);
T.steer_zero = firstLogical(rawPerSpeed,rawSpeedCount,{'steer_zero'}, ...
    finite(steerRad) & abs(steerRad) <= zeroTolerance);
T.lat_velocity_zero = firstLogical(rawPerSpeed,rawSpeedCount, ...
    {'lat_velocity_zero'},finite(latVelocity) & abs(latVelocity) <= zeroTolerance);
T.yaw_rate_zero = firstLogical(rawPerSpeed,rawSpeedCount,{'yaw_rate_zero'}, ...
    finite(yawRate) & abs(yawRate) <= zeroTolerance);
T.rear_slip_ratio = firstNumeric(rawPerSpeed,rawSpeedCount,{'rear_slip_ratio'},1);
rearKappa3 = firstNumeric(rawPerSpeed,rawSpeedCount,{'kappa_3','kappa_RL'},1);
rearKappa4 = firstNumeric(rawPerSpeed,rawSpeedCount,{'kappa_4','kappa_RR'},1);
rearKappa = meanAvailable(rearKappa3,rearKappa4);
T.rear_slip_ratio = fillMissing(T.rear_slip_ratio,rearKappa);
T.throttle_upper_active = firstLogical(rawPerSpeed,rawSpeedCount, ...
    {'throttle_upper_active'},finite(T.throttle) & T.throttle >= 1-1e-9);
T.rear_slip_upper_active = firstLogical(rawPerSpeed,rawSpeedCount, ...
    {'rear_slip_upper_active'},finite(T.rear_slip_ratio) & ...
    T.rear_slip_ratio >= 0.2-1e-9);
T.traction_limited = firstLogical(rawPerSpeed,rawSpeedCount, ...
    {'traction_limited'},T.rear_slip_upper_active);
T.lateral_metrics_applicable(:) = type == "lateral";

perSpeedInequality = firstNumeric(rawPerSpeed,rawSpeedCount, ...
    {'max_inequality_violation','max_constraint_violation'},1);
perSpeedInequality = fillMissing(perSpeedInequality, ...
    boundViolationFromMinLoad(T.min_Fz_N));
[T.valid,T.status,T.reason] = gateValidity(rawPerSpeed,T.valid,T.status, ...
    T.reason,sourcePresent,settings, ...
    firstNumeric(rawPerSpeed,rawSpeedCount,{'exitflag'},1), ...
    T.max_constraint_residual,NaN(n,1),perSpeedInequality, ...
    T.n_exitflag1,T.n_exitflag2,matchReasons);

if type == "longitudinal"
    T = clearLateralMetrics(T);
end
run.perSpeed = T;
run.points = normalizePoints(block(raw,'points'),type,settings);

sourceStatus = firstString(raw,1,{'status'},"");
if isfield(runMeta,'status') && ~isempty(runMeta.status)
    run.status = string(runMeta.status);
elseif strlength(strtrim(sourceStatus)) > 0
    run.status = sourceStatus(1);
elseif n > 0 && any(T.valid)
    run.status = "completed";
else
    run.status = "failed";
end
end

function [valid,status,reason] = gateValidity(data,valid,status,reason, ...
        sourcePresent,settings,exitflag,residual,equalityResidual,inequality, ...
        nExitflag1,nExitflag2,matchReasons)
n = numel(valid);
valid = fitLogical(valid,n);
status = fitString(status,n);
reason = fitString(reason,n);
sourcePresent = fitLogical(sourcePresent,n);
exitflag = fitRows(exitflag,n);
residual = fitRows(residual,n);
equalityResidual = fitRows(equalityResidual,n);
inequality = fitRows(inequality,n);
nExitflag1 = fitRows(nExitflag1,n);
nExitflag2 = fitRows(nExitflag2,n);
matchReasons = fitString(matchReasons,n);

[rawValid,hasValid] = readField(data,'valid',n);
if hasValid
    explicitValid = toLogical(rawValid,n);
else
    explicitValid = false(n,1);
end
rawStatus = firstString(data,n,{'status'},strings(n,1));

residualTolerance = getNumericSetting(settings, ...
    {'ceqTol','residualTolerance','constraintTolerance'},1e-2);
inequalityTolerance = getNumericSetting(settings, ...
    {'inequalityTolerance','constraintTolerance','ceqTol'}, ...
    residualTolerance);

hasExit = finite(exitflag);
[rawNear,hasNear] = readField(data,'accepted_near_feasible',n);
if hasNear
    acceptedNear = toLogical(rawNear,n);
else
    acceptedNear = false(n,1);
end
exitBad = hasExit & ~(exitflag == 1 | exitflag == 2) & ~acceptedNear;
countsKnown = finite(nExitflag1) | finite(nExitflag2);
countGood = (finite(nExitflag1) & nExitflag1 > 0) | ...
    (finite(nExitflag2) & nExitflag2 > 0);
countBad = countsKnown & ~countGood;
residualKnown = finite(residual) | finite(equalityResidual);
residualBad = (finite(residual) & residual > residualTolerance) | ...
    (finite(equalityResidual) & equalityResidual > residualTolerance);
inequalityKnown = finite(inequality);
inequalityBad = inequalityKnown & inequality > inequalityTolerance;
rawStatusLower = lower(strtrim(rawStatus));
statusBad = ismember(rawStatusLower,["failed","invalid","cancelled"]);
canonicalStatus = ismember(rawStatusLower, ...
    ["planned","running","converged","near_feasible","infeasible", ...
    "solver_failed","cancelled","complete","completed"]);
canonicalSuccess = ismember(rawStatusLower, ...
    ["converged","near_feasible","complete","completed"]);
explicitBad = hasValid & ~explicitValid;

bad = sourcePresent & (exitBad | countBad | residualBad | ...
    inequalityBad | statusBad | explicitBad);
evidence = residualKnown | inequalityKnown | countsKnown;
good = sourcePresent & ~bad & evidence;
unknown = sourcePresent & ~bad & ~evidence;
missing = ~sourcePresent;

valid = good;
status(good) = "complete";
status(bad) = "invalid";
status(unknown) = "unknown";
status(missing) = "missing";
status(canonicalStatus & sourcePresent) = rawStatusLower( ...
    canonicalStatus & sourcePresent);
valid(canonicalStatus & sourcePresent) = ...
    canonicalSuccess(canonicalStatus & sourcePresent);

for i = 1:n
    if strlength(strtrim(matchReasons(i))) > 0
        reason(i) = matchReasons(i);
    elseif bad(i) && strlength(strtrim(reason(i))) == 0
        if exitBad(i) || countBad(i)
            reason(i) = "solver exit flag is not feasible";
        elseif residualBad(i)
            reason(i) = "constraint residual exceeds tolerance";
        elseif inequalityBad(i)
            reason(i) = "bound violation exceeds tolerance";
        elseif statusBad(i) || explicitBad(i)
            reason(i) = "source marked invalid";
        else
            reason(i) = "solver feasibility gate rejected row";
        end
    elseif unknown(i) && strlength(strtrim(reason(i))) == 0
        reason(i) = "feasibility evidence unavailable";
    end
    if canonicalStatus(i) && strlength(strtrim(reason(i))) == 0 && ...
            ~canonicalSuccess(i)
        reason(i) = "source status: " + rawStatusLower(i);
    end
end
end

function T = clearLateralMetrics(T)
numericNames = {'aLat_free_mps2','aLat_sustainable_mps2', ...
    'ramp_complete_fraction','K_linear_rad_per_mps2','K_linear_r2', ...
    'K_at_limit_rad_per_mps2','cuo_steer_linear_rad','cuo_steer_limit_rad', ...
    'mechanical_balance_front','grip_balance_mid','grip_balance_limit', ...
    'alpha_balance_mid_rad','alpha_balance_limit_rad', ...
    'LLT_norm_balance_mid','LLT_norm_balance_limit', ...
    'front_Fz_fraction_mid','front_Fz_fraction_limit','min_Fz_limit_N', ...
    'n_exitflag1','n_exitflag2'};
for i = 1:numel(numericNames)
    T.(numericNames{i}) = NaN(height(T),1);
end
T.truncated(:) = false;
end

function P = normalizePoints(blockData,type,settings)
n = blockHeight(blockData);
P = typedTable(n,pointNames(),pointTypes());
if n == 0
    return
end
g = 9.80665;
zeroTolerance = getNumericSetting(settings, ...
    {'zeroStateTolerance','stateTolerance'},1e-12);
P.speed_mps = firstNumeric(blockData,n, ...
    {'speed_mps','vCar','speed','long_vel'},1);
P.speed_index = firstNumeric(blockData,n,{'speed_index','speedIndex'},1);
missingIndex = ~finite(P.speed_index);
if any(missingIndex)
    inferredSpeedIndex = inferStableSpeedIndices(P.speed_mps,settings);
    P.speed_index(missingIndex) = inferredSpeedIndex(missingIndex);
end
P.point_index = firstNumeric(blockData,n,{'point_index','pointIndex'},1);
missingPoint = ~finite(P.point_index);
if any(missingPoint)
    inferredPointIndex = inferWithinSpeedPointIndices(P.speed_index);
    P.point_index(missingPoint) = inferredPointIndex(missingPoint);
end
P.valid = firstLogical(blockData,n,{'valid'},true(n,1));
P.status = firstString(blockData,n,{'status'},strings(n,1));
missingStatus = strlength(strtrim(P.status)) == 0;
P.status(missingStatus & P.valid) = "complete";
P.status(missingStatus & ~P.valid) = "failed";
P.exitflag = firstNumeric(blockData,n,{'exitflag'},1);
P.max_constraint_residual = firstNumeric(blockData,n, ...
    {'max_constraint_residual','max_ceq'},1);
P.max_equality_residual = firstNumeric(blockData,n, ...
    {'max_equality_residual','max_ceq'},1);
P.max_inequality_violation = firstNumeric(blockData,n, ...
    {'max_inequality_violation'},1);
minFz = firstNumeric(blockData,n,{'min_Fz_N','min_Fz'},1);
P.max_inequality_violation = fillMissing(P.max_inequality_violation, ...
    boundViolationFromMinLoad(minFz));
P.aLat_mps2 = firstNumeric(blockData,n,{'aLat_mps2','lat_accel'},1);
P.aLat_mps2 = fillMissing(P.aLat_mps2,firstNumeric(blockData,n,{'gLat'},g));
P.aLong_mps2 = firstNumeric(blockData,n,{'aLong_mps2','long_accel'},1);
P.aLong_mps2 = fillMissing(P.aLong_mps2,firstNumeric(blockData,n,{'gLong'},g));
P.steer_rad = firstNumeric(blockData,n,{'steer_rad'},1);
P.steer_rad = fillMissing(P.steer_rad,firstNumeric(blockData,n, ...
    {'steer_angle','steer_avg'},pi/180));
P.lat_velocity_mps = firstNumeric(blockData,n, ...
    {'lat_velocity_mps','lat_vel'},1);
P.yaw_rate_rps = firstNumeric(blockData,n,{'yaw_rate_rps','yaw_rate'},1);
P.engine_rpm = firstNumeric(blockData,n,{'engine_rpm'},1);
P.current_gear = firstNumeric(blockData,n,{'current_gear'},1);
P.throttle = firstNumeric(blockData,n,{'throttle'},1);
P.downforce_N = firstNumeric(blockData,n,{'downforce_N','downforce'},1);
P.drag_N = firstNumeric(blockData,n,{'drag_N','drag'},1);
P.ClA_m2 = firstNumeric(blockData,n,{'ClA_m2','ClA'},1);
P.CdA_m2 = firstNumeric(blockData,n,{'CdA_m2','CdA'},1);
P.LoD = firstNumeric(blockData,n,{'LoD'},1);
P.aero_balance_front = firstNumeric(blockData,n, ...
    {'aero_balance_front','CoP','aero_balance'},1);
P.aero_outside_map = firstLogical(blockData,n,{'aero_outside_map'},false(n,1));
P.aero_residual_m = firstNumeric(blockData,n,{'aero_residual_m'},1);
P.aero_residual_m = fillMissing(P.aero_residual_m, ...
    firstNumeric(blockData,n,{'aero_residual_in'},0.0254));
P.Fz_front_axle_N = firstNumeric(blockData,n,{'Fz_front_axle_N','Fz_front_axle'},1);
P.Fz_rear_axle_N = firstNumeric(blockData,n,{'Fz_rear_axle_N','Fz_rear_axle'},1);
P.min_Fz_N = minFz;
P.wheel_lift = firstLogical(blockData,n,{'wheel_lift'},finite(minFz) & minFz <= 0);
P.LLTD = firstNumeric(blockData,n,{'LLTD'},1);
P.LLT_front_N = firstNumeric(blockData,n,{'LLT_front_N','LLT_front'},1);
P.LLT_rear_N = firstNumeric(blockData,n,{'LLT_rear_N','LLT_rear'},1);
P.long_load_transfer_N = firstNumeric(blockData,n, ...
    {'long_load_transfer_N','long_load_transfer'},1);
P.aLong_max_mps2 = firstNumeric(blockData,n, ...
    {'aLong_max_mps2','aLong_max','long_accel'},1);
P.aLat_achieved_mps2 = firstNumeric(blockData,n, ...
    {'aLat_achieved_mps2','aLat_achieved','lat_accel'},1);
P.aLat_achieved_mps2 = fillMissing(P.aLat_achieved_mps2,P.aLat_mps2);
P.aLat_force_residual_mps2 = firstNumeric(blockData,n, ...
    {'aLat_force_residual_mps2','aLat_force_residual', ...
    'lat_accel_residual_mps2','lat_accel_residual'},1);
P.pure_ay0 = firstLogical(blockData,n,{'pure_ay0'}, ...
    finite(P.aLat_achieved_mps2) & abs(P.aLat_achieved_mps2) <= zeroTolerance);
P.steer_zero = firstLogical(blockData,n,{'steer_zero'}, ...
    finite(P.steer_rad) & abs(P.steer_rad) <= zeroTolerance);
P.lat_velocity_zero = firstLogical(blockData,n,{'lat_velocity_zero'}, ...
    finite(P.lat_velocity_mps) & abs(P.lat_velocity_mps) <= zeroTolerance);
P.yaw_rate_zero = firstLogical(blockData,n,{'yaw_rate_zero'}, ...
    finite(P.yaw_rate_rps) & abs(P.yaw_rate_rps) <= zeroTolerance);
P.lateral_metrics_applicable(:) = type == "lateral";

corners = {'FL','FR','RL','RR'};
for i = 1:4
    code = corners{i};
    index = num2str(i);
    P.(['Fz_' code '_N']) = firstNumeric(blockData,n, ...
        {['Fz_' code '_N'],['Fz_' code],['Fz_' index]},1);
    P.(['Fx_' code '_N']) = firstNumeric(blockData,n, ...
        {['Fx_' code '_N'],['Fx_' code],['Fx_' index]},1);
    P.(['Fy_' code '_N']) = firstNumeric(blockData,n, ...
        {['Fy_' code '_N'],['Fy_' code],['Fy_' index]},1);
    P.(['alpha_' code '_rad']) = firstNumeric(blockData,n, ...
        {['alpha_' code '_rad'],['alpha_' code],['alpha_' index]},pi/180);
    P.(['gamma_' code '_rad']) = firstNumeric(blockData,n, ...
        {['gamma_' code '_rad'],['gamma_' code],['gamma_' index]},pi/180);
    P.(['kappa_' code]) = firstNumeric(blockData,n, ...
        {['kappa_' code],['kappa_' index]},1);
    P.(['T_' code '_Nm']) = firstNumeric(blockData,n, ...
        {['T_' code '_Nm'],['T_' code],['T_' index]},1);
    P.(['omega_' code '_rps']) = firstNumeric(blockData,n, ...
        {['omega_' code '_rps'],['omega_' code],['omega_' index]},1);
end
[P.valid,P.status] = gateValidity(blockData,P.valid,P.status, ...
    strings(n,1),true(n,1),settings,P.exitflag, ...
    P.max_constraint_residual,P.max_equality_residual, ...
    P.max_inequality_violation, ...
    NaN(n,1),NaN(n,1),strings(n,1));
end

function indices = inferStableSpeedIndices(speeds,settings)
speeds = double(speeds(:));
indices = NaN(size(speeds));
tolerance = getNumericSetting(settings, ...
    {'speedMatchTolerance_mps','speedTolerance_mps'},1e-9);
if ~isfinite(tolerance) || tolerance < 0
    tolerance = 1e-9;
end
uniqueSpeeds = zeros(0,1);
for i = 1:numel(speeds)
    if ~isfinite(speeds(i))
        continue
    end
    match = find(abs(uniqueSpeeds - speeds(i)) <= tolerance,1,"first");
    if isempty(match)
        uniqueSpeeds(end+1,1) = speeds(i); %#ok<AGROW>
        match = numel(uniqueSpeeds);
    end
    indices(i) = match;
end
end

function indices = inferWithinSpeedPointIndices(speedIndices)
speedIndices = double(speedIndices(:));
indices = NaN(size(speedIndices));
groups = unique(speedIndices(isfinite(speedIndices)),"stable");
for i = 1:numel(groups)
    rows = find(speedIndices == groups(i));
    indices(rows) = (1:numel(rows)).';
end
unknownRows = find(~isfinite(indices));
indices(unknownRows) = (1:numel(unknownRows)).';
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

function T = typedTable(n,names,types)
columns = cell(1,numel(names));
for i = 1:numel(names)
    if strcmp(types{i},'logical')
        columns{i} = false(n,1);
    elseif strcmp(types{i},'string')
        columns{i} = strings(n,1);
    else
        columns{i} = NaN(n,1);
    end
end
T = table(columns{:},'VariableNames',names);
end

function S = mergeStruct(S,extra)
if ~isstruct(extra) || ~isscalar(extra)
    return
end
names = fieldnames(extra);
for i = 1:numel(names)
    S.(names{i}) = extra.(names{i});
end
end

function value = getString(s,name,default)
if isstruct(s) && isfield(s,name) && ~isempty(s.(name))
    value = string(s.(name));
else
    value = string(default);
end
value = value(1);
end

function speeds = requestedSpeeds(settings)
if isstruct(settings) && isfield(settings,'speeds') && ~isempty(settings.speeds)
    speeds = double(settings.speeds(:));
else
    speeds = zeros(0,1);
end
end

function value = getNumericSetting(settings,names,default)
value = default;
for i = 1:numel(names)
    if isfield(settings,names{i}) && ~isempty(settings.(names{i}))
        candidate = double(settings.(names{i}));
        if isscalar(candidate)
            value = candidate;
            return
        end
    end
end
end

function data = block(raw,name)
if istable(raw)
    if nargin < 2 || isempty(name) || strcmp(name,'perSpeed')
        data = raw;
    else
        data = [];
    end
elseif isstruct(raw) && isscalar(raw) && isfield(raw,name)
    data = raw.(name);
else
    data = [];
end
end

function [aligned,sourcePresent,matchReasons,unusedRawSpeeds,tolerance] = ...
        alignPerSpeed(data,requested,settings)
rawCount = blockHeight(data);
n = numel(requested);
tolerance = getNumericSetting(settings, ...
    {'speedMatchTolerance_mps','speedTolerance_mps'},1e-9);
if ~isfinite(tolerance) || tolerance < 0
    error('rampSpeed:invalidSpeedTolerance', ...
        'speedMatchTolerance_mps must be a finite non-negative scalar.');
end

rawSpeed = firstNumeric(data,rawCount, ...
    {'speed_mps','vCar','speed','long_vel'},1);
sourceIndex = zeros(n,1);
used = false(rawCount,1);
matchReasons = strings(n,1);
for i = 1:n
    if ~isfinite(requested(i))
        matchReasons(i) = "requested speed is not finite";
        continue
    end
    candidates = find(~used & isfinite(rawSpeed) & ...
        abs(rawSpeed-requested(i)) <= tolerance);
    if isempty(candidates)
        matchReasons(i) = "requested speed not returned by solver";
    else
        sourceIndex(i) = candidates(1);
        used(candidates(1)) = true;
    end
end
sourcePresent = sourceIndex > 0;
unusedRawSpeeds = rawSpeed(~used & isfinite(rawSpeed));

aligned = struct('rampSpeedRowCount',n);
names = blockNames(data);
for i = 1:numel(names)
    [values,found] = readField(data,names{i},rawCount);
    if found
        aligned.(names{i}) = reindexValues(values,sourceIndex,n);
    end
end
end

function names = blockNames(data)
if istable(data)
    names = data.Properties.VariableNames;
elseif isstruct(data)
    names = fieldnames(data).';
    names(strcmp(names,'rampSpeedRowCount')) = [];
else
    names = {};
end
end

function values = reindexValues(values,sourceIndex,n)
if ischar(values)
    values = string(values);
end
if isstring(values)
    aligned = strings(n,1);
elseif islogical(values)
    aligned = false(n,1);
elseif isnumeric(values)
    aligned = NaN(n,1);
elseif isdatetime(values)
    aligned = NaT(n,1);
elseif iscell(values)
    aligned = cell(n,1);
else
    aligned = strings(n,1);
end
target = find(sourceIndex > 0);
if isempty(target)
    values = aligned;
    return
end
if iscell(aligned)
    aligned(target) = values(sourceIndex(target));
else
    aligned(target) = values(sourceIndex(target));
end
values = aligned;
end

function n = blockHeight(data)
if isempty(data)
    n = 0;
elseif istable(data)
    n = height(data);
elseif isstruct(data)
    if isfield(data,'rampSpeedRowCount')
        n = double(data.rampSpeedRowCount);
        return
    end
    if numel(data) > 1
        n = numel(data);
    else
        names = fieldnames(data);
        n = 1;
        for i = 1:numel(names)
            value = data.(names{i});
            if (isnumeric(value) || islogical(value) || isstring(value)) && ...
                    numel(value) > 1
                n = numel(value);
                break
            end
        end
    end
else
    n = numel(data);
end
end

function [value,found] = readField(data,name,n)
found = false;
value = [];
if isempty(data)
    value = NaN(n,1);
    return
end
if istable(data)
    found = any(strcmp(data.Properties.VariableNames,name));
    if found
        value = data.(name);
    end
elseif isstruct(data)
    found = isfield(data,name);
    if found
        if numel(data) == 1
            value = data.(name);
        else
            value = {data.(name)}.';
        end
    end
end
if ~found
    value = NaN(n,1);
    return
end
actualN = blockHeight(data);
value = fitValue(value,actualN);
value = padRows(value,n);
end

function value = fitValue(value,n)
if iscell(value)
    if isempty(value)
        value = NaN(n,1);
    elseif all(cellfun(@(x)isnumeric(x) && isscalar(x),value))
        value = cell2mat(value(:));
    else
        value = string(value(:));
    end
end
if ischar(value) && (n > 1 || size(value,1) > 1)
    value = string(cellstr(value));
end
if isrow(value) && ~(ischar(value) && n == 1)
    value = value(:);
end
if n == 0
    value = value([]);
elseif numel(value) == 1 && n > 1
    value = repmat(value,n,1);
elseif numel(value) < n
    if isstring(value)
        value(end+1:n,1) = "";
    elseif islogical(value)
        value(end+1:n,1) = false;
    elseif isdatetime(value)
        value(end+1:n,1) = NaT;
    else
        value(end+1:n,1) = NaN;
    end
elseif numel(value) > n
    value = value(1:n);
end
end

function value = padRows(value,n)
if n == 0
    value = value([]);
elseif numel(value) < n
    if isstring(value)
        value(end+1:n,1) = "";
    elseif islogical(value)
        value(end+1:n,1) = false;
    elseif isdatetime(value)
        value(end+1:n,1) = NaT;
    else
        value(end+1:n,1) = NaN;
    end
elseif numel(value) > n
    value = value(1:n);
end
end

function value = firstNumeric(data,n,names,scale)
value = NaN(n,1);
for i = 1:numel(names)
    [candidate,found] = readField(data,names{i},n);
    if found
        value = toNumeric(candidate,n)*scale;
        return
    end
end
end

function value = firstLogical(data,n,names,default)
value = fitLogical(default,n);
for i = 1:numel(names)
    [candidate,found] = readField(data,names{i},n);
    if found
        value = toLogical(candidate,n);
        return
    end
end
end

function value = firstString(data,n,names,default)
value = fitString(default,n);
for i = 1:numel(names)
    [candidate,found] = readField(data,names{i},n);
    if found
        value = string(candidate);
        value = fitString(value,n);
        return
    end
end
end

function value = toNumeric(value,n)
if isnumeric(value) || islogical(value)
    value = double(value);
else
    value = str2double(string(value));
end
value = fitRows(value,n);
end

function value = toLogical(value,n)
if islogical(value)
    value = value;
elseif isnumeric(value)
    value = value ~= 0 & ~isnan(value);
else
    text = lower(strtrim(string(value)));
    value = ismember(text,["true","yes","on","valid","complete","completed"]);
end
value = fitLogical(value,n);
end

function value = fitRows(value,n)
value = value(:);
if n == 0
    value = value([]);
elseif numel(value) == 1 && n > 1
    value = repmat(value,n,1);
elseif numel(value) < n
    value(end+1:n,1) = NaN;
elseif numel(value) > n
    value = value(1:n);
end
end

function value = fitLogical(value,n)
value = logical(value(:));
if n == 0
    value = value([]);
elseif numel(value) == 1 && n > 1
    value = repmat(value,n,1);
elseif numel(value) < n
    value(end+1:n,1) = false;
elseif numel(value) > n
    value = value(1:n);
end
end

function value = fitString(value,n)
value = string(value(:));
if n == 0
    value = value([]);
elseif numel(value) == 1 && n > 1
    value = repmat(value,n,1);
elseif numel(value) < n
    value(end+1:n,1) = "";
elseif numel(value) > n
    value = value(1:n);
end
end

function value = fillMissing(primary,fallback)
primary = primary(:);
fallback = fitRows(fallback,numel(primary));
missing = isnan(primary) & ~isnan(fallback);
primary(missing) = fallback(missing);
value = primary;
end

function value = meanAvailable(varargin)
n = numel(varargin{1});
value = NaN(n,1);
sumValue = zeros(n,1);
count = zeros(n,1);
for i = 1:nargin
    candidate = fitRows(varargin{i},n);
    available = finite(candidate);
    sumValue(available) = sumValue(available) + candidate(available);
    count(available) = count(available) + 1;
end
hasValue = count > 0;
value(hasValue) = sumValue(hasValue)./count(hasValue);
end

function value = finite(x)
value = ~isnan(x) & ~isinf(x);
end

function violation = boundViolationFromMinLoad(minLoad)
minLoad = minLoad(:);
violation = NaN(numel(minLoad),1);
known = finite(minLoad);
violation(known) = max(-minLoad(known),0);
end
