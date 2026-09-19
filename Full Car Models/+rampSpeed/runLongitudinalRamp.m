function run = runLongitudinalRamp(car,settings,caseInfo,callbacks)
%RUNLONGITUDINALRAMP Run pure-longitudinal acceleration at requested speeds.

if nargin < 2 || isempty(settings)
    settings = struct();
end
if nargin < 3 || isempty(caseInfo)
    caseInfo = struct();
end
if nargin < 4 || isempty(callbacks)
    callbacks = struct();
end
if ~isstruct(settings) || ~isscalar(settings)
    error("rampSpeed:invalidSettings", ...
        "settings must be a scalar struct.");
end
if ~isstruct(caseInfo) || ~isscalar(caseInfo)
    error("rampSpeed:invalidCaseInfo", ...
        "caseInfo must be a scalar struct.");
end
if ~isstruct(callbacks) || ~isscalar(callbacks)
    error("rampSpeed:invalidCallbacks", ...
        "callbacks must be a scalar struct.");
end

speeds = requestedSpeeds(settings);
if isempty(speeds)
    speeds = (5:2.5:30).';
end
settings.speeds = speeds.';
solverOptions = solverOptionsFromSettings(settings);
stateTolerance = numericSetting(settings, ...
    {"zeroStateTolerance","stateTolerance"},1e-12);

progressFcn = [];
cancelFcn = [];
if isfield(callbacks,"onProgress") && ~isempty(callbacks.onProgress)
    progressFcn = callbacks.onProgress;
end
if isfield(callbacks,"isCancelled") && ~isempty(callbacks.isCancelled)
    cancelFcn = callbacks.isCancelled;
end

started = datetime('now');
runMeta = struct("source","rampSpeed.runLongitudinalRamp", ...
    "started",started,"completed",datetime.empty, ...
    "warnings",strings(0,1),"errors",strings(0,1), ...
    "solver",solverOptions, ...
    "requestedSpeeds_mps",speeds, ...
    "lateralMetricsApplicable",false);

n = numel(speeds);
rows = repmat(blankRow(NaN,0,"not solved"),n,1);
diagnosticRows = repmat(blankDiagnostic(),n,1);
speedErrors = emptySpeedErrors();
previousState = [];
cancelled = false;
completedSpeeds = 0;

for i = 1:n
    speed = speeds(i);
    if cancellationRequested(cancelFcn)
        cancelled = true;
        [rows,diagnosticRows] = fillCancelledRows(rows,diagnosticRows, ...
            speeds,i);
        break
    end

    if ~isempty(progressFcn)
        progressFcn(struct("phase","speed","speedIndex",i, ...
            "speed_mps",speed,"completedSpeeds",completedSpeeds, ...
            "requestedSpeeds",n));
    end

    try
        if isempty(previousState)
            [~,longAccel,candidateState,diagnostics] = ...
                max_long_accel(speed,car,[],solverOptions);
        else
            [~,longAccel,candidateState,diagnostics] = ...
                max_long_accel(speed,car,previousState,solverOptions);
        end
        [isFeasible,failureMessage] = solverDiagnosticsFeasible( ...
            diagnostics,solverOptions);
        if isFeasible
            previousState = candidateState;
            rows(i) = solvedRow(speed,i,longAccel,diagnostics,car, ...
                stateTolerance);
            diagnosticRows(i) = solvedDiagnostic(speed,i,diagnostics);
            completedSpeeds = completedSpeeds + 1;
            if ~isempty(progressFcn)
                progressFcn(struct("phase","speed","speedIndex",i, ...
                    "speed_mps",speed, ...
                    "completedSpeeds",completedSpeeds, ...
                    "requestedSpeeds",n));
            end
        else
            previousState = [];
            failureIdentifier = "rampSpeed:nonconverged";
            rows(i) = blankRow(speed,i, ...
                failureIdentifier + ": " + failureMessage);
            speedErrors(end+1) = makeDiagnosticSpeedError(i,speed, ...
                failureIdentifier,failureMessage); %#ok<AGROW>
            diagnosticRows(i) = failedDiagnostic(speed,i,diagnostics, ...
                failureIdentifier,failureMessage,[]);
            runMeta.errors(end+1,1) = failureMessage;
        end
    catch ME
        rows(i) = blankRow(speed,i, ...
            "speed solve failed: " + string(ME.message));
        speedErrors(end+1) = makeSpeedError(i,speed,ME); %#ok<AGROW>
        diagnosticRows(i) = failedDiagnostic(speed,i,[], ...
            string(ME.identifier),string(ME.message),ME.stack);
        runMeta.errors(end+1,1) = string(ME.message);
        previousState = [];
    end

    if ~cancelled && cancellationRequested(cancelFcn)
        cancelled = true;
        [rows,diagnosticRows] = fillCancelledRows(rows,diagnosticRows, ...
            speeds,i+1);
        break
    end
end

validRows = [rows.valid].';
if cancelled
    runStatus = "cancelled";
elseif isempty(validRows)
    runStatus = "completed";
elseif ~any(validRows)
    runStatus = "failed";
elseif any(~validRows)
    runStatus = "partial";
else
    runStatus = "completed";
end
pointRows = rows;
runMeta.status = runStatus;
runMeta.completed = datetime('now');
runMeta.speedErrors = speedErrors;
if ~isempty(speedErrors)
    for i = 1:numel(speedErrors)
        runMeta.warnings(end+1,1) = sprintf( ...
            "speed %g m/s failed: %s",speedErrors(i).speed_mps, ...
            speedErrors(i).message);
    end
end

raw = struct();
raw.perSpeed = struct2table(rows);
raw.points = struct2table(pointRows);
raw.settings = settings;
raw.status = runStatus;
raw.diagnostics = diagnosticRows;
raw.speedErrors = speedErrors;

run = rampSpeed.normalizeRampResult(raw,"longitudinal",settings, ...
    caseInfo,runMeta);
run.settings = settings;
run.settings.speeds = speeds.';
run.runMeta.requestedSpeeds_mps = speeds;
run.runMeta.lateralMetricsApplicable = false;
run.runMeta.speedErrors = speedErrors;
invalid = ~([rows.valid].');
if any(invalid)
    reasons = strings(n,1);
    for i = 1:n
        reasons(i) = string(rows(i).reason);
    end
    run.perSpeed.speed_mps = speeds;
    run.perSpeed.valid(invalid) = false;
    run.perSpeed.status(invalid) = "invalid";
    run.perSpeed.reason(invalid) = reasons(invalid);
    run.points.speed_mps = speeds;
    run.points.valid(invalid) = false;
    run.points.status(invalid) = "invalid";
end
run.status = runStatus;
end

function speeds = requestedSpeeds(settings)
if isfield(settings,"speeds") && ~isempty(settings.speeds)
    speeds = double(settings.speeds(:));
else
    speeds = zeros(0,1);
end
end

function options = solverOptionsFromSettings(settings)
options = struct();
if isfield(settings,"solverOptions") && ~isempty(settings.solverOptions)
    if ~isstruct(settings.solverOptions) || ~isscalar(settings.solverOptions)
        error("rampSpeed:invalidSolverOptions", ...
            "settings.solverOptions must be a scalar struct.");
    end
    options = settings.solverOptions;
end
names = ["maxFunctionEvaluations","constraintTolerance", ...
    "stepTolerance","display"];
for i = 1:numel(names)
    name = char(names(i));
    if isfield(settings,name) && ~isfield(options,name)
        options.(name) = settings.(name);
    end
end
end

function value = numericSetting(settings,names,default)
value = default;
for i = 1:numel(names)
    name = char(names{i});
    if isfield(settings,name) && ~isempty(settings.(name))
        value = double(settings.(name));
        value = value(1);
        return
    end
end
end

function tf = cancellationRequested(cancelFcn)
tf = false;
if ~isempty(cancelFcn)
    value = cancelFcn();
    if ~isscalar(value)
        error("rampSpeed:invalidCancellation", ...
            "callbacks.isCancelled must return a scalar logical.");
    end
    tf = logical(value);
end
end

function [rows,diagnosticRows] = fillCancelledRows(rows,diagnosticRows, ...
    speeds,startIndex)
for j = startIndex:numel(speeds)
    rows(j) = blankRow(speeds(j),j, ...
        "cancelled before this speed was solved");
    diagnosticRows(j) = failedDiagnostic(speeds(j),j,[], ...
        "rampSpeed:cancelled", ...
        "cancelled before this speed was solved",[]);
end
end

function [tf,message] = solverDiagnosticsFeasible(diagnostics,solverOptions)
exitflag = NaN;
maxInequalityViolation = NaN;
maxEqualityResidual = NaN;
state = [];
if isstruct(diagnostics) && isscalar(diagnostics)
    if isfield(diagnostics,"exitflag")
        exitflag = diagnostics.exitflag;
    end
    if isfield(diagnostics,"max_inequality_violation")
        maxInequalityViolation = diagnostics.max_inequality_violation;
    end
    if isfield(diagnostics,"max_equality_residual")
        maxEqualityResidual = diagnostics.max_equality_residual;
    end
    if isfield(diagnostics,"state")
        state = diagnostics.state;
    end
end
if ~isscalar(exitflag)
    exitflag = NaN;
end
if ~isscalar(maxInequalityViolation)
    maxInequalityViolation = NaN;
end
if ~isscalar(maxEqualityResidual)
    maxEqualityResidual = NaN;
end
constraintTolerance = numericSetting(solverOptions, ...
    {"constraintTolerance"},1e-2);
validExitflag = isnumeric(exitflag) && isfinite(exitflag) && ...
    any(exitflag == [1 2]);
finiteState = isnumeric(state) && ~isempty(state) && ...
    all(isfinite(state(:)));
finiteResiduals = isnumeric(maxInequalityViolation) && ...
    isfinite(maxInequalityViolation) && ...
    isnumeric(maxEqualityResidual) && isfinite(maxEqualityResidual);
validTolerance = isnumeric(constraintTolerance) && ...
    isscalar(constraintTolerance) && isfinite(constraintTolerance) && ...
    constraintTolerance >= 0;
tf = validExitflag && finiteState && finiteResiduals && ...
    validTolerance && ...
    maxInequalityViolation <= constraintTolerance && ...
    maxEqualityResidual <= constraintTolerance;
if tf
    message = "";
else
    message = sprintf(['nonconverged solver result: exitflag=%g, ' ...
        'max equality residual=%g, max inequality violation=%g, ' ...
        'tolerance=%g'],exitflag,maxEqualityResidual, ...
        maxInequalityViolation,constraintTolerance);
end
end

function row = blankRow(speed,speedIndex,reason)
row = struct();
row.speed_mps = speed;
row.speed_index = speedIndex;
row.point_index = 1;
row.valid = false;
row.status = "failed";
row.reason = string(reason);
row.exitflag = NaN;
row.max_constraint_residual = NaN;
row.max_equality_residual = NaN;
row.max_inequality_violation = NaN;
row.aLat_mps2 = NaN;
row.aLong_mps2 = NaN;
row.steer_rad = NaN;
row.lat_velocity_mps = NaN;
row.yaw_rate_rps = NaN;
row.engine_rpm = NaN;
row.current_gear = NaN;
row.throttle = NaN;
row.downforce_N = NaN;
row.drag_N = NaN;
row.ClA_m2 = NaN;
row.CdA_m2 = NaN;
row.LoD = NaN;
row.aero_balance_front = NaN;
row.aero_outside_map = false;
row.aero_residual_m = NaN;
row.Fz_front_axle_N = NaN;
row.Fz_rear_axle_N = NaN;
row.min_Fz_N = NaN;
row.wheel_lift = false;
row.LLTD = NaN;
row.LLT_front_N = NaN;
row.LLT_rear_N = NaN;
row.long_load_transfer_N = NaN;
row.aLong_max_mps2 = NaN;
row.aLat_achieved_mps2 = NaN;
row.aLat_force_residual_mps2 = NaN;
row.pure_ay0 = false;
row.steer_zero = false;
row.lat_velocity_zero = false;
row.yaw_rate_zero = false;
row.rear_slip_ratio = NaN;
row.throttle_upper_active = false;
row.rear_slip_upper_active = false;
row.power_limited = false;
row.traction_limited = false;
row.lateral_metrics_applicable = false;
row.aLat_free_mps2 = NaN;
row.aLat_sustainable_mps2 = NaN;
row.ramp_complete_fraction = NaN;
row.truncated = false;
row.K_linear_rad_per_mps2 = NaN;
row.K_linear_r2 = NaN;
row.K_at_limit_rad_per_mps2 = NaN;
row.cuo_steer_linear_rad = NaN;
row.cuo_steer_limit_rad = NaN;
row.mechanical_balance_front = NaN;
row.grip_balance_mid = NaN;
row.grip_balance_limit = NaN;
row.alpha_balance_mid_rad = NaN;
row.alpha_balance_limit_rad = NaN;
row.LLT_norm_balance_mid = NaN;
row.LLT_norm_balance_limit = NaN;
row.front_Fz_fraction_mid = NaN;
row.front_Fz_fraction_limit = NaN;
row.front_downforce_N = NaN;
row.rear_downforce_N = NaN;
row.front_ride_height_m = NaN;
row.rear_ride_height_m = NaN;
row.pitch_rad = NaN;
row.front_shock_travel_m = NaN;
row.rear_shock_travel_m = NaN;
row.front_camber_rad = NaN;
row.rear_camber_rad = NaN;
row.min_Fz_limit_N = NaN;
row.n_exitflag1 = NaN;
row.n_exitflag2 = NaN;

corners = ["FL","FR","RL","RR"];
for i = 1:numel(corners)
    corner = char(corners(i));
    row.(['Fz_' corner '_N']) = NaN;
    row.(['Fx_' corner '_N']) = NaN;
    row.(['Fy_' corner '_N']) = NaN;
    row.(['alpha_' corner '_rad']) = NaN;
    row.(['gamma_' corner '_rad']) = NaN;
    row.(['kappa_' corner]) = NaN;
    row.(['T_' corner '_Nm']) = NaN;
    row.(['omega_' corner '_rps']) = NaN;
end
end

function row = solvedRow(speed,speedIndex,longAccel,diagnostics,car, ...
    stateTolerance)
row = blankRow(speed,speedIndex,"");
x = diagnostics.state;
m = car.metrics(x);
row.valid = true;
row.status = "complete";
row.reason = "";
row.exitflag = diagnostics.exitflag;
row.max_equality_residual = diagnostics.max_equality_residual;
row.max_inequality_violation = diagnostics.max_inequality_violation;
row.max_constraint_residual = max([ ...
    diagnostics.max_equality_residual, ...
    diagnostics.max_inequality_violation]);
row.aLat_mps2 = m.gLat*car.g;
row.aLong_mps2 = longAccel;
row.steer_rad = x(1)*pi/180;
row.lat_velocity_mps = m.lat_vel;
row.yaw_rate_rps = m.yaw_rate;
row.engine_rpm = m.engine_rpm;
row.current_gear = m.current_gear;
row.throttle = m.throttle;
row.downforce_N = m.downforce;
row.drag_N = m.drag;
row.ClA_m2 = m.ClA;
row.CdA_m2 = m.CdA;
row.LoD = m.LoD;
row.aero_balance_front = m.CoP;
row.aero_outside_map = logical(m.aero_outside_map);
row.aero_residual_m = m.aero_residual_in*0.0254;
row.Fz_front_axle_N = m.Fz_front_axle;
row.Fz_rear_axle_N = m.Fz_rear_axle;
row.min_Fz_N = m.min_Fz;
row.wheel_lift = m.min_Fz <= 0;
row.LLTD = m.LLTD;
row.LLT_front_N = m.LLT_front;
row.LLT_rear_N = m.LLT_rear;
row.long_load_transfer_N = m.long_load_transfer;
row.aLong_max_mps2 = longAccel;
row.aLat_achieved_mps2 = m.gLat*car.g;
row.aLat_force_residual_mps2 = m.lat_accel_residual;
row.pure_ay0 = diagnostics.pure_ay0;
row.steer_zero = abs(row.steer_rad) <= stateTolerance;
row.lat_velocity_zero = abs(row.lat_velocity_mps) <= stateTolerance;
row.yaw_rate_zero = abs(row.yaw_rate_rps) <= stateTolerance;
row.rear_slip_ratio = meanFinite([m.kappa_3,m.kappa_4]);
row.throttle_upper_active = row.throttle >= 1-1e-9;
row.rear_slip_upper_active = row.rear_slip_ratio >= 0.2-1e-9;
row.traction_limited = row.rear_slip_upper_active;
row.power_limited = row.throttle_upper_active && ~row.traction_limited;
row.front_downforce_N = m.aero_downforce_front_N;
row.rear_downforce_N = m.aero_downforce_rear_N;
row.front_ride_height_m = m.front_ride_height_in*0.0254;
row.rear_ride_height_m = m.rear_ride_height_in*0.0254;
row.pitch_rad = atan2((m.front_ride_height_in-m.rear_ride_height_in) ...
    *0.0254,car.W_b);
row.front_shock_travel_m = shockTravelM(car.rideHeightAero, ...
    m.front_ride_height_in,"front");
row.rear_shock_travel_m = shockTravelM(car.rideHeightAero, ...
    m.rear_ride_height_in,"rear");
row.front_camber_rad = meanFinite([m.gamma_1,m.gamma_2])*pi/180;
row.rear_camber_rad = meanFinite([m.gamma_3,m.gamma_4])*pi/180;

corners = ["FL","FR","RL","RR"];
for i = 1:4
    corner = char(corners(i));
    index = num2str(i);
    row.(['Fz_' corner '_N']) = m.(['Fz_' index]);
    row.(['Fx_' corner '_N']) = m.(['Fx_' index]);
    row.(['Fy_' corner '_N']) = m.(['Fy_' index]);
    row.(['alpha_' corner '_rad']) = m.(['alpha_' index])*pi/180;
    row.(['gamma_' corner '_rad']) = m.(['gamma_' index])*pi/180;
    row.(['kappa_' corner]) = m.(['kappa_' index]);
    row.(['T_' corner '_Nm']) = m.(['T_' index]);
    row.(['omega_' corner '_rps']) = m.(['omega_' index]);
end
end

function value = shockTravelM(config,rideHeightIn,axle)
staticName = ['static_' char(axle) '_ride_height_in'];
motionName = ['motion_ratio_' char(axle)];
if ~isstruct(config) || ~isfield(config,staticName) || ...
        ~isfield(config,motionName)
    value = NaN;
    return
end
value = (config.(staticName)-rideHeightIn)*config.(motionName)*0.0254;
end

function value = meanFinite(values)
values = values(isfinite(values));
if isempty(values)
    value = NaN;
else
    value = mean(values);
end
end

function record = blankDiagnostic()
record = struct("speed_index",NaN,"speed_mps",NaN,"success",false, ...
    "state",[],"exitflag",NaN,"c",[],"ceq",[], ...
    "max_inequality_violation",NaN,"max_equality_residual",NaN, ...
    "metrics",struct(),"error_identifier","","error_message","", ...
    "error_stack",[]);
end

function record = solvedDiagnostic(speed,speedIndex,diagnostics)
record = blankDiagnostic();
record.speed_index = speedIndex;
record.speed_mps = speed;
record.success = true;
record.state = diagnostics.state;
record.exitflag = diagnostics.exitflag;
record.c = diagnostics.c;
record.ceq = diagnostics.ceq;
record.max_inequality_violation = diagnostics.max_inequality_violation;
record.max_equality_residual = diagnostics.max_equality_residual;
record.metrics = diagnostics.metrics;
end

function record = failedDiagnostic(speed,speedIndex,diagnostics, ...
    identifier,message,stack)
record = blankDiagnostic();
record.speed_index = speedIndex;
record.speed_mps = speed;
if ~isempty(diagnostics)
    record.state = diagnostics.state;
    record.exitflag = diagnostics.exitflag;
    record.c = diagnostics.c;
    record.ceq = diagnostics.ceq;
    record.max_inequality_violation = ...
        diagnostics.max_inequality_violation;
    record.max_equality_residual = diagnostics.max_equality_residual;
    record.metrics = diagnostics.metrics;
end
record.error_identifier = string(identifier);
record.error_message = string(message);
record.error_stack = stack;
end

function errors = emptySpeedErrors()
errors = struct('speed_index',{},'speed_mps',{},'identifier',{}, ...
    'message',{},'stack',{});
end

function entry = makeSpeedError(speedIndex,speed,ME)
entry = struct('speed_index',speedIndex,'speed_mps',speed, ...
    'identifier',string(ME.identifier),'message',string(ME.message), ...
    'stack',ME.stack);
end

function entry = makeDiagnosticSpeedError(speedIndex,speed,identifier,message)
entry = struct('speed_index',speedIndex,'speed_mps',speed, ...
    'identifier',string(identifier),'message',string(message), ...
    'stack',[]);
end
