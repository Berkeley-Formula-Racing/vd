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
% The canonical executor owns both fixed and adaptive scheduling. The legacy
% helpers below remain temporarily available for saved-study compatibility, but
% no longer recursively call this entry point.
run = rampSpeed.runCanonicalLongitudinalRamp(car,settings,caseInfo,callbacks);
return
[solverProfile,settings] = rampSpeed.resolveSolverProfileFromSettings(settings);
car = rampSpeed.applySolverProfileToCar(car,solverProfile);
if adaptiveRequested(settings)
    run = runAdaptiveLongitudinalRamp(car,settings,caseInfo,callbacks);
    return
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
    "retrySpeeds_mps",zeros(0,1), ...
    "solver",solverOptions, ...
    "solverProfile",rampSpeed.serializeSolverProfile(solverProfile), ...
    "requestedSpeeds_mps",speeds, ...
    "lateralMetricsApplicable",false);
if solverProfile.approximate
    runMeta.warnings(end+1,1) = ...
        "Approximate aero preview: ride-height aero iteration is disabled.";
end

n = numel(speeds);
rows = repmat(blankRow(NaN,0,"not solved"),n,1);
diagnosticRows = repmat(blankDiagnostic(),n,1);
speedErrors = emptySpeedErrors();
previousState = [];
if isfield(settings,"adaptiveInitialState") && ...
        ~isempty(settings.adaptiveInitialState)
    candidateState = double(settings.adaptiveInitialState(:).');
    if numel(candidateState) ~= 9 || any(~isfinite(candidateState))
        error("rampSpeed:invalidAdaptiveInitialState", ...
            "adaptiveInitialState must contain nine finite state values.");
    end
    previousState = candidateState;
end
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
        initialState = previousState;
        if isempty(initialState)
            [~,longAccel,candidateState,diagnostics] = ...
                max_long_accel(speed,car,[],solverOptions);
        else
            [~,longAccel,candidateState,diagnostics] = ...
                max_long_accel(speed,car,initialState,solverOptions);
        end
        [isFeasible,failureMessage] = solverDiagnosticsFeasible( ...
            diagnostics,solverOptions);
        retriedFromDefault = false;
        tightRetry = false;
        if ~isFeasible && ~isempty(initialState)
            retryDiagnostics = [];
            try
                [~,retryLongAccel,retryState,retryDiagnostics] = ...
                    max_long_accel(speed,car,[],solverOptions);
                [retryFeasible,retryMessage] = solverDiagnosticsFeasible( ...
                    retryDiagnostics,solverOptions);
            catch retryError
                retryFeasible = false;
                retryMessage = "fresh-state retry failed: " + ...
                    string(retryError.message);
            end
            if retryFeasible
                longAccel = retryLongAccel;
                candidateState = retryState;
                diagnostics = retryDiagnostics;
                isFeasible = true;
                failureMessage = "";
                retriedFromDefault = true;
            elseif ~isempty(retryDiagnostics)
                diagnostics = retryDiagnostics;
                failureMessage = retryMessage;
            end
        end
        if ~isFeasible && tightRetryEligible(diagnostics,solverOptions)
            [tightFeasible,tightLongAccel,tightState,tightDiagnostics, ...
                tightMessage] = tightConstraintRetry(speed,car,solverOptions);
            if tightFeasible
                longAccel = tightLongAccel;
                candidateState = tightState;
                diagnostics = tightDiagnostics;
                isFeasible = true;
                failureMessage = "";
                retriedFromDefault = true;
                tightRetry = true;
            elseif ~isempty(tightDiagnostics)
                diagnostics = tightDiagnostics;
                failureMessage = tightMessage;
            end
        end
        if isFeasible
            previousState = candidateState;
            rows(i) = solvedRow(speed,i,longAccel,diagnostics,car, ...
                stateTolerance);
            if tightRetry
                rows(i).accepted_near_feasible = true;
            end
            diagnosticRows(i) = solvedDiagnostic(speed,i,diagnostics);
            completedSpeeds = completedSpeeds + 1;
            if retriedFromDefault
                runMeta.retrySpeeds_mps(end+1,1) = speed;
                if tightRetry
                    warning = "speed " + string(speed) + ...
                        " m/s required a tight-tolerance cold-state solver retry.";
                else
                    warning = "speed " + string(speed) + ...
                        " m/s required a fresh-state solver retry.";
                end
                runMeta.warnings(end+1,1) = warning;
            end
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

normalizationSettings = longitudinalNormalizationSettings(settings,solverOptions);
run = rampSpeed.normalizeRampResult(raw,"longitudinal",normalizationSettings, ...
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

function normalizationSettings = longitudinalNormalizationSettings(settings,solverOptions)
% The shared app residualTolerance is for lateral ramp bisection. Pure
% longitudinal feasibility is gated by the optimizer constraint tolerance.
normalizationSettings = settings;
constraintTolerance = numericSetting(solverOptions, ...
    {"constraintTolerance"},1e-2);
normalizationSettings.residualTolerance = constraintTolerance;
normalizationSettings.ceqTol = constraintTolerance;
normalizationSettings.constraintTolerance = constraintTolerance;
normalizationSettings.inequalityTolerance = constraintTolerance;
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

function [tf,message] = solverDiagnosticsFeasible(diagnostics,solverOptions,requireConverged)
if nargin < 3, requireConverged = true; end
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
if requireConverged
    validExitflag = isnumeric(exitflag) && isfinite(exitflag) && ...
        any(exitflag == [1 2]);
else
    validExitflag = isnumeric(exitflag) && isfinite(exitflag);
end
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

function tf = tightRetryEligible(diagnostics,solverOptions)
% Only spend the expensive cold retry on a numerically near-feasible result.
% A deliberately starved optimizer must remain a failed diagnostic rather
% than being rescued by a fresh retry budget.
tf = false;
if ~isstruct(diagnostics) || ~isscalar(diagnostics)
    return
end
requiredFields = ["exitflag","state", ...
    "max_equality_residual","max_inequality_violation"];
if ~all(isfield(diagnostics,requiredFields))
    return
end
constraintTolerance = numericSetting(solverOptions, ...
    {"constraintTolerance"},1e-2);
if ~isscalar(constraintTolerance) || ~isfinite(constraintTolerance) || ...
        constraintTolerance <= 0
    return
end
state = diagnostics.state;
eqResidual = diagnostics.max_equality_residual;
ineqViolation = diagnostics.max_inequality_violation;
tf = isnumeric(state) && ~isempty(state) && all(isfinite(state(:))) && ...
    isnumeric(diagnostics.exitflag) && isscalar(diagnostics.exitflag) && ...
    isfinite(diagnostics.exitflag) && isnumeric(eqResidual) && ...
    isscalar(eqResidual) && isfinite(eqResidual) && ...
    isnumeric(ineqViolation) && isscalar(ineqViolation) && ...
    isfinite(ineqViolation) && eqResidual <= 2*constraintTolerance && ...
    ineqViolation <= 2*constraintTolerance;
end

function [tf,longAccel,state,diagnostics,message] = ...
        tightConstraintRetry(speed,car,solverOptions)
tf = false;
longAccel = NaN;
state = [];
diagnostics = [];
message = "";
if ~isscalar(speed) || ~isfinite(speed)
    message = "tight-state retry skipped for a nonfinite speed.";
    return
end
constraintTolerance = numericSetting(solverOptions, ...
    {"constraintTolerance"},1e-2);
if ~isscalar(constraintTolerance) || ~isfinite(constraintTolerance) || ...
        constraintTolerance <= 0
    message = "tight-state retry skipped for an invalid constraint tolerance.";
    return
end
options = solverOptions;
options.constraintTolerance = max(1e-8,min(1e-4, ...
    0.01*constraintTolerance));
baseEvaluations = numericSetting(solverOptions, ...
    {"maxFunctionEvaluations"},2000);
options.maxFunctionEvaluations = min(4000,max(baseEvaluations, ...
    2*baseEvaluations));
options.stepTolerance = min(1e-12,numericSetting( ...
    solverOptions,{"stepTolerance"},1e-10));
try
    [~,longAccel,state,diagnostics] = ...
        max_long_accel(speed,car,[],options);
    acceptOptions = solverOptions;
    acceptOptions.constraintTolerance = 0.5*constraintTolerance;
    [tf,message] = solverDiagnosticsFeasible(diagnostics, ...
        acceptOptions,false);
    if ~tf && isempty(message)
        message = "tight-state retry did not meet the acceptance tolerance.";
    end
catch ME
    message = "tight-state retry failed: " + string(ME.message);
end
end
function row = blankRow(speed,speedIndex,reason)
row = struct();
row.speed_mps = speed;
row.speed_index = speedIndex;
row.point_index = 1;
row.valid = false;
row.status = "failed";
row.accepted_near_feasible = false;
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
function tf = adaptiveRequested(settings)
mode = "fixed";
if isfield(settings,"speedGrid") && ~isempty(settings.speedGrid)
    if isstruct(settings.speedGrid) && isfield(settings.speedGrid,"mode")
        mode = lower(strtrim(string(settings.speedGrid.mode)));
    elseif isstring(settings.speedGrid) || ischar(settings.speedGrid)
        mode = lower(strtrim(string(settings.speedGrid)));
    end
elseif isfield(settings,"adaptiveSpeedMode") && ~isempty(settings.adaptiveSpeedMode)
    mode = lower(strtrim(string(settings.adaptiveSpeedMode)));
end
tf = isscalar(mode) && ~ismissing(mode) && ...
    any(mode == ["preview","accurate","highaccuracy"]);
end

function run = runAdaptiveLongitudinalRamp(car,settings,caseInfo,callbacks)
% Execute a deterministic seed grid, then refine only where the response
% changes enough to affect balance or where a solver result needs recovery.
policy = adaptivePolicyFromSettings(settings);
requested = requestedSpeeds(settings);
grid = speedGridSettings(settings);
if isfield(grid,"range_mps") && ~isempty(grid.range_mps)
    requested = double(grid.range_mps(:));
end
if isempty(requested)
    requested = (5:2.5:30).';
end
plan = rampSpeed.planAdaptiveSpeeds(requested,policy);
originalGrid = grid;
fixedSettings = settings;
fixedSettings.speedGrid = struct("mode","fixed");
fixedSettings.speeds = plan.seedSpeeds_mps.';
run = rampSpeed.runLongitudinalRamp(car,fixedSettings,caseInfo,callbacks);
exactRetrySpeeds = zeros(0,1);
[run,exactRetrySpeeds] = retryInvalidAdaptiveSpeeds( ...
    run,car,settings,caseInfo,callbacks,exactRetrySpeeds);
notifyAdaptiveProgress(callbacks,sprintf( ...
    'Adaptive seed scan complete: %d speeds.',height(run.perSpeed)));
if string(run.status) == "cancelled"
    plan.stopReason = "cancelled";
    plan.status = "limited";
else
while plan.pass < policy.maxPasses
    [plan,report] = rampSpeed.refineAdaptiveSpeeds(plan,run);
    if isempty(report.insertedSpeeds_mps)
        if ~isempty(report.invalidSpeeds_mps)
            plan.stopReason = "invalid_speeds";
            plan.status = "limited";
            if ~isempty(plan.refinementHistory)
                plan.refinementHistory(end).stopReason = plan.stopReason;
            end
        end
        notifyAdaptiveProgress(callbacks,sprintf( ...
            'Adaptive refinement stopped after %d speeds (%s).', ...
            height(run.perSpeed),string(plan.stopReason)));
        break
    end
    notifyAdaptiveProgress(callbacks,sprintf( ...
        'Adaptive pass %d: adding %d speeds.', ...
        plan.pass,numel(report.insertedSpeeds_mps)));
    batchSettings = settings;
    batchSettings.speedGrid = struct("mode","fixed");
    batchSettings.speeds = report.insertedSpeeds_mps.';
    seedState = adaptiveSeedState(run,min(report.insertedSpeeds_mps));
    if ~isempty(seedState)
        batchSettings.adaptiveInitialState = seedState;
    end
    batch = rampSpeed.runLongitudinalRamp(car,batchSettings,caseInfo,callbacks);
    run = mergeAdaptiveBatch(run,batch);
    [run,exactRetrySpeeds] = retryInvalidAdaptiveSpeeds( ...
        run,car,settings,caseInfo,callbacks,exactRetrySpeeds);
    notifyAdaptiveProgress(callbacks,sprintf( ...
        'Adaptive pass %d complete: %d speeds.', ...
        plan.pass,height(run.perSpeed)));
    if string(run.status) == "cancelled", break, end
end
end
if plan.pass >= policy.maxPasses && ...
        any(string(plan.stopReason) == ["", "refine"])
    plan.stopReason = "max_passes";
    plan.status = "limited";
    if ~isempty(plan.refinementHistory)
        plan.refinementHistory(end).stopReason = plan.stopReason;
    end
end
plan.retrySpeeds_mps = exactRetrySpeeds;
run.settings = settings;
run.settings.speeds = plan.speeds_mps.';
run.settings.speedGrid = originalGrid;
run.runMeta.requestedSpeeds_mps = plan.requestedSpeeds_mps;
run.runMeta.speedGrid = struct( ...
    "policy",plan.policy, ...
    "requestedSpeeds_mps",plan.requestedSpeeds_mps, ...
    "seedSpeeds_mps",plan.seedSpeeds_mps, ...
    "finalSpeeds_mps",plan.speeds_mps, ...
    "passes",plan.pass, ...
    "provenance",plan.provenance, ...
    "refinementHistory",plan.refinementHistory, ...
    "stopReason",string(plan.stopReason), ...
    "exactRetrySpeeds_mps",exactRetrySpeeds);
run.raw.settings = run.settings;
end

function notifyAdaptiveProgress(callbacks,message)
if ~isstruct(callbacks) || ~isscalar(callbacks) || ...
        ~isfield(callbacks,"onProgress") || isempty(callbacks.onProgress)
    return
end
event = struct("phase","adaptive","speedIndex",NaN,"speed_mps",NaN, ...
    "completedSpeeds",NaN,"requestedSpeeds",NaN,"message",string(message));
try
    callbacks.onProgress(event);
catch
    % Progress is advisory and must not abort an adaptive solve.
end
end

function [run,attempted] = retryInvalidAdaptiveSpeeds( ...
        run,car,settings,caseInfo,callbacks,attempted)
if ~isstruct(run) || ~isscalar(run) || ~isfield(run,"perSpeed") || ...
        ~istable(run.perSpeed) || ~ismember("valid",run.perSpeed.Properties.VariableNames) || ...
        ~ismember("speed_mps",run.perSpeed.Properties.VariableNames) || ...
        string(run.status) == "cancelled"
    return
end
invalid = double(run.perSpeed.speed_mps(~logical(run.perSpeed.valid)));
invalid = sort(invalid(isfinite(invalid)));
retrySpeeds = zeros(0,1);
for value = invalid(:).'
    if ~hasAdaptiveSpeed(attempted,value)
        retrySpeeds(end+1,1) = value; %#ok<AGROW>
    end
end
if isempty(retrySpeeds)
    return
end
notifyAdaptiveProgress(callbacks,sprintf( ...
    'Adaptive exact retry: %d invalid speeds.',numel(retrySpeeds)));
for value = retrySpeeds(:).'
    retrySettings = settings;
    retrySettings.speedGrid = struct("mode","fixed");
    retrySettings.speeds = value;
    retrySettings.adaptiveInitialState = [];
    batch = rampSpeed.runLongitudinalRamp(car,retrySettings,caseInfo,callbacks);
    run = mergeAdaptiveBatch(run,batch);
    attempted(end+1,1) = value; %#ok<AGROW>
    if string(run.status) == "cancelled"
        break
    end
end
notifyAdaptiveProgress(callbacks,sprintf( ...
    'Adaptive exact retry complete: %d speeds.',numel(retrySpeeds)));
end

function tf = hasAdaptiveSpeed(values,speed)
tf = false;
for value = double(values(:).')
    if isfinite(value) && isfinite(speed) && ...
            abs(value-speed) <= 32*eps(max([1,abs(value),abs(speed)]))
        tf = true;
        return
    end
end
end

function state = adaptiveSeedState(run,nextSpeed)
state = [];
if ~isfield(run,"perSpeed") || ~istable(run.perSpeed) || ...
        ~ismember("valid",run.perSpeed.Properties.VariableNames) || ...
        ~ismember("speed_mps",run.perSpeed.Properties.VariableNames)
    return
end
valid = logical(run.perSpeed.valid) & ...
    double(run.perSpeed.speed_mps) < nextSpeed;
if ~any(valid) || ~isfield(run,"raw") || ...
        ~isfield(run.raw,"diagnostics")
    return
end
priorSpeed = max(double(run.perSpeed.speed_mps(valid)));
diagnostics = run.raw.diagnostics;
if isempty(diagnostics) || ~isfield(diagnostics,"speed_mps")
    return
end
index = find(abs([diagnostics.speed_mps]-priorSpeed) <= ...
    32*eps(max(1,abs(priorSpeed))),1);
if ~isempty(index) && isfield(diagnostics,"state") && ...
        numel(diagnostics(index).state) == 9 && ...
        all(isfinite(diagnostics(index).state))
    state = diagnostics(index).state;
end
end

function run = mergeAdaptiveBatch(run,batch)
run.perSpeed = mergeSpeedTables(run.perSpeed,batch.perSpeed);
run.points = mergeSpeedTables(run.points,batch.points);
run.raw.perSpeed = mergeSpeedTables(run.raw.perSpeed,batch.raw.perSpeed);
run.raw.points = mergeSpeedTables(run.raw.points,batch.raw.points);
run.raw.diagnostics = mergeDiagnosticRecords( ...
    run.raw.diagnostics,batch.raw.diagnostics);
batchSpeeds = double(batch.perSpeed.speed_mps(:));
run.raw.speedErrors = mergeSpeedErrors(run.raw.speedErrors, ...
    batch.raw.speedErrors,batchSpeeds);
run.runMeta.warnings = [string(run.runMeta.warnings(:)); ...
    string(batch.runMeta.warnings(:))];
run.runMeta.errors = [string(run.runMeta.errors(:)); ...
    string(batch.runMeta.errors(:))];
run.runMeta.retrySpeeds_mps = unique([run.runMeta.retrySpeeds_mps(:); ...
    batch.runMeta.retrySpeeds_mps(:)]);
run.runMeta.speedErrors = mergeSpeedErrors(run.runMeta.speedErrors, ...
    batch.runMeta.speedErrors,batchSpeeds);
if string(run.status) == "cancelled" || string(batch.status) == "cancelled"
    run.status = "cancelled";
elseif ~any(run.perSpeed.valid)
    run.status = "failed";
elseif any(~run.perSpeed.valid)
    run.status = "partial";
else
    run.status = "completed";
end
run.runMeta.status = run.status;
run.raw.status = run.status;
end

function merged = mergeSpeedTables(primary,added)
if isempty(primary)
    merged = added;
    return
elseif isempty(added)
    merged = primary;
end
merged = sortrows([primary;added],"speed_mps");
speeds = double(merged.speed_mps(:));
keep = true(height(merged),1);
for i = 1:numel(speeds)-1
    if sameAdaptiveSpeed(speeds(i),speeds(i+1))
        keep(i) = false;
    end
end
merged = merged(keep,:);
end

function records = mergeDiagnosticRecords(primary,added)
records = [primary(:);added(:)];
if isempty(records)
    return
end
[~,order] = sort([records.speed_mps]);
records = records(order);
keep = true(numel(records),1);
for i = 1:numel(records)-1
    if sameAdaptiveSpeed(records(i).speed_mps,records(i+1).speed_mps)
        keep(i) = false;
    end
end
records = records(keep);
end

function merged = mergeSpeedErrors(primary,added,replacedSpeeds)
primary = primary(:);
added = added(:);
keep = true(numel(primary),1);
for i = 1:numel(primary)
    if isfield(primary,"speed_mps")
        for speed = replacedSpeeds(:).'
            if sameAdaptiveSpeed(primary(i).speed_mps,speed)
                keep(i) = false;
                break
            end
        end
    end
end
merged = [primary(keep);added];
end

function tf = sameAdaptiveSpeed(a,b)
tf = isfinite(a) && isfinite(b) && ...
    abs(double(a)-double(b)) <= 32*eps(max([1,abs(double(a)),abs(double(b))]));
end

function policy = adaptivePolicyFromSettings(settings)
grid = speedGridSettings(settings);
mode = "accurate";
if isfield(grid,"mode") && ~isempty(grid.mode)
    mode = string(grid.mode);
end
overrides = struct();
if isfield(grid,"policy") && isstruct(grid.policy)
    overrides = grid.policy;
end
allowed = ["baseSpacing_mps","minRefinementSpacing_mps", ...
    "balanceTolerance_fraction","relativeTolerance", ...
    "forceAbsoluteTolerance_N","accelerationAbsoluteTolerance_mps2", ...
    "residualRefinementFraction","maxPasses","maxPoints","stateFields"];
for name = allowed
    fieldName = char(name);
    if isfield(grid,fieldName)
        overrides.(fieldName) = grid.(fieldName);
    end
end
policy = rampSpeed.adaptiveSpeedPolicy(mode,overrides);
end

function grid = speedGridSettings(settings)
grid = struct("mode","fixed");
if isfield(settings,"speedGrid") && ~isempty(settings.speedGrid)
    if isstruct(settings.speedGrid)
        grid = settings.speedGrid;
    else
        grid.mode = settings.speedGrid;
    end
elseif isfield(settings,"adaptiveSpeedMode") && ~isempty(settings.adaptiveSpeedMode)
    grid.mode = settings.adaptiveSpeedMode;
end
if ~isstruct(grid) || ~isscalar(grid)
    error("rampSpeed:invalidAdaptivePolicy", ...
        "settings.speedGrid must be a scalar struct or mode name.");
end
end
