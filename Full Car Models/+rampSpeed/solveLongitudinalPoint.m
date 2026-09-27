function result = solveLongitudinalPoint(car,speed_mps,seed,profile,control)
%SOLVELONGITUDINALPOINT Solve one pure-longitudinal ramp-speed point.
%
% The independent variables are throttle and one symmetric rear slip ratio.
% Each gearbox branch is solved explicitly so a gear transition cannot be
% hidden inside the optimizer's discontinuous automatic gear selector.

if nargin < 3, seed = []; end
if nargin < 4 || isempty(profile), profile = "accurate"; end
if nargin < 5 || isempty(control), control = struct(); end
validateattributes(speed_mps,{'numeric'},{'real','finite','scalar','positive'}, ...
    mfilename,'speed_mps');
if ~isobject(car) || ~isscalar(car) || ~isprop(car,'powertrain')
    error('rampSpeed:invalidCar','A scalar solver-ready Car is required.');
end
if ~isstruct(control) || ~isscalar(control)
    error('rampSpeed:invalidControl','control must be a scalar struct.');
end

profile = rampSpeed.resolveSolverProfile(profile);
if isfield(profile,'powertrainModel') && ...
        string(profile.powertrainModel) == "continuousEnvelope"
    result = rampSpeed.solveContinuousEnvelopePoint(car,speed_mps,seed, ...
        profile,control);
    return
end
task = struct('speedIndex',1,'speed_mps',double(speed_mps), ...
    'origin',"requested",'passIndex',1);
emptyMetrics = struct();

if shouldCancel(control)
    result = rampSpeed.makeSpeedResult(task,"cancelled",emptyMetrics, ...
        struct('gearAttempts',struct([]),'evaluationCount',0, ...
        'reducedState',[NaN NaN], ...
        'reason',"cancelled before first evaluation"));
    result.state = NaN(1,9);
    return
end

candidates = rampSpeed.candidateGears(car,speed_mps);
attempts = repmat(attemptTemplate(),numel(candidates),1);
best = struct('found',false,'state',NaN(1,9),'data',struct(), ...
    'gear',NaN,'status',"infeasible",'maxEqualityResidual',Inf, ...
    'maxInequalityViolation',Inf);
evaluationCount = 0;
anySolverFailure = false;
anyEvaluated = false;

for i = 1:numel(candidates)
    candidate = candidates(i);
    attempts(i).gear = candidate.gear;
    attempts(i).predictedEngineRpm = candidate.predictedEngineRpm;
    attempts(i).withinRedline = candidate.withinRedline;
    attempts(i).reason = string(candidate.reason);

    if shouldCancel(control)
        result = cancelledResult(task,attempts,evaluationCount, ...
            "cancelled during gear evaluation",best.state);
        return
    end
    if ~candidate.withinRedline
        attempts(i).status = "infeasible";
        continue
    end

    anyEvaluated = true;
    try
        [state,data,exitflag,exitMessage,solverOutput] = solveGear( ...
            candidate.gear,seed);
        attempts(i).exitflag = exitflag;
        attempts(i).exitMessage = string(exitMessage);
        attempts(i).solverOutput = solverOutput;
        attempts(i).engineRpm = data.engineRpm;
        attempts(i).maxInequalityViolation = data.maxInequalityViolation;
        attempts(i).maxEqualityResidual = data.maxEqualityResidual;
        attempts(i).long_accel_mps2 = data.long_accel_mps2;
        attempts(i).state = state;
        attempts(i).feasible = isAcceptable(exitflag,data, ...
            constraintTolerance(profile));
        if attempts(i).feasible
            if data.maxEqualityResidual <= constraintTolerance(profile) && ...
                    data.maxInequalityViolation <= constraintTolerance(profile)
                attempts(i).status = "converged";
            else
                attempts(i).status = "near_feasible";
            end
            attempts(i).reason = "accepted";
            if ~best.found || data.long_accel_mps2 > best.data.long_accel_mps2
                best.found = true;
                best.state = state;
                best.data = data;
                best.gear = candidate.gear;
                best.status = attempts(i).status;
                best.maxEqualityResidual = data.maxEqualityResidual;
                best.maxInequalityViolation = data.maxInequalityViolation;
            end
        elseif exitflag <= 0
            anySolverFailure = true;
            attempts(i).status = "solver_failed";
            attempts(i).reason = "optimizer did not converge to a feasible point";
        else
            attempts(i).status = "infeasible";
            attempts(i).reason = "physical or equilibrium residual exceeded tolerance";
        end
    catch ME
        if strcmp(ME.identifier,'rampSpeed:cancelled')
            result = cancelledResult(task,attempts,evaluationCount, ...
                "cancelled during gear evaluation",best.state);
            return
        end
        anySolverFailure = true;
        attempts(i).status = "solver_failed";
        attempts(i).reason = string(ME.message);
        attempts(i).exitMessage = string(ME.message);
    end
end

diagnostics = struct();
diagnostics.gearAttempts = attempts;
diagnostics.evaluationCount = evaluationCount;
diagnostics.solverProfile = profile;
diagnostics.referenceFallback = referenceFallback(car,speed_mps,seed, ...
    profile,anyEvaluated,best.found,control);

if best.found
    state = best.state;
    metrics = car.metrics(state,struct('gearOverride',best.gear));
    metrics.aLong_mps2 = best.data.long_accel_mps2;
    metrics.aLong_max_mps2 = best.data.long_accel_mps2;
    metrics.aLat_mps2 = speed_mps*state(5);
    metrics.longitudinal_acceleration_mps2 = best.data.long_accel_mps2;
    metrics.lat_accel_residual = best.data.lat_accel;
    metrics.yaw_accel_residual = best.data.yaw_accel;
    metrics.max_equality_residual = best.data.maxEqualityResidual;
    metrics.max_inequality_violation = best.data.maxInequalityViolation;
    metrics.current_gear = best.gear;
    metrics.rear_slip = state(8);
    diagnostics.reducedState = [state(2),state(8)];
    diagnostics.selectedGear = best.gear;
    diagnostics.selectedEngineRpm = best.data.engineRpm;
    diagnostics.maxEqualityResidual = best.data.maxEqualityResidual;
    diagnostics.maxInequalityViolation = best.data.maxInequalityViolation;
    status = best.status;
else
    state = best.state;
    metrics = emptyMetrics;
    diagnostics.reducedState = [NaN NaN];
    diagnostics.selectedGear = NaN;
    diagnostics.selectedEngineRpm = NaN;
    diagnostics.maxEqualityResidual = NaN;
    diagnostics.maxInequalityViolation = NaN;
    if anySolverFailure
        status = "solver_failed";
    else
        status = "infeasible";
    end
end

result = rampSpeed.makeSpeedResult(task,status,metrics,diagnostics,attempts);
result.state = state;

    function [state,data,exitflag,exitMessage,solverOutput] = solveGear(gear,seedState)
        bounds = [0 1; 0 0.2];
        x0 = initialReducedState(seedState,bounds);
        lastX = NaN(1,2);
        lastData = struct();
        rideHeightContext = [];

        function dataOut = evaluateReduced(z)
            z = double(z(:).');
            if isequal(z,lastX)
                dataOut = lastData;
                return
            end
            if shouldCancel(control)
                error('rampSpeed:cancelled', ...
                    'cancelled during longitudinal point evaluation');
            end
            evaluationCount = evaluationCount + 1;
            stateLocal = rampSpeed.makePureLongState(car,speed_mps,z(1),z(2));
            options = struct('gearOverride',gear);
            [engineRpm,beta,latAccel,longAccel,yawAccel,wheelAccel, ...
                    omega,currentGear,Fzvirtual,Fz,alpha,T,Fy,gamma,Fx, ...
                    ss,rideHeightContext] = car.equations(stateLocal, ...
                    rideHeightContext,options);
            inequality = [engineRpm-car.powertrain.redline,abs(beta)-20, ...
                -Fzvirtual(:).'];
            equality = [latAccel,yawAccel,wheelAccel(:).'];
            if any(~isfinite([engineRpm,beta,longAccel, ...
                    Fzvirtual(:).',equality]))
                maxIneq = Inf;
                maxEq = Inf;
            else
                maxIneq = max([inequality,0]);
                maxEq = max(abs(equality));
            end
            dataOut = struct('state',stateLocal,'engineRpm',engineRpm, ...
                'beta',beta,'lat_accel',latAccel, ...
                'long_accel_mps2',longAccel,'yaw_accel',yawAccel, ...
                'wheel_accel',wheelAccel,'omega',omega, ...
                'currentGear',currentGear,'Fzvirtual',Fzvirtual,'Fz',Fz, ...
                'alpha',alpha,'T',T,'Fy',Fy,'gamma',gamma,'Fx',Fx, ...
                'ssInfo',ss,'inequality',inequality,'equality',equality, ...
                'maxInequalityViolation',maxIneq, ...
                'maxEqualityResidual',maxEq);
            lastX = z;
            lastData = dataOut;
        end

        objective = @objectiveReduced;
        nonlinear = @constraintsReduced;
        options = fminconOptions(profile);
        [z,~,exitflag,solverOutput] = fmincon(objective,x0,[],[],[],[], ...
            bounds(:,1),bounds(:,2),nonlinear,options);
        data = evaluateReduced(z);
        state = data.state;
        exitMessage = solverExitMessage(exitflag,solverOutput);

        function value = objectiveReduced(zLocal)
            dataLocal = evaluateReduced(zLocal);
            if isfinite(dataLocal.long_accel_mps2)
                value = -dataLocal.long_accel_mps2;
            else
                value = realmax('double')/100;
            end
        end

        function [c,ceq] = constraintsReduced(zLocal)
            dataLocal = evaluateReduced(zLocal);
            c = dataLocal.inequality(:);
            % Both rear wheels are constrained to the same slip ratio and
            % have the same pure-longitudinal state, so their equilibrium
            % residuals are duplicates. One scalar equality avoids a
            % singular KKT system while the full pair remains diagnostic.
            ceq = mean(dataLocal.wheel_accel(3:4));
            if any(~isfinite(c)) || any(~isfinite(ceq))
                c = ones(numel(c),1)*realmax('double')/100;
                ceq = ones(numel(ceq),1)*realmax('double')/100;
            end
        end
    end
end

function attempt = attemptTemplate()
attempt = struct('gear',NaN,'predictedEngineRpm',NaN, ...
    'withinRedline',false,'exitflag',NaN,'exitMessage',"", ...
    'solverOutput',struct(),'engineRpm',NaN, ...
    'maxInequalityViolation',NaN,'maxEqualityResidual',NaN, ...
    'long_accel_mps2',NaN,'state',NaN(1,9),'feasible',false, ...
    'status',"not_attempted",'reason',"");
end

function x0 = initialReducedState(seed,bounds)
x0 = [1 0.03];
if isnumeric(seed) && numel(seed) == 9 && all(isfinite(seed(:)))
    seedState = double(seed(:).');
    x0 = [seedState(2),seedState(8)];
elseif isstruct(seed) && isfield(seed,'state') && ...
        isnumeric(seed.state) && numel(seed.state) == 9 && ...
        all(isfinite(seed.state(:)))
    seedState = double(seed.state(:).');
    x0 = [seedState(2),seedState(8)];
end
x0 = min(max(x0,bounds(:,1).'),bounds(:,2).');
end

function options = fminconOptions(profile)
solver = profile.solverOptions;
options = optimoptions('fmincon', ...
    'Algorithm','interior-point', ...
    'MaxFunctionEvaluations',double(solver.maxFunctionEvaluations), ...
    'ConstraintTolerance',double(solver.constraintTolerance), ...
    'StepTolerance',double(solver.stepTolerance), ...
    'Display',char(solver.display));
end

function message = solverExitMessage(exitflag,solverOutput)
message = "";
if isstruct(solverOutput) && isfield(solverOutput,'message')
    message = string(solverOutput.message);
end
if strlength(message) == 0
    message = "exitflag " + string(exitflag);
end
end

function tf = isAcceptable(exitflag,data,tol)
tf = ismember(exitflag,[1 2]) && ...
    isfinite(data.long_accel_mps2) && ...
    data.maxEqualityResidual <= 5*tol && ...
    data.maxInequalityViolation <= 5*tol;
end

function tol = constraintTolerance(profile)
tol = double(profile.solverOptions.constraintTolerance);
end

function value = shouldCancel(control)
value = false;
if isfield(control,'shouldCancel') && ~isempty(control.shouldCancel)
    callback = control.shouldCancel;
    if isa(callback,'function_handle')
        value = callback();
    else
        value = callback;
    end
    if ~isscalar(value)
        error('rampSpeed:invalidCancellation', ...
            'control.shouldCancel must return a scalar logical.');
    end
    value = logical(value);
end
end

function result = cancelledResult(task,attempts,evaluationCount,reason,state)
diagnostics = struct('gearAttempts',attempts,'evaluationCount',evaluationCount, ...
    'reducedState',[NaN NaN],'reason',string(reason));
result = rampSpeed.makeSpeedResult(task,"cancelled",struct(),diagnostics,attempts);
result.state = state;
end

function fallback = referenceFallback(car,speed_mps,seed,profile, ...
        anyEvaluated,reducedFound,control)
fallback = struct('attempted',false,'status',"not_attempted", ...
    'reason',"",'exitflag',NaN,'state',NaN(1,9), ...
    'long_accel_mps2',NaN);
if reducedFound || shouldCancel(control)
    return
end
fallback.attempted = true;
try
    seedState = [];
    if isnumeric(seed) && numel(seed) == 9 && all(isfinite(seed(:)))
        seedState = seed;
    elseif isstruct(seed) && isfield(seed,'state')
        seedState = seed.state;
    end
    [~,longAccel,state,diagnostics] = max_long_accel( ...
        speed_mps,car,seedState,profile.solverOptions);
    fallback.exitflag = diagnostics.exitflag;
    fallback.state = state;
    fallback.long_accel_mps2 = longAccel;
    if diagnostics.exitflag > 0 && diagnostics.max_equality_residual <= ...
            profile.solverOptions.constraintTolerance
        fallback.status = "reference_feasible";
        fallback.reason = "reduced solver did not accept a candidate";
    else
        fallback.status = "reference_fallback";
        fallback.reason = "reference solver also failed feasibility checks";
    end
catch ME
    fallback.status = "reference_fallback";
    fallback.reason = string(ME.message);
end
end
