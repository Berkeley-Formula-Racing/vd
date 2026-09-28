function result = solveContinuousEnvelopePoint(car,speed_mps,seed,profile,control)
%SOLVECONTINUOUSENVELOPEPOINT Solve one pure-longitudinal point.
%
% The powertrain ratio is compiled once before fmincon. Throttle and one
% symmetric rear slip ratio remain the only independent solver variables.

if nargin < 3, seed = []; end
if nargin < 4 || isempty(profile), profile = "accurate"; end
if nargin < 5 || isempty(control), control = struct(); end
profile = rampSpeed.resolveSolverProfile(profile);
task = struct('speedIndex',1,'speed_mps',double(speed_mps), ...
    'origin',"requested",'passIndex',1);

if shouldCancel(control)
    result = rampSpeed.makeSpeedResult(task,"cancelled",struct(), ...
        struct('powertrainModel',"continuousEnvelope", ...
        'continuousEnvelope',struct(),'gearAttempts',struct([]), ...
        'evaluationCount',0,'reducedState',[NaN NaN], ...
        'envelopeSource',"not_compiled", ...
        'reason',"cancelled before first evaluation"),struct([]));
    result.state = NaN(1,9);
    return
end

envelopeSource = "pointCompile";
try
    [envelope,cacheHit] = rampSpeed.getContinuousEnvelope( ...
        getRampModel(control),speed_mps,car);
    if cacheHit
        envelopeSource = "cachedRampModel";
    else
        envelope = rampSpeed.buildContinuousEnvelope(car,speed_mps);
    end
catch ME
    status = "solver_failed";
    if startsWith(string(ME.identifier),"rampSpeed:invalid")
        status = "infeasible";
    end
    diagnostics = struct('powertrainModel',"continuousEnvelope", ...
        'continuousEnvelope',struct(),'gearAttempts',struct([]), ...
        'evaluationCount',0,'reducedState',[NaN NaN], ...
        'envelopeSource',envelopeSource, ...
        'reason',string(ME.message),'errorIdentifier',string(ME.identifier));
    result = rampSpeed.makeSpeedResult(task,status,struct(),diagnostics,struct([]));
    result.state = NaN(1,9);
    return
end

attempt = attemptTemplate(envelope);
evaluationCount = 0;
best = struct('found',false,'state',NaN(1,9),'data',struct(), ...
    'status',"infeasible");

if shouldCancel(control)
    result = cancelledResult(task,attempt,evaluationCount,envelope, ...
        "cancelled before optimization",best.state,envelopeSource);
    return
end

try
    [state,data,exitflag,exitMessage,solverOutput] = solveEnvelope();
    attempt.exitflag = exitflag;
    attempt.exitMessage = string(exitMessage);
    attempt.solverOutput = solverOutput;
    attempt.engineRpm = data.engineRpm;
    attempt.maxInequalityViolation = data.maxInequalityViolation;
    attempt.maxEqualityResidual = data.maxEqualityResidual;
    attempt.long_accel_mps2 = data.long_accel_mps2;
    attempt.state = state;
    attempt.feasible = isAcceptable(exitflag,data,profile);
    if attempt.feasible
        if data.maxEqualityResidual <= constraintTolerance(profile) && ...
                data.maxInequalityViolation <= constraintTolerance(profile)
            attempt.status = "converged";
        else
            attempt.status = "near_feasible";
        end
        attempt.reason = "accepted";
        best.found = true;
        best.state = state;
        best.data = data;
        best.status = attempt.status;
    elseif exitflag <= 0
        attempt.status = "solver_failed";
        attempt.reason = "optimizer did not converge to a feasible point";
    else
        attempt.status = "infeasible";
        attempt.reason = "physical or equilibrium residual exceeded tolerance";
    end
catch ME
    if strcmp(ME.identifier,'rampSpeed:cancelled')
        result = cancelledResult(task,attempt,evaluationCount,envelope, ...
            "cancelled during continuous-envelope evaluation",best.state, ...
            envelopeSource);
        return
    end
    attempt.status = "solver_failed";
    attempt.reason = string(ME.message);
    attempt.exitMessage = string(ME.message);
    state = best.state;
end

diagnostics = struct( ...
    'powertrainModel',"continuousEnvelope", ...
    'continuousEnvelope',envelope, ...
    'envelopeSource',envelopeSource, ...
    'gearAttempts',struct([]), ...
    'evaluationCount',evaluationCount, ...
    'referenceFallback',struct('attempted',false,'status',"not_attempted", ...
        'reason',"continuousEnvelope is the primary ramp powertrain model"), ...
    'reducedState',[NaN NaN]);

if best.found
    state = best.state;
    evaluationOptions = struct('continuousRatio',envelope.drivetrainReduction);
    metrics = car.metrics(state,evaluationOptions);
    metrics.aLong_mps2 = best.data.long_accel_mps2;
    metrics.aLong_max_mps2 = best.data.long_accel_mps2;
    metrics.aLat_mps2 = speed_mps*state(5);
    metrics.longitudinal_acceleration_mps2 = best.data.long_accel_mps2;
    metrics.lat_accel_residual = best.data.lat_accel;
    metrics.yaw_accel_residual = best.data.yaw_accel;
    metrics.max_equality_residual = best.data.maxEqualityResidual;
    metrics.max_inequality_violation = best.data.maxInequalityViolation;
    metrics.current_gear = NaN;
    metrics.powertrain_model = "continuousEnvelope";
    metrics.drivetrain_reduction = envelope.drivetrainReduction;
    metrics.continuous_envelope_wheel_force_N = envelope.wheelForce_N;
    metrics.rear_slip = state(8);
    diagnostics.reducedState = [state(2),state(8)];
    diagnostics.selectedEngineRpm = envelope.engineRpm;
    diagnostics.drivetrainReduction = envelope.drivetrainReduction;
    diagnostics.maxEqualityResidual = best.data.maxEqualityResidual;
    diagnostics.maxInequalityViolation = best.data.maxInequalityViolation;
    status = best.status;
else
    metrics = struct();
    if attempt.status == "solver_failed"
        status = "solver_failed";
    else
        status = "infeasible";
    end
    diagnostics.reducedState = [NaN NaN];
    diagnostics.selectedEngineRpm = NaN;
    diagnostics.drivetrainReduction = envelope.drivetrainReduction;
    diagnostics.maxEqualityResidual = attempt.maxEqualityResidual;
    diagnostics.maxInequalityViolation = attempt.maxInequalityViolation;
end

result = rampSpeed.makeSpeedResult(task,status,metrics,diagnostics,attempt);
result.state = state;

    function [stateOut,dataOut,exitflag,exitMessage,solverOutput] = solveEnvelope()
        bounds = [0 1; 0 0.2];
        x0 = initialReducedState(seed,bounds);
        lastX = NaN(1,2);
        lastData = struct();
        rideHeightContext = [];

        function dataLocal = evaluateReduced(z)
            z = double(z(:).');
            if isequal(z,lastX)
                dataLocal = lastData;
                return
            end
            if shouldCancel(control)
                error('rampSpeed:cancelled', ...
                    'cancelled during longitudinal point evaluation');
            end
            evaluationCount = evaluationCount + 1;
            stateLocal = rampSpeed.makePureLongState(car,speed_mps,z(1),z(2));
            options = struct('continuousRatio',envelope.drivetrainReduction);
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
            dataLocal = struct('state',stateLocal,'engineRpm',engineRpm, ...
                'beta',beta,'lat_accel',latAccel, ...
                'long_accel_mps2',longAccel,'yaw_accel',yawAccel, ...
                'wheel_accel',wheelAccel,'omega',omega, ...
                'currentGear',currentGear,'Fzvirtual',Fzvirtual,'Fz',Fz, ...
                'alpha',alpha,'T',T,'Fy',Fy,'gamma',gamma,'Fx',Fx, ...
                'ssInfo',ss,'inequality',inequality,'equality',equality, ...
                'maxInequalityViolation',maxIneq, ...
                'maxEqualityResidual',maxEq);
            lastX = z;
            lastData = dataLocal;
        end

        objective = @objectiveReduced;
        nonlinear = @constraintsReduced;
        options = fminconOptions(profile);
        [z,~,exitflag,solverOutput] = fmincon(objective,x0,[],[],[],[], ...
            bounds(:,1),bounds(:,2),nonlinear,options);
        dataOut = evaluateReduced(z);
        stateOut = dataOut.state;
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
            ceq = mean(dataLocal.wheel_accel(3:4));
            if any(~isfinite(c)) || any(~isfinite(ceq))
                c = ones(numel(c),1)*realmax('double')/100;
                ceq = ones(numel(ceq),1)*realmax('double')/100;
            end
        end
    end
end

function attempt = attemptTemplate(envelope)
attempt = struct('powertrainModel',"continuousEnvelope", ...
    'drivetrainReduction',envelope.drivetrainReduction, ...
    'engineRpm',envelope.engineRpm,'exitflag',NaN,'exitMessage',"", ...
    'solverOutput',struct(),'maxInequalityViolation',NaN, ...
    'maxEqualityResidual',NaN,'long_accel_mps2',NaN, ...
    'state',NaN(1,9),'feasible',false,'status',"not_attempted", ...
    'reason',"");
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

function tf = isAcceptable(exitflag,data,profile)
tol = constraintTolerance(profile);
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

function result = cancelledResult(task,attempt,evaluationCount,envelope,reason, ...
        state,envelopeSource)
if nargin < 7 || isempty(envelopeSource)
    envelopeSource = "pointCompile";
end
diagnostics = struct('powertrainModel',"continuousEnvelope", ...
    'continuousEnvelope',envelope,'gearAttempts',struct([]), ...
    'evaluationCount',evaluationCount,'reducedState',[NaN NaN], ...
    'envelopeSource',string(envelopeSource), ...
    'reason',string(reason));
result = rampSpeed.makeSpeedResult(task,"cancelled",struct(),diagnostics,attempt);
result.state = state;
end

function model = getRampModel(control)
model = struct();
if isstruct(control) && isscalar(control) && isfield(control,'rampModel') && ...
        isstruct(control.rampModel) && isscalar(control.rampModel)
    model = control.rampModel;
end
end
