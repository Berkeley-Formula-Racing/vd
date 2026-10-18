function report = benchmarkQssPerformance(car,options)
%BENCHMARKQSSPERFORMANCE Repeatable QSS timing and validity snapshot.
%   REPORT = BENCHMARKQSSPERFORMANCE(CAR,OPTIONS) records direct equation
%   timing, profiler call counts for the ride-height/aero/camber path, one
%   lateral optimizer result, and optionally a fast g-g grid.  It returns data
%   rather than writing files so callers can compare identical cars/options.

if nargin < 1 || isempty(car), car = carConfigBaseline(); end
if nargin < 2 || isempty(options), options = struct(); end
if ~isstruct(options) || ~isscalar(options)
    error('benchmarkQssPerformance:badOptions','options must be a scalar struct.');
end
sampleCount = optionOr(options,'equationSamples',20);
warmupSamples = optionOr(options,'warmupSamples',2);
timingRepeats = optionOr(options,'timingRepeats',3);
runGG = logical(optionOr(options,'runGG',false));
validateattributes(sampleCount,{'numeric'},{'scalar','integer','positive'});
validateattributes(warmupSamples,{'numeric'},{'scalar','integer','nonnegative'});
validateattributes(timingRepeats,{'numeric'},{'scalar','integer','positive'});

P = [4,0.15,20,0.2,0.45,0,0,0.02,0.02];
for k = 1:warmupSamples
    P(4) = 0.2 + 0.001*k;
    car.equations(P);
end

equationSeconds = zeros(1,timingRepeats);
for repeat = 1:timingRepeats
    timer = tic;
    runEquationSamples(car,P,sampleCount);
    equationSeconds(repeat) = toc(timer);
end

% Keep profiler overhead out of the headline timing while retaining one
% repeatable call-count snapshot for the measured hot path.
profile('off');
profile('clear');
profile('on','-history');
timer = tic;
runEquationSamples(car,P,sampleCount);
profiledEquationSeconds = toc(timer);
profile('off');
profileInfo = profile('info');

report = struct();
report.timestamp = datetime('now');
report.equations = struct('samples',sampleCount,'warmupSamples',warmupSamples, ...
    'timingRepeats',timingRepeats,'seconds',median(equationSeconds), ...
    'secondsPerSample',median(equationSeconds)/sampleCount, ...
    'secondsByRepeat',equationSeconds,'minimumSeconds',min(equationSeconds), ...
    'maximumSeconds',max(equationSeconds), ...
    'profiledSeconds',profiledEquationSeconds);
report.camber = struct('preloaded',isfield(car.camberKinematics,'camberModelData'), ...
    'modelDataFields',{{}},'modelDataLoadCount',0, ...
    'ratioDataPreloaded',~isempty(car.tire.camberInterpolator), ...
    'evaluationCalls',0,'fallbackHotPath',false);
if report.camber.preloaded
    report.camber.modelDataFields = fieldnames(car.camberKinematics.camberModelData);
    report.camber.modelDataLoadCount = 1;
else
    report.camber.fallbackHotPath = true;
end
report.profile = summarizeProfile(profileInfo);
report.camber.evaluationCalls = report.profile.camberEvaluation.calls;
report.aero = struct('rideHeightResidualCalls', ...
    report.profile.aeroRideResidual.calls, ...
    'rideHeightResidualBatchCalls', ...
    report.profile.aeroRideResidualBatch.calls, ...
    'mapEvaluationCalls',report.profile.aeroMapEvaluateNumeric.calls, ...
    'rideHeightSolverCalls',report.profile.solveRideHeightAero.calls);

[~,~,~,acceptedState,optimizer] = max_lat_accel(20,car);
report.optimizer = optimizer;
report.optimizer.acceptedState = acceptedState;
report.optimizer.accepted = acceptedStateEvidence(car,acceptedState);

if runGG
    gridOptions = struct('fastScreening',true,'continuation',true);
    if isfield(options,'ggOptions') && ~isempty(options.ggOptions)
        gridOptions = options.ggOptions;
    end
    timer = tic;
    paramArr = gg2(car,0,gridOptions);
    ggCar = makeGG(paramArr,car);
    report.gg = struct('requestedSize',size(paramArr),'seconds',toc(timer), ...
        'lateralValidCount',nnz(ggCar.ggMask.lateral), ...
        'accelerationValidCount',nnz(ggCar.ggMask.acceleration), ...
        'brakingValidCount',nnz(ggCar.ggMask.braking), ...
        'pointCount',size(ggCar.ggPoints,1), ...
        'points',ggCar.ggPoints,'missingVelocityCount', ...
        nnz(~rowMask(ggCar.ggMask,'lateral')));
else
    report.gg = struct('requestedSize',[0 0], ...
        'lateralValidCount',0,'accelerationValidCount',0, ...
        'pointCount',0,'points',zeros(0,3),'seconds',0, ...
        'missingVelocityCount',0);
end
end

function runEquationSamples(car,P,sampleCount)
for k = 1:sampleCount
    P(4) = 0.2 + 0.001*k;
    car.equations(P);
end
end

function evidence = acceptedStateEvidence(car,P)
evidence = struct('state',P,'constraints',[],'equalities',[], ...
    'residual',Inf,'aeroConverged',false,'aeroOutsideMap',true, ...
    'aeroCoverageValid',false,'aeroResidualIn',Inf, ...
    'minFzvirtual',-Inf);
if isempty(P)
    return
end
[engineRpm,beta,latAccel,~,yawAccel,wheelAccel,~,~,Fzvirtual,~,~,~,~,~,~,ssInfo,~] = ...
    car.equations(P);
c = car.qssInequalities(engineRpm,beta,Fzvirtual,ssInfo);
ceq = car.qssFiniteResidual([latAccel,yawAccel,wheelAccel(1:4)]);
evidence.constraints = c;
evidence.equalities = ceq;
evidence.residual = qssConstraintResidual(c,ceq);
evidence.aeroConverged = diagnosticValue(ssInfo,'aero_converged',false);
evidence.aeroOutsideMap = diagnosticValue(ssInfo,'aero_outside_map',true);
evidence.aeroCoverageValid = diagnosticValue(ssInfo,'aero_coverage_valid',false);
evidence.aeroResidualIn = diagnosticValue(ssInfo,'aero_residual_in',Inf);
evidence.minFzvirtual = min(Fzvirtual);
end

function mask = rowMask(ggMask,name)
mask = false(0,1);
if ~isstruct(ggMask) || ~isfield(ggMask,name) || ...
        ~isfield(ggMask,'velocity') || isempty(ggMask.velocity)
    return
end
velocity = ggMask.velocity(:);
values = logical(ggMask.(name)(:));
uniqueVelocity = unique(round(double(velocity),6));
mask = false(size(uniqueVelocity));
for k = 1:numel(uniqueVelocity)
    mask(k) = any(values(abs(round(double(velocity),6)-uniqueVelocity(k)) <= 1e-9));
end
end

function summary = summarizeProfile(info)
summary = struct('aeroRideResidual',emptyEntry(), ...
    'aeroRideResidualBatch',emptyEntry(), ...
    'aeroMapEvaluateNumeric',emptyEntry(), ...
    'camberEvaluation',emptyEntry(), ...
    'solveRideHeightAero',emptyEntry(), ...
    'tireForce',emptyEntry(), ...
    'tireFx',emptyEntry(), ...
    'tireFy',emptyEntry());
names = fieldnames(summary);
for k = 1:numel(names)
    key = names{k};
    entry = emptyEntry();
    if isfield(info,'FunctionTable')
        table = info.FunctionTable;
        for j = 1:numel(table)
            functionName = string(table(j).FunctionName);
            lowerName = lower(functionName);
            matches = (strcmp(key,'aeroRideResidual') && ...
                contains(lowerName,'aerorideresidual') && ...
                ~contains(lowerName,'aerorideresidualbatch')) || ...
                (strcmp(key,'aeroRideResidualBatch') && ...
                contains(lowerName,'aerorideresidualbatch')) || ...
                (strcmp(key,'aeroMapEvaluateNumeric') && ...
                contains(lowerName,'aeromap') && contains(lowerName,'evaluatenumeric')) || ...
                (strcmp(key,'camberEvaluation') && ...
                contains(lowerName,'camber_evaluation')) || ...
                (strcmp(key,'solveRideHeightAero') && ...
                contains(lowerName,'solverideheightaero')) || ...
                (strcmp(key,'tireForce') && contains(lowerName,'tireforce')) || ...
                (strcmp(key,'tireFx') && contains(lowerName,'tire2') && ...
                contains(lowerName,'f_x')) || ...
                (strcmp(key,'tireFy') && contains(lowerName,'tire2') && ...
                contains(lowerName,'f_y'));
            if matches
                entry.calls = entry.calls + table(j).NumCalls;
                entry.seconds = entry.seconds + table(j).TotalTime;
            end
        end
    end
    summary.(key) = entry;
end
end

function value = diagnosticValue(info,name,defaultValue)
if isstruct(info) && isfield(info,name) && ~isempty(info.(name))
    value = info.(name);
else
    value = defaultValue;
end
end

function entry = emptyEntry()
entry = struct('calls',0,'seconds',0);
end

function value = optionOr(options,name,defaultValue)
if isfield(options,name) && ~isempty(options.(name))
    value = options.(name);
else
    value = defaultValue;
end
end
