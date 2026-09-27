function run = assembleRun(tasks,results,request,execution)
%ASSEMBLERUN Assemble one canonical raw CaseRun row for every planned task.
if nargin < 3 || isempty(request), request = struct(); end
if nargin < 4 || isempty(execution), execution = struct(); end
if ~istable(tasks) || ~iscell(results) || numel(results) ~= height(tasks)
    error("rampSpeed:invalidExecutionResults", ...
        "tasks and results must have matching table/cell dimensions.");
end
required = {'speedIndex','speed_mps','origin','passIndex'};
if ~all(ismember(required,tasks.Properties.VariableNames))
    error("rampSpeed:invalidSpeedTasks", ...
        "tasks must include speedIndex, speed_mps, origin, and passIndex.");
end
indices = double(tasks.speedIndex(:));
if numel(unique(indices)) ~= numel(indices)
    error("rampSpeed:duplicateSpeedIndex", ...
        "Stable speedIndex values must be unique in a CaseRun.");
end

n = height(tasks);
origins = string(tasks.origin(:));
passes = double(tasks.passIndex(:));
if ismember('refinementReason',tasks.Properties.VariableNames)
    refinementReasons = string(tasks.refinementReason(:));
else
    refinementReasons = strings(n,1);
end
statuses = strings(n,1);
valid = false(n,1);
reasons = strings(n,1);
metricNames = strings(0,1);
for i = 1:n
    result = results{i};
    if isempty(result)
        continue
    end
    statuses(i) = string(result.status);
    valid(i) = ismember(statuses(i),["converged","near_feasible"]);
    reasons(i) = resultReason(result);
    if isfield(result,'metrics') && isstruct(result.metrics) && isscalar(result.metrics)
        names = string(fieldnames(result.metrics));
        metricNames = [metricNames; names(~ismember(names,metricNames))]; %#ok<AGROW>
    end
end
statuses(strlength(statuses)==0) = "planned";
reasons(strlength(reasons)==0 & ~valid) = ...
    "no result returned for planned speed";

run = struct();
run.type = getText(request,"rampType","longitudinal");
run.perSpeed = table(indices,double(tasks.speed_mps(:)),origins,passes, ...
    refinementReasons,valid,statuses,reasons, ...
    'VariableNames',{'speedIndex','speed_mps','origin','passIndex', ...
    'refinementReason','valid','status','reason'});
for name = metricNames(:).'
    field = matlab.lang.makeValidName(char(name));
    if ismember(field,run.perSpeed.Properties.VariableNames)
        continue
    end
    values = NaN(n,1);
    for i = 1:n
        if isempty(results{i}) || ~isfield(results{i},'metrics') || ...
                ~isfield(results{i}.metrics,char(name))
            continue
        end
        value = results{i}.metrics.(char(name));
        if isnumeric(value) && isscalar(value) && isfinite(value)
            values(i) = double(value);
        end
    end
    run.perSpeed.(field) = values;
end
if string(run.type) == "lateral"
    lateralPoints = collectLateralPoints(tasks,results);
    if isempty(lateralPoints)
        run.points = run.perSpeed;
    else
        run.points = lateralPoints;
    end
else
    run.points = run.perSpeed;
end
run.raw = struct("perSpeed",run.perSpeed,"points",run.points, ...
    "diagnostics",makeDiagnostics(tasks,results), ...
    "speedErrors",makeSpeedErrors(tasks,results), ...
    "status",aggregateStatus(run.perSpeed,execution));
% Keep diagnostics and speed errors available at the CaseRun root for legacy
% consumers while retaining the nested raw payload used by the new contract.
run.diagnostics = run.raw.diagnostics;
run.speedErrors = run.raw.speedErrors;
run.settings = getStructField(request,"settings",struct());
run.status = run.raw.status;
run.runMeta = struct("source","rampSpeed.assembleRun", ...
    "status",run.status,"request",request, ...
    "completedTasks",sum(~ismember(run.perSpeed.status,["planned","running"])), ...
    "plannedTasks",n);
end

function diagnostics = makeDiagnostics(tasks,results)
template = struct("speedIndex",NaN,"speed_mps",NaN,"status","", ...
    "state",NaN(1,9),"metrics",struct(),"diagnostics",struct(), ...
    "attempts",struct([]),"reason","","success",false, ...
    "exitflag",NaN,"c",NaN,"ceq",NaN, ...
    "max_equality_residual",NaN,"max_inequality_violation",NaN, ...
    "error_identifier","","error_message","");
diagnostics = repmat(template,height(tasks),1);
for i = 1:height(tasks)
    diagnostics(i).speedIndex = tasks.speedIndex(i);
    diagnostics(i).speed_mps = tasks.speed_mps(i);
    result = results{i};
    if isempty(result), continue, end
    diagnostics(i).status = string(result.status);
    if isfield(result,'state') && isnumeric(result.state) && ...
            numel(result.state)==9
        diagnostics(i).state = double(result.state(:).');
    end
    if isfield(result,'metrics'), diagnostics(i).metrics = result.metrics; end
    if isfield(result,'diagnostics'), diagnostics(i).diagnostics = result.diagnostics; end
    if isfield(result,'attempts'), diagnostics(i).attempts = result.attempts; end
    diagnostics(i).reason = resultReason(result);
    diagnostics(i).success = ismember(diagnostics(i).status, ...
        ["converged","near_feasible"]);
    d = diagnostics(i).diagnostics;
    if isfield(d,'maxEqualityResidual') && isfinite(d.maxEqualityResidual)
        diagnostics(i).max_equality_residual = d.maxEqualityResidual;
    end
    if isfield(d,'maxInequalityViolation') && isfinite(d.maxInequalityViolation)
        diagnostics(i).max_inequality_violation = d.maxInequalityViolation;
    end
    if isfield(d,'errorIdentifier')
        diagnostics(i).error_identifier = string(d.errorIdentifier);
    end
    attempts = diagnostics(i).attempts;
    if isempty(attempts) && isfield(d,'gearAttempts')
        attempts = d.gearAttempts;
    end
    if ~isempty(attempts) && isfield(attempts,'exitflag')
        validAttempt = find(isfinite([attempts.exitflag]),1,'last');
        if ~isempty(validAttempt)
            selected = attempts(validAttempt);
            if isfield(selected,'state') && isnumeric(selected.state) && ...
                    numel(selected.state) == 9 && ...
                    all(isfinite(selected.state(:)))
                diagnostics(i).state = double(selected.state(:).');
            end
            diagnostics(i).exitflag = selected.exitflag;
            if isfield(selected,'maxEqualityResidual') && ...
                    isfinite(selected.maxEqualityResidual)
                diagnostics(i).max_equality_residual = selected.maxEqualityResidual;
            end
            if isfield(selected,'maxInequalityViolation') && ...
                    isfinite(selected.maxInequalityViolation)
                diagnostics(i).max_inequality_violation = ...
                    selected.maxInequalityViolation;
            end
        end
    end
    if isfinite(diagnostics(i).max_inequality_violation)
        diagnostics(i).c = diagnostics(i).max_inequality_violation;
    end
    if isfinite(diagnostics(i).max_equality_residual)
        diagnostics(i).ceq = diagnostics(i).max_equality_residual;
    end
    if diagnostics(i).status == "solver_failed" && ...
            isfinite(diagnostics(i).exitflag)
        diagnostics(i).error_identifier = "rampSpeed:nonconverged";
        diagnostics(i).reason = nonconvergedReason(diagnostics(i));
        diagnostics(i).error_message = diagnostics(i).reason;
    elseif strlength(diagnostics(i).error_identifier) == 0 && ...
            diagnostics(i).status == "solver_failed"
        diagnostics(i).error_identifier = "rampSpeed:solver_failed";
    end
    if strlength(diagnostics(i).error_message) == 0
        diagnostics(i).error_message = diagnostics(i).reason;
    end
end
end

function points = collectLateralPoints(tasks,results)
parts = cell(0,1);
for i = 1:height(tasks)
    result = results{i};
    if isempty(result) || ~isfield(result,"diagnostics") || ...
            ~isstruct(result.diagnostics) || ...
            ~isfield(result.diagnostics,"raw")
        continue
    end
    raw = result.diagnostics.raw;
    if ~isstruct(raw) || ~isfield(raw,"points") || ~istable(raw.points) || ...
            height(raw.points) == 0
        continue
    end
    part = raw.points;
    if ismember("speed_index",part.Properties.VariableNames)
        part.speed_index(:) = tasks.speedIndex(i);
    else
        part.speed_index = repmat(tasks.speedIndex(i),height(part),1);
    end
    if ~ismember("speed_mps",part.Properties.VariableNames)
        part.speed_mps = repmat(tasks.speed_mps(i),height(part),1);
    end
    parts{end+1,1} = part; %#ok<AGROW>
end
if isempty(parts)
    points = table();
else
    points = vertcat(parts{:});
end
end

function errors = makeSpeedErrors(tasks,results)
template = struct("speed_index",NaN,"speed_mps",NaN, ...
    "identifier","","message","","stack",[]);
errors = repmat(template,0,1);
for i = 1:height(tasks)
    result = results{i};
    if isempty(result) || ismember(string(result.status),["converged","near_feasible"])
        continue
    end
    entry = template;
    entry.speed_index = tasks.speedIndex(i);
    entry.speed_mps = tasks.speed_mps(i);
    entry.identifier = "rampSpeed:"+string(result.status);
    entry.message = resultReason(result);
    if isfield(result,'diagnostics') && isfield(result.diagnostics,'errorIdentifier')
        entry.identifier = string(result.diagnostics.errorIdentifier);
    end
    attempts = resultAttempts(result);
    if string(result.status) == "solver_failed" && ~isempty(attempts) && ...
            isfield(attempts,'exitflag') && any(isfinite([attempts.exitflag]))
        entry.identifier = "rampSpeed:nonconverged";
        entry.message = nonconvergedReasonFromResult(result,attempts);
    end
    errors(end+1,1) = entry; %#ok<AGROW>
end
end

function reason = resultReason(result)
reason = "";
if isfield(result,'diagnostics') && isstruct(result.diagnostics) && ...
        isfield(result.diagnostics,'reason') && ...
        ~isempty(result.diagnostics.reason)
    reason = string(result.diagnostics.reason);
elseif isfield(result,'status')
    reason = string(result.status);
end
reason = reason(1);
if isfield(result,'diagnostics') && isstruct(result.diagnostics) && ...
        isfield(result.diagnostics,'errorIdentifier') && ...
        strlength(string(result.diagnostics.errorIdentifier)) > 0 && ...
        ~ismember(string(result.status),["converged","near_feasible"])
    reason = string(result.diagnostics.errorIdentifier) + ": " + reason;
end
if string(result.status) == "solver_failed"
    attempts = resultAttempts(result);
    if ~isempty(attempts) && isfield(attempts,'exitflag') && ...
            any(isfinite([attempts.exitflag]))
        reason = nonconvergedReasonFromResult(result,attempts);
    end
end
end

function attempts = resultAttempts(result)
attempts = struct([]);
if isfield(result,'attempts') && ~isempty(result.attempts)
    attempts = result.attempts;
elseif isfield(result,'diagnostics') && isstruct(result.diagnostics) && ...
        isfield(result.diagnostics,'gearAttempts')
    attempts = result.diagnostics.gearAttempts;
end
end

function reason = nonconvergedReason(diagnostic)
reason = "nonconverged solver result: exitflag=" + ...
    string(diagnostic.exitflag);
end

function reason = nonconvergedReasonFromResult(result,attempts)
index = find(isfinite([attempts.exitflag]),1,'last');
reason = "nonconverged solver result: exitflag=" + ...
    string(attempts(index).exitflag);
if isfield(attempts,'exitMessage') && ...
        strlength(string(attempts(index).exitMessage)) > 0
    reason = reason + ": " + string(attempts(index).exitMessage);
elseif isfield(result,'diagnostics') && isfield(result.diagnostics,'reason') && ...
        strlength(string(result.diagnostics.reason)) > 0
    reason = reason + ": " + string(result.diagnostics.reason);
end
end

function status = aggregateStatus(perSpeed,execution)
if isfield(execution,'status') && string(execution.status) == "cancelled"
    status = "cancelled";
elseif isempty(perSpeed)
    status = "failed";
elseif all(perSpeed.valid)
    status = "completed";
elseif any(perSpeed.valid)
    status = "partial";
else
    status = "failed";
end
end

function value = getText(s,name,default)
value = string(default);
if isstruct(s) && isfield(s,name) && ~isempty(s.(name))
    value = string(s.(name));
end
value = value(1);
end

function value = getStructField(s,name,default)
value = default;
if isstruct(s) && isfield(s,name) && isstruct(s.(name))
    value = s.(name);
end
end
