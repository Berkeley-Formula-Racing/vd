function [ok,issues] = validateStudy(study)
%VALIDATESTUDY Check the versioned study shape and canonical SI columns.

issues = strings(0,1);
if ~isstruct(study) || ~isscalar(study)
    issues(end+1) = "study must be a scalar struct";
    ok = false;
    return
end
requiredStudy = {'schemaVersion','created','appVersion','cases','runs','displayUnits'};
for i = 1:numel(requiredStudy)
    if ~isfield(study,requiredStudy{i})
        issues(end+1) = "missing study field: " + requiredStudy{i};
    end
end
if isfield(study,'schemaVersion') && ~isequal(study.schemaVersion,1)
    issues(end+1) = "unsupported study schema version";
end
if isfield(study,'setupSpecifications')
    if ~isstruct(study.setupSpecifications)
        issues(end+1) = "study.setupSpecifications must be a struct array";
    elseif ~isempty(study.setupSpecifications)
        setupIds = strings(numel(study.setupSpecifications),1);
        for i = 1:numel(study.setupSpecifications)
            if ~isfield(study.setupSpecifications(i),'id') || ...
                    strlength(strtrim(string( ...
                    study.setupSpecifications(i).id))) == 0
                issues(end+1) = sprintf( ...
                    'setup specification %d is missing id',i);
            else
                setupIds(i) = string(study.setupSpecifications(i).id);
            end
        end
        nonempty = strlength(strtrim(setupIds)) > 0;
        if any(nonempty) && numel(unique(setupIds(nonempty))) ~= nnz(nonempty)
            issues(end+1) = "duplicate setup specification IDs";
        end
    end
end
if isfield(study,'cases')
    if ~isstruct(study.cases)
        issues(end+1) = "study.cases must be a struct array";
    else
        ids = strings(numel(study.cases),1);
        for i = 1:numel(study.cases)
            if ~isfield(study.cases(i),'id')
                issues(end+1) = sprintf('missing case field: id (case %d)',i);
            else
                ids(i) = string(study.cases(i).id);
            end
        end
        nonempty = strlength(strtrim(ids)) > 0;
        if any(nonempty) && numel(unique(ids(nonempty))) ~= nnz(nonempty)
            issues(end+1) = "duplicate case IDs";
        end
    end
end

if isfield(study,'runs')
    if ~isstruct(study.runs)
        issues(end+1) = "study.runs must be a struct array";
    else
        for i = 1:numel(study.runs)
            issues = validateRun(study.runs(i),i,issues);
        end
    end
end

issues = unique(issues,'stable');
ok = isempty(issues);
end

function issues = validateRun(run,index,issues)
requiredRun = {'schemaVersion','caseId','type','mode','settings', ...
    'perSpeed','points','runMeta','status','raw'};
for i = 1:numel(requiredRun)
    if ~isfield(run,requiredRun{i})
        issues(end+1) = sprintf('missing run field: %s (run %d)', ...
            requiredRun{i},index);
    end
end
if isfield(run,'schemaVersion') && ~isequal(run.schemaVersion,1)
    issues(end+1) = sprintf('unsupported run schema version (run %d)',index);
end
if isfield(run,'type')
    type = lower(string(run.type));
    if ~isscalar(type) || ~any(type == ["lateral","longitudinal"])
        issues(end+1) = "unsupported run type: " + string(run.type);
    end
end
if isfield(run,'status')
    status = lower(string(run.status));
    allowed = ["pending","planned","running","converged","near_feasible", ...
        "completed","complete","failed","solver_failed","infeasible", ...
        "cancelled","partial","warning","unknown","missing","invalid"];
    if ~isscalar(status) || ~any(status == allowed)
        issues(end+1) = "invalid run status: " + string(run.status);
    end
end
if isfield(run,'perSpeed')
    if ~istable(run.perSpeed)
        issues(end+1) = sprintf('run.perSpeed must be a table (run %d)',index);
    else
        catalog = schemaCatalog();
        issues = validateColumns(run.perSpeed.Properties.VariableNames, ...
            catalog.perSpeed,index,issues);
        if ~any(strcmp(run.perSpeed.Properties.VariableNames,'speed_mps'))
            issues(end+1) = sprintf('missing canonical perSpeed column: speed_mps (run %d)',index);
        end
    end
end
if isfield(run,'points')
    if ~istable(run.points)
        issues(end+1) = sprintf('run.points must be a table (run %d)',index);
    else
        catalog = schemaCatalog();
        issues = validateColumns(run.points.Properties.VariableNames, ...
            catalog.points,index,issues);
    end
end
if isfield(run,'runMeta') && (~isstruct(run.runMeta) || ~isscalar(run.runMeta))
    issues(end+1) = sprintf('run.runMeta must be a scalar struct (run %d)',index);
end
end

function issues = validateColumns(names,catalog,index,issues)
names = string(names);
unknown = setdiff(names,string(catalog),'stable');
unknown = setdiff(unknown,stablePerSpeedColumns(),'stable');
for i = 1:numel(unknown)
    issues(end+1) = sprintf('non-SI canonical column: %s (run %d)', ...
        unknown(i),index);
end
end

function names = stablePerSpeedColumns()
names = ["speedIndex","origin","passIndex","refinementReason", ...
    "solver_status","solver_reason"];
end

function catalog = schemaCatalog()
template = rampSpeed.makeRun("lateral","coast",struct(),struct());
catalog.perSpeed = string(template.perSpeed.Properties.VariableNames);
catalog.points = string(template.points.Properties.VariableNames);
end
