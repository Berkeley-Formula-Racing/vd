classdef RampSpeedSession < handle
    %RAMPSPEEDSESSION UI-independent Ramp Speed setup and run controller.

    properties (SetAccess=private)
        config = []
        cars = []
        cases = struct.empty(0,1)
        setupSpecifications = struct.empty(0,1)
        designTable = table()
        selectedCaseIds = strings(0,1)
        rampType = "lateral"
        settings = struct()
        options = struct()
        state = "idle"
        stateHistory = "idle"
        progress = struct("completedCases",0,"totalCases",0, ...
            "speedIndex",NaN,"speed_mps",NaN,"status","idle","message","")
        study = struct()
        job = []
        readOnly = false
    end

    methods
        function obj = RampSpeedSession(first,second,third)
            if nargin < 1
                error("rampSpeed:invalidSessionInput", ...
                    "A setup configuration or cars/cases pair is required.");
            end

            isConfig = isstruct(first) && isscalar(first) && ...
                isfield(first,"defaultSetup");
            if isConfig
                obj.config = first;
                if nargin >= 2 && ~isempty(second)
                    options = second;
                else
                    options = struct();
                end
                [obj.cars,obj.cases,obj.designTable] = ...
                    rampSpeed.buildSetupCatalog(obj.config, ...
                    obj.config.defaultSetup);
                obj.setupSpecifications = reshape([obj.cases.setupSpec],[],1);
            else
                if nargin < 2
                    error("rampSpeed:invalidSessionInput", ...
                        "Fake/injected sessions require cars and cases.");
                end
                obj.cars = first;
                obj.cases = normalizeInjectedCases(second);
                obj.setupSpecifications = reshape([obj.cases.setupSpec],[],1);
                obj.designTable = makeSetupTable(obj.setupSpecifications);
                if nargin >= 3
                    options = third;
                else
                    options = struct();
                end
            end

            if isempty(obj.cases)
                obj.selectedCaseIds = strings(0,1);
            else
                obj.selectedCaseIds = string({obj.cases.id}).';
            end
            obj.options = normalizeOptions(options);
            obj.rampType = normalizeRampType(fieldOr( ...
                obj.options,"rampType","lateral"));
            obj.settings = fieldOr(obj.options,"settings",struct());
            obj.study = rampSpeed.makeStudy(fieldOr( ...
                obj.options,"appVersion","dev"));
            obj.study.cases = obj.cases;
            obj.study.setupSpecifications = obj.setupSpecifications;
            obj.study.baselineVersion = baselineVersion( ...
                obj.setupSpecifications);
        end

        function model = viewModel(obj)
            obj.syncJob();
            model = struct();
            model.state = obj.state;
            model.stateHistory = obj.stateHistory;
            model.setupTable = obj.designTable;
            model.setupSpecifications = obj.setupSpecifications;
            model.cases = obj.cases;
            model.selectedCaseIds = obj.selectedCaseIds;
            model.rampType = obj.rampType;
            model.settings = obj.settings;
            model.study = obj.study;
            model.progress = obj.progress;
            model.readOnly = obj.readOnly;
            if isempty(obj.job)
                model.jobState = "idle";
            else
                model.jobState = string(obj.job.state);
            end
        end

        function duplicateSetup(obj,sourceId,newId,newLabel)
            obj.ensureEditable();
            sourceIndex = obj.setupIndex(sourceId);
            source = obj.setupSpecifications(sourceIndex);
            if nargin < 3 || isempty(newId)
                newId = string(source.id) + "-copy";
            end
            if nargin < 4 || isempty(newLabel)
                newLabel = string(source.label) + " copy";
            end
            newId = string(newId);
            if any(string({obj.cases.id}) == newId)
                error("rampSpeed:duplicateSetupId", ...
                    "Setup ID %s already exists.",newId);
            end
            copy = rampSpeed.duplicateSetup(source,newId,newLabel);
            specs = [obj.setupSpecifications(:);copy];
            obj.commitSetupCatalog(specs,sourceIndex);
        end

        function editSetup(obj,setupId,changes)
            obj.ensureEditable();
            index = obj.setupIndex(setupId);
            if obj.isBaseline(index)
                error("rampSpeed:baselineProtected", ...
                    "The baseline setup cannot be edited.");
            end
            if ~isstruct(changes) || ~isscalar(changes)
                error("rampSpeed:invalidSetupChanges", ...
                    "changes must be a scalar struct.");
            end
            current = obj.setupSpecifications(index);
            edited = mergeStruct(current,changes);
            edited.id = current.id;
            edited.isBaseline = false;
            specs = obj.setupSpecifications;
            specs(index) = edited;
            obj.commitSetupCatalog(specs,index);
        end

        function deleteSetup(obj,setupId)
            obj.ensureEditable();
            index = obj.setupIndex(setupId);
            if obj.isBaseline(index)
                error("rampSpeed:baselineProtected", ...
                    "The baseline setup cannot be deleted.");
            end
            if numel(obj.setupSpecifications) <= 1
                error("rampSpeed:lastSetup", ...
                    "The baseline setup is the minimum catalog.");
            end
            specs = obj.setupSpecifications;
            specs(index) = [];
            obj.commitSetupCatalog(specs,max(1,min(index,numel(specs))));
        end

        function selectCases(obj,caseIds)
            if obj.isRunning()
                error("rampSpeed:sessionBusy", ...
                    "Cases cannot be changed while a study is running.");
            end
            ids = string(caseIds(:));
            ids = ids(strlength(strtrim(ids)) > 0);
            available = string({obj.cases.id}).';
            if any(~ismember(ids,available)) || numel(unique(ids)) ~= numel(ids)
                error("rampSpeed:invalidSelection", ...
                    "Selected case IDs must be unique entries in the setup catalog.");
            end
            obj.selectedCaseIds = ids;
        end

        function setRampType(obj,value)
            if obj.isRunning()
                error("rampSpeed:sessionBusy", ...
                    "Ramp type cannot be changed while a study is running.");
            end
            obj.rampType = normalizeRampType(value);
        end

        function job = start(obj,requestOverride)
            if nargin < 2 || isempty(requestOverride)
                requestOverride = struct();
            end
            if obj.isRunning()
                error("rampSpeed:sessionBusy", ...
                    "A Ramp Speed study is already running.");
            end
            if obj.readOnly
                error("rampSpeed:readOnlySession", ...
                    "A read-only study cannot be started.");
            end
            if isempty(obj.selectedCaseIds)
                error("rampSpeed:noCases", ...
                    "At least one setup must be selected.");
            end

            indices = find(ismember(string({obj.cases.id}), ...
                obj.selectedCaseIds));
            if isempty(indices)
                error("rampSpeed:noCases", ...
                    "No selected setup is available.");
            end
            [cars,cases] = obj.selectedInputs(indices);
            request = obj.makeRequest(requestOverride);
            obj.setState("running");
            obj.progress = struct("completedCases",0, ...
                "totalCases",numel(cases),"speedIndex",NaN, ...
                "speed_mps",NaN,"status","running","message","starting");

            userProgress = fieldOr(obj.options,"onProgress",[]);
            userJobCreated = fieldOr(obj.options,"onJobCreated",[]);
            callbacks = struct();
            callbacks.onJobCreated = @(created)obj.jobCreated( ...
                created,userJobCreated);
            callbacks.onProgress = @(event)obj.progressReceived( ...
                event,userProgress);
            if isfield(obj.options,"isCancelled")
                callbacks.isCancelled = obj.options.isCancelled;
            end

            try
                job = rampSpeed.StudyExecutor.start( ...
                    cars,cases,request,callbacks);
            catch ME
                obj.setState("failed");
                rethrow(ME);
            end
            obj.job = job;
            obj.syncJob();
        end

        function cancel(obj)
            if isempty(obj.job)
                return
            end
            rampSpeed.StudyExecutor.cancel(obj.job);
            obj.syncJob();
        end

        function poll(obj)
            obj.syncJob();
        end

        function acceptResult(obj,study)
            if ~isstruct(study) || ~isscalar(study)
                error("rampSpeed:invalidStudy", ...
                    "study must be a scalar struct.");
            end
            obj.study = study;
            obj.readOnly = fieldOr(study,"readOnly",false);
            obj.setState(mapStudyState(fieldOr(study,"status","complete")));
        end

        function clearResults(obj)
            if obj.isRunning()
                error("rampSpeed:sessionBusy", ...
                    "Results cannot be cleared while a study is running.");
            end
            obj.study = rampSpeed.makeStudy(fieldOr( ...
                obj.options,"appVersion","dev"));
            obj.study.cases = obj.cases;
            obj.study.setupSpecifications = obj.setupSpecifications;
            obj.study.baselineVersion = baselineVersion( ...
                obj.setupSpecifications);
            obj.job = [];
            obj.progress = struct("completedCases",0, ...
                "totalCases",0,"speedIndex",NaN,"speed_mps",NaN, ...
                "status","idle","message","");
            obj.setState("idle");
        end

        function save(obj,fileName)
            candidate = obj.study;
            candidate.cases = obj.cases;
            candidate.setupSpecifications = obj.setupSpecifications;
            candidate.baselineVersion = baselineVersion( ...
                obj.setupSpecifications);
            candidate.readOnly = obj.readOnly;
            rampSpeed.saveStudy(fileName,candidate);
        end

        function loaded = load(obj,fileName)
            loaded = rampSpeed.loadStudy(fileName,fieldOr( ...
                obj.options,"appVersion","dev"));
            if isfield(loaded,"setupSpecifications") && ...
                    ~isempty(loaded.setupSpecifications)
                obj.setupSpecifications = loaded.setupSpecifications;
                if ~isempty(obj.config)
                    [obj.cars,obj.cases,obj.designTable] = ...
                        rampSpeed.buildSetupCatalog(obj.config, ...
                        obj.setupSpecifications);
                else
                    obj.cases = normalizeInjectedCases(loaded.cases);
                    obj.designTable = makeSetupTable(obj.setupSpecifications);
                end
                obj.selectedCaseIds = string({obj.cases.id}).';
            end
            obj.study = loaded;
            obj.readOnly = fieldOr(loaded,"readOnly",false);
            obj.job = [];
            obj.setState(mapStudyState(fieldOr( ...
                loaded,"status","complete")));
        end

        function result = export(obj,outputDirectory,options)
            if nargin < 2
                outputDirectory = pwd;
            end
            if nargin < 3
                options = struct();
            end
            candidate = obj.study;
            candidate.cases = obj.cases;
            candidate.setupSpecifications = obj.setupSpecifications;
            candidate.baselineVersion = baselineVersion( ...
                obj.setupSpecifications);
            result = rampSpeed.exportStudy(candidate, ...
                outputDirectory,options);
        end
    end

    methods (Access=private)
        function value = isRunning(obj)
            value = string(obj.state) == "running";
            if ~isempty(obj.job)
                value = value || any(string(obj.job.state) == ...
                    ["running","queued"]);
            end
        end

        function ensureEditable(obj)
            if obj.isRunning()
                error("rampSpeed:sessionBusy", ...
                    "Setup edits are locked while a study is running.");
            end
            if obj.readOnly
                error("rampSpeed:readOnlySession", ...
                    "This study is read-only.");
            end
        end

        function index = setupIndex(obj,setupId)
            id = string(setupId);
            index = find(string({obj.cases.id}) == id,1,"first");
            if isempty(index)
                error("rampSpeed:unknownSetup", ...
                    "Unknown setup ID %s.",id);
            end
        end

        function value = isBaseline(obj,index)
            value = false;
            if isfield(obj.cases,"isBaseline")
                value = logical(obj.cases(index).isBaseline);
            end
            if ~value && isfield(obj.setupSpecifications,"isBaseline")
                value = logical(obj.setupSpecifications(index).isBaseline);
            end
        end

        function commitSetupCatalog(obj,specs,anchorIndex)
            priorSelection = obj.selectedCaseIds;
            if ~isempty(obj.config)
                [cars,cases,design] = rampSpeed.buildSetupCatalog( ...
                    obj.config,specs);
            else
                [cars,cases,design] = obj.rebuildInjectedCatalog(specs);
            end
            obj.cars = cars;
            obj.cases = cases;
            obj.designTable = design;
            obj.setupSpecifications = reshape([obj.cases.setupSpec],[],1);
            available = string({obj.cases.id}).';
            keep = priorSelection(ismember(priorSelection,available));
            if isempty(keep)
                anchorIndex = max(1,min(anchorIndex,numel(available)));
                keep = available(anchorIndex);
            end
            obj.selectedCaseIds = keep;
            obj.study.cases = obj.cases;
            obj.study.setupSpecifications = obj.setupSpecifications;
            obj.study.baselineVersion = baselineVersion( ...
                obj.setupSpecifications);
        end

        function [cars,cases,design] = rebuildInjectedCatalog(obj,specs)
            count = numel(specs);
            cases = repmat(emptySessionCase(),count,1);
            for i = 1:count
                oldIndex = find(string({obj.cases.id}) == ...
                    string(specs(i).id),1,"first");
                if isempty(oldIndex)
                    oldIndex = min(i,numel(obj.cases));
                end
                if oldIndex > 0 && oldIndex <= numel(obj.cases)
                    cases(i) = obj.cases(oldIndex);
                end
                cases(i).id = string(specs(i).id);
                cases(i).label = string(specs(i).label);
                cases(i).source = fieldOr(specs(i),"source","session");
                cases(i).designRow = i;
                cases(i).sourceIndex = i;
                cases(i).carRole = "auto";
                cases(i).carColumn = 1;
                cases(i).setupSpec = specs(i);
                cases(i).derived = fieldOr(specs(i),"derived",struct());
                cases(i).isBaseline = fieldOr(specs(i),"isBaseline",i == 1);
            end
            if iscell(obj.cars)
                oldCars = obj.cars(:);
                cars = cell(count,1);
                for i = 1:count
                    oldIndex = find(string({obj.cases.id}) == ...
                        string(specs(i).id),1,"first");
                    if ~isempty(oldIndex) && oldIndex <= numel(oldCars)
                        cars{i} = oldCars{oldIndex};
                    elseif i <= numel(oldCars)
                        cars{i} = oldCars{i};
                    else
                        cars{i} = oldCars{1};
                    end
                end
            else
                cars = obj.cars;
            end
            design = makeSetupTable(specs);
        end

        function [cars,cases] = selectedInputs(obj,indices)
            cases = obj.cases(indices);
            if iscell(obj.cars)
                cars = obj.cars(indices);
            else
                cars = obj.cars(indices);
            end
            for i = 1:numel(cases)
                cases(i).designRow = i;
                cases(i).sourceIndex = i;
                cases(i).carColumn = 1;
            end
        end

        function request = makeRequest(obj,override)
            request = struct( ...
                "rampType",obj.rampType, ...
                "settings",obj.settings, ...
                "parallelRequested",false, ...
                "numWorkers",0, ...
                "checkpointPath","", ...
                "appVersion",fieldOr(obj.options,"appVersion","dev"), ...
                "runCaseFcn",@rampSpeed.runCase);
            optionFields = ["parallelRequested","numWorkers", ...
                "checkpointPath","appVersion","runCaseFcn"];
            for name = optionFields
                if isfield(obj.options,char(name)) && ...
                        ~isempty(obj.options.(char(name)))
                    request.(char(name)) = obj.options.(char(name));
                end
            end
            if isfield(obj.options,"request") && ...
                    isstruct(obj.options.request)
                request = mergeStruct(request,obj.options.request);
            end
            if isstruct(override)
                request = mergeStruct(request,override);
            end
            request.rampType = obj.rampType;
            if ~isfield(request,"settings") || isempty(request.settings)
                request.settings = obj.settings;
            end
        end

        function jobCreated(obj,job,userCallback)
            obj.job = job;
            obj.setState("running");
            if ~isempty(userCallback)
                userCallback(job);
            end
        end

        function progressReceived(obj,event,userCallback)
            if isstruct(event) && isscalar(event)
                obj.progress = event;
                obj.progress.status = "running";
            end
            if ~isempty(userCallback)
                userCallback(event);
            end
        end

        function syncJob(obj)
            if isempty(obj.job)
                return
            end
            obj.job = rampSpeed.StudyExecutor.poll(obj.job);
            if any(string(obj.job.state) == ["completed","cancelled","failed"])
                if isstruct(obj.job.study) && isscalar(obj.job.study) && ...
                        ~isempty(fieldnames(obj.job.study))
                    obj.study = obj.job.study;
                end
                obj.setState(obj.job.state);
                if isstruct(obj.progress)
                    obj.progress.status = obj.job.state;
                end
            end
        end

        function setState(obj,value)
            value = string(value);
            if ~isscalar(value)
                value = value(1);
            end
            if string(obj.state) ~= value
                obj.state = value;
                obj.stateHistory(end+1,1) = value;
            end
        end
    end
end

function cases = normalizeInjectedCases(cases)
if isempty(cases)
    cases = repmat(emptySessionCase(),0,1);
    return
end
if ~isstruct(cases)
    error("rampSpeed:invalidCases","cases must be a struct array.");
end
normalized = repmat(emptySessionCase(),numel(cases),1);
for i = 1:numel(cases)
    normalized(i).id = fieldOr(cases(i),"id","case-" + i);
    normalized(i).label = fieldOr(cases(i),"label",normalized(i).id);
    normalized(i).source = fieldOr(cases(i),"source","session");
    normalized(i).designRow = fieldOr(cases(i),"designRow",i);
    normalized(i).sourceIndex = fieldOr(cases(i),"sourceIndex",i);
    normalized(i).carRole = fieldOr(cases(i),"carRole","auto");
    normalized(i).carColumn = fieldOr(cases(i),"carColumn",1);
    normalized(i).derived = fieldOr(cases(i),"derived",struct());
    normalized(i).isBaseline = fieldOr(cases(i),"isBaseline",i == 1);
    if isfield(cases(i),"setupSpec") && ...
            isstruct(cases(i).setupSpec) && ...
            isscalar(cases(i).setupSpec)
        normalized(i).setupSpec = cases(i).setupSpec;
    else
        normalized(i).setupSpec = struct( ...
            "id",string(normalized(i).id), ...
            "label",string(normalized(i).label), ...
            "source",string(normalized(i).source), ...
            "baselineVersion","injected", ...
            "isBaseline",logical(normalized(i).isBaseline));
    end
    normalized(i).setupSpec.id = string(normalized(i).id);
    normalized(i).setupSpec.label = string(normalized(i).label);
    normalized(i).setupSpec.isBaseline = logical(normalized(i).isBaseline);
end
cases = normalized;
end

function value = emptySessionCase()
value = struct("id","","label","","source","session", ...
    "designRow",NaN,"sourceIndex",NaN,"carRole","auto","carColumn",1, ...
    "setupSpec",struct(),"derived",struct(),"isBaseline",false);
end

function tableValue = makeSetupTable(specs)
if isempty(specs)
    tableValue = table(strings(0,1),strings(0,1),strings(0,1), ...
        false(0,1),'VariableNames', ...
        {'id','label','source','isBaseline'});
    return
end
count = numel(specs);
ids = strings(count,1);
labels = strings(count,1);
sources = strings(count,1);
baseline = false(count,1);
for i = 1:count
    ids(i) = fieldOr(specs(i),"id","");
    labels(i) = fieldOr(specs(i),"label",ids(i));
    sources(i) = fieldOr(specs(i),"source","session");
    baseline(i) = logical(fieldOr(specs(i),"isBaseline",i == 1));
end
tableValue = table(ids,labels,sources,baseline, ...
    'VariableNames',{'id','label','source','isBaseline'});
end

function options = normalizeOptions(options)
if nargin < 1 || isempty(options)
    options = struct();
end
if ~isstruct(options) || ~isscalar(options)
    error("rampSpeed:invalidSessionOptions", ...
        "options must be a scalar struct.");
end
end

function value = fieldOr(record,name,default)
value = default;
if isstruct(record) && isfield(record,name) && ~isempty(record.(name))
    value = record.(name);
end
end

function value = baselineVersion(specs)
value = "";
if isempty(specs) || ~isfield(specs,"baselineVersion")
    return
end
values = string({specs.baselineVersion});
values = values(strlength(strtrim(values)) > 0);
if ~isempty(values)
    value = values(1);
end
end

function value = mergeStruct(base,updates)
value = base;
names = fieldnames(updates);
for i = 1:numel(names)
    value.(names{i}) = updates.(names{i});
end
end

function value = normalizeRampType(value)
value = lower(strtrim(string(value)));
if ~isscalar(value)
    error("rampSpeed:unsupportedType", ...
        "rampType must be scalar text.");
end
switch value
    case {"lateral","lateral-limit","lateral limit"}
        value = "lateral";
    case {"longitudinal","pure-longitudinal","pure longitudinal", ...
            "pure_longitudinal"}
        value = "longitudinal";
    otherwise
        error("rampSpeed:unsupportedType", ...
            "Unsupported ramp type: %s",value);
end
end

function value = mapStudyState(status)
status = lower(string(status));
if status == "cancelled"
    value = "cancelled";
elseif status == "failed"
    value = "failed";
else
    value = "completed";
end
end
