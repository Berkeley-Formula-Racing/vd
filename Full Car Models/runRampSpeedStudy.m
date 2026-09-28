function [study,results,figures,events] = runRampSpeedStudy(options)
%RUNRAMPSPEEDSTUDY Compatibility entry point for the canonical Ramp Speed API.
%
%   STUDY = runRampSpeedStudy() runs one protected ramp-speed baseline setup.
%   [STUDY,RESULTS,FIGURES,EVENTS] also returns the normalized run, plot
%   handles, and executor progress events. OPTIONS is a scalar struct; the
%   fields below preserve the useful knobs from the former script:
%
%     rampType, speeds, nRamp, nBisect, mode, solverProfile
%     parallelRequested, numWorkers, checkpointPath
%     loadFromPrev, cacheName, makeFigures, figurePath, figureOptions, outputs
%
% The compatibility entry point owns no solver or legacy result contract. It
% builds the frozen baseline catalog and delegates execution to
% rampSpeed.StudyExecutor/rampSpeed.runStudy.

if nargin < 1 || isempty(options)
    options = struct();
end
options = normalizeOptions(options);

if options.loadFromPrev
    if ~isfile(options.cacheName)
        error("runRampSpeedStudy:missingCache", ...
            "Requested cache does not exist: %s",options.cacheName);
    end
    study = rampSpeed.loadStudy(options.cacheName,options.appVersion);
    results = study.runs;
    events = repmat(emptyEvent(),0,1);
else
    [~,config] = carConfigBaseline();
    [cars,cases] = rampSpeed.buildSetupCatalog(config,config.defaultSetup);
    request = makeRequest(options);

    job = rampSpeed.StudyExecutor.start(cars,cases,request, ...
        struct("onProgress",@printProgress));
    rampSpeed.StudyExecutor.wait(job,Inf);
    if string(job.state) == "failed"
        if ~isempty(job.error)
            rethrow(job.error);
        end
        error("runRampSpeedStudy:failed","The Ramp Speed study failed.");
    end
    study = job.study;
    results = study.runs;
    events = job.progressEvents;
    if strlength(options.cacheName) > 0
        rampSpeed.saveStudy(options.cacheName,study);
    end
end

figures = gobjects(0,1);
if options.makeFigures
    figures = renderStudyFigures(study,options);
    if strlength(options.figurePath) > 0
        saveFigures(options.figurePath,figures,options.figureOptions);
    end
end

if nargout == 0
    clear study results figures events
end
end

function request = makeRequest(options)
settings = struct( ...
    "speeds",double(options.speeds(:).'), ...
    "nRamp",options.nRamp, ...
    "nBisect",options.nBisect, ...
    "mode",options.mode, ...
    "verbose",options.verbose, ...
    "solverProfile",options.solverProfile, ...
    "speedGrid",struct("mode",options.speedGridMode));
request = struct( ...
    "rampType",options.rampType, ...
    "settings",settings, ...
    "parallelRequested",options.parallelRequested, ...
    "numWorkers",options.numWorkers, ...
    "checkpointPath",options.checkpointPath, ...
    "appVersion",options.appVersion, ...
    "runCaseFcn",@rampSpeed.runCase);
end

function figures = renderStudyFigures(study,options)
runs = study.runs;
labels = string({study.cases.label});
metricIds = string(options.outputs(:));
figures = gobjects(0,1);
for i = 1:numel(metricIds)
    metricId = canonicalMetricId(metricIds(i),options.rampType);
    data = rampSpeed.buildPlotData(runs,metricId,struct("labels",labels));
    figureHandle = figure("Name",char(metricId),"NumberTitle","off", ...
        "Visible",options.figureVisibility);
    ax = axes(figureHandle);
    rampSpeed.renderMetric(ax,data,struct( ...
        "showLegend",true,"showWarnings",true));
    figures(end+1,1) = figureHandle; %#ok<AGROW>
end
end

function id = canonicalMetricId(value,rampType)
id = lower(strtrim(string(value)));
if id == "capability"
    if string(rampType) == "longitudinal"
        id = "capability_longitudinal";
    else
        id = "capability_free";
    end
end
end

function printProgress(event)
if ~isstruct(event) || ~isscalar(event) || ...
        ~isfield(event,"message")
    return
end
message = string(event.message);
if isfield(event,"speed_mps") && isfinite(double(event.speed_mps))
    fprintf("[rampSpeed] vCar %.2f m/s: %s\n", ...
        double(event.speed_mps),message);
elseif strlength(strtrim(message)) > 0
    fprintf("[rampSpeed] %s\n",message);
end
end

function options = normalizeOptions(options)
if ~isstruct(options) || ~isscalar(options)
    error("runRampSpeedStudy:invalidOptions", ...
        "options must be a scalar struct.");
end
defaults = struct( ...
    "rampType","longitudinal", ...
    "speeds",[5 10 15 17.5 20 22.5 25], ...
    "nRamp",4, ...
    "nBisect",0, ...
    "mode","coast", ...
    "verbose",false, ...
    "solverProfile","accurate", ...
    "speedGridMode","fixed", ...
    "parallelRequested",false, ...
    "numWorkers",0, ...
    "checkpointPath","", ...
    "appVersion","ramp-speed-app-1", ...
    "loadFromPrev",false, ...
    "cacheName","ramp_speed_study.mat", ...
    "makeFigures",true, ...
    "figureVisibility","on", ...
    "figurePath","figures/ramp_speed_study", ...
    "figureOptions",struct("formats",{{"png","fig"}}, ...
        "resolution",200,"stamp",false), ...
    "outputs",["capability","downforce","drag","aero_balance", ...
        "mechanical_balance","handling_balance","front_camber", ...
        "rear_camber","front_ride_height","rear_ride_height"]);
names = fieldnames(options);
for i = 1:numel(names)
    defaults.(names{i}) = options.(names{i});
end
options = defaults;
options.rampType = lower(strtrim(string(options.rampType)));
if ~ismember(options.rampType,["lateral","longitudinal"])
    error("runRampSpeedStudy:invalidRampType", ...
        "rampType must be lateral or longitudinal.");
end
options.speeds = double(options.speeds(:).');
if isempty(options.speeds) || any(~isfinite(options.speeds)) || ...
        any(options.speeds <= 0) || any(diff(options.speeds) <= 0)
    error("runRampSpeedStudy:invalidSpeeds", ...
        "speeds must be a strictly increasing positive vector.");
end
options.nRamp = requireNonnegativeInteger(options.nRamp,"nRamp",2);
options.nBisect = requireNonnegativeInteger(options.nBisect,"nBisect",0);
options.numWorkers = requireNonnegativeInteger(options.numWorkers, ...
    "numWorkers",0);
options.parallelRequested = logicalScalar(options.parallelRequested, ...
    "parallelRequested");
options.loadFromPrev = logicalScalar(options.loadFromPrev,"loadFromPrev");
options.makeFigures = logicalScalar(options.makeFigures,"makeFigures");
options.checkpointPath = string(options.checkpointPath);
options.cacheName = string(options.cacheName);
options.figurePath = string(options.figurePath);
options.appVersion = string(options.appVersion);
options.figureVisibility = string(options.figureVisibility);
options.outputs = string(options.outputs(:));
end

function value = requireNonnegativeInteger(value,name,minimum)
value = double(value);
if ~isscalar(value) || ~isfinite(value) || value < minimum || ...
        value ~= round(value)
    error("runRampSpeedStudy:invalidOption", ...
        "%s must be an integer >= %d.",name,minimum);
end
end

function value = logicalScalar(value,name)
if ~isscalar(value)
    error("runRampSpeedStudy:invalidOption","%s must be scalar.",name);
end
value = logical(value);
end

function event = emptyEvent()
event = struct("phase","","caseId","","speedIndex",NaN, ...
    "speed_mps",NaN,"completedCases",0,"totalCases",0,"message","");
end
