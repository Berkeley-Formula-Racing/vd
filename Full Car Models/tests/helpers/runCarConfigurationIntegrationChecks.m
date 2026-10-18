function report = runCarConfigurationIntegrationChecks(options)
%RUNCARCONFIGURATIONINTEGRATIONCHECKS Run explicitly selected bounded checks.
% Calling with no enabled phase is an error; test discovery cannot start solvers.
% G-G is never run here. The Explicit phase reuses accepted 28 m/s baseline
% data only after exact lap configuration and production-option equality.

if nargin < 1 || isempty(options), options = struct(); end
if ~isstruct(options) || ~isscalar(options)
    error('runCarConfigurationIntegrationChecks:badOptions', ...
        'options must be a scalar struct.');
end
defaults = struct('runExplicitEventCaseStudy',false,'runRampPair',false, ...
    'evidenceDirectory',"",'explicitDriverDeltaKg',5, ...
    'rampEditedRearSpringRateLbIn',300,'rampSpeedMps',20);
unknown = setdiff(fieldnames(options),fieldnames(defaults));
if ~isempty(unknown)
    error('runCarConfigurationIntegrationChecks:unknownOption', ...
        'Unknown option: %s',strjoin(unknown,', '));
end
names = fieldnames(options);
for i = 1:numel(names), defaults.(names{i}) = options.(names{i}); end
options = defaults;
options.runExplicitEventCaseStudy = logicalScalar(options.runExplicitEventCaseStudy);
options.runRampPair = logicalScalar(options.runRampPair);
if ~options.runExplicitEventCaseStudy && ~options.runRampPair
    error('runCarConfigurationIntegrationChecks:noPhase', ...
        'Enable runExplicitEventCaseStudy or runRampPair explicitly.');
end
validateattributes(options.explicitDriverDeltaKg,{'numeric'}, ...
    {'scalar','real','finite','positive'});
validateattributes(options.rampEditedRearSpringRateLbIn,{'numeric'}, ...
    {'scalar','real','finite','positive'});
validateattributes(options.rampSpeedMps,{'numeric'}, ...
    {'scalar','real','finite','positive'});

helperDir = fileparts(mfilename('fullpath'));
fullCarModelsDir = fileparts(fileparts(helperDir));
repoRoot = fileparts(fullCarModelsDir);
oldPath = path;
oldFolder = pwd;
restoreSession = onCleanup(@() restoreSessionState(oldPath,oldFolder)); %#ok<NASGU>

% Bootstrap all native classes before any MAT load or save.
run(fullfile(fullCarModelsDir,'entrypoints','bootstrap.m'));
[excludedArtifactPaths,sourceOrigins] = filterTask8ArtifactPaths();

sourceOrigins = assertTask8ModelOrigins(fullCarModelsDir);
addpath(helperDir);
if strlength(string(options.evidenceDirectory)) == 0
    evidenceParent = fullfile(repoRoot,'artifacts','car_configuration_refactor');
else
    evidenceParent = char(string(options.evidenceDirectory));
end
if ~isfolder(evidenceParent), mkdir(evidenceParent); end
runId = char(datetime('now','Format',"yyyyMMdd'T'HHmmss"));
outputDirectory = uniqueDirectory(fullfile(evidenceParent, ...
    ['task8-integration-' runId]));
[ok,message] = mkdir(outputDirectory);
if ~ok
    error('runCarConfigurationIntegrationChecks:outputDirectory', ...
        'Could not create evidence folder: %s',message);
end
ignoreId = fopen(fullfile(outputDirectory,'.gitignore'),'w');
if ignoreId < 0
    error('runCarConfigurationIntegrationChecks:ignoreFile', ...
        'Could not create per-evidence-folder .gitignore.');
end
fprintf(ignoreId,'*\n');
fclose(ignoreId);

report = struct('runId',string(runId),'created',datetime('now'), ...
    'repositoryRoot',string(repoRoot),'matlabVersion',string(version), ...
    'options',options,'outputDirectory',string(outputDirectory), ...
    'phases',strings(0,1),'phaseStatus',strings(0,1), ...
    'phaseElapsedSeconds',zeros(0,1),'phaseDetails',strings(0,1));
report.excludedArtifactPathEntries = excludedArtifactPaths;
report.sourceOrigins = sourceOrigins;
if options.runExplicitEventCaseStudy
    t = tic;
    try
        report.explicitCaseStudy = runExplicitEventCaseStudy(repoRoot,options);
        phaseStatus = "pass";
        phaseDetails = report.explicitCaseStudy.acceptanceDetails;
        if ~report.explicitCaseStudy.accepted
            phaseStatus = "partial";
        end
    catch exception
        report.explicitCaseStudyError = exceptionRecord(exception);
        phaseStatus = "error-captured";
        phaseDetails = "Explicit case study raised " + string(exception.identifier) + ...
            ": " + string(exception.message);
    end
    report = appendPhase(report,"explicit-two-row-events",phaseStatus,toc(t),phaseDetails);
    report = persistIntegrationReport(report,outputDirectory);
end
if options.runRampPair
    t = tic;
    try
        report.rampPair = runRampPair(options);
        phaseStatus = "pass";
        phaseDetails = report.rampPair.acceptanceDetails;
        if ~report.rampPair.accepted, phaseStatus = "partial"; end
    catch exception
        report.rampPairError = exceptionRecord(exception);
        phaseStatus = "error-captured";
        phaseDetails = "Ramp pair raised " + string(exception.identifier) + ...
            ": " + string(exception.message);
    end
    report = appendPhase(report,"ramp-baseline-and-edited-setup",phaseStatus,toc(t),phaseDetails);
    report = persistIntegrationReport(report,outputDirectory);
end
report = persistIntegrationReport(report,outputDirectory);
fprintf('TASK8_INTEGRATION_REPORT=%s\n',report.reportPath);
fprintf('TASK8_INTEGRATION_SUMMARY=%s\n',report.summaryPath);
for i = 1:height(report.phaseTable)
    fprintf('TASK8_PHASE %s status=%s elapsed_s=%.3f details=%s\n', ...
        report.phases(i),report.phaseStatus(i), ...
        report.phaseElapsedSeconds(i),report.phaseDetails(i));
end
end

function evidence = runExplicitEventCaseStudy(repoRoot,options)
afterPath = fullfile(repoRoot,'artifacts','car_configuration_refactor', ...
    '2026-10-05-after-production-28mps-retry1','numerical-reference.mat');
if ~isfile(afterPath)
    error('runCarConfigurationIntegrationChecks:missingAfterReference', ...
        'Accepted after numerical reference is missing: %s',afterPath);
end
n = load(afterPath,'numerical');
n = n.numerical;
if ~isfield(n,'qss') || ~strcmp(n.qss.ggStatus,'pass') || ...
        ~strcmp(n.qss.eventStatus,'pass') || ~strcmp(n.qss.telemetryStatus,'pass')
    error('runCarConfigurationIntegrationChecks:invalidAfterReference', ...
        'Accepted after numerical reference is incomplete.');
end
dependencyGate = verifyAcceptedPhysicsDependencies(repoRoot);
constructionPath = char(string(n.metadata.constructionReferenceFile));
if ~isfile(constructionPath)
    error('runCarConfigurationIntegrationChecks:missingConstructionReference', ...
        'Construction reference named by the numerical result is missing.');
end
r = load(constructionPath,'reference');
r = r.reference;
ggOptions = ggProductionOptions();
if ~isequaln(ggOptions,n.input.ggOptions) || n.input.ggWorkerCount ~= 0
    error('runCarConfigurationIntegrationChecks:ggOptionsMismatch', ...
        'Cached G-G options differ from the accepted serial production run.');
end

[~,eventParams,baseDesign] = carConfig();
if height(baseDesign) ~= 1
    error('runCarConfigurationIntegrationChecks:baseDesignShape', ...
        'Expected one scalar default design row; got %d.',height(baseDesign));
end
rows = [baseDesign(1,:);baseDesign(1,:)];
rows.accel_driver_weight(2) = rows.accel_driver_weight(1) + ...
    options.explicitDriverDeltaKg;
changed = strings(0,1);
for k = 1:width(rows)
    name = string(rows.Properties.VariableNames{k});
    if ~isequaln(rows{1,k},rows{2,k}), changed(end+1,1) = name; end %#ok<AGROW>
end
if ~isequal(changed,"accel_driver_weight")
    error('runCarConfigurationIntegrationChecks:unexpectedDesignChange', ...
        'Only accel_driver_weight may change; changed: %s',strjoin(changed,', '));
end
[cases,eventParams2,design] = carConfig( ...
    "Explicit",rows,"legacy","nominal","cases");
if numel(cases) ~= 2 || height(design) ~= 2
    error('runCarConfigurationIntegrationChecks:explicitShape', ...
        'Expected exactly two Explicit cases.');
end
caseIds = string({cases.id});
if numel(unique(caseIds)) ~= 2
    error('runCarConfigurationIntegrationChecks:duplicateCaseIds', ...
        'Explicit cases must have distinct case IDs.');
end
if ~isequaln(eventParams2,eventParams) || ...
        ~isequaln(eventParams2,r.qss.default.eventParams)
    error('runCarConfigurationIntegrationChecks:eventParamsMismatch', ...
        'Explicit event parameters differ from the accepted baseline.');
end
snapshot1 = carConfigurationSnapshot(cases(1).lapCar);
snapshot2 = carConfigurationSnapshot(cases(2).lapCar);
if ~isequaln(snapshot1,snapshot2) || ...
        ~isequaln(snapshot1,r.qss.default.lapCar)
    error('runCarConfigurationIntegrationChecks:lapCacheGate', ...
        'Lap configuration is not exactly equal to the validated baseline.');
end
massDelta = cases(2).accelCar.M-cases(1).accelCar.M;
if abs(massDelta-options.explicitDriverDeltaKg) > 1e-12
    error('runCarConfigurationIntegrationChecks:accelMassDelta', ...
        'Acceleration car mass did not change by the requested driver delta.');
end

solved = n.qss.solvedLapCar;
assertSolvedCarMatchesConfiguration(cases(1).lapCar,solved);
torqueCurveExact = isequaln(cases(1).lapCar.powertrain.torque_fn, ...
    n.qss.solvedLapCar.powertrain.torque_fn);
if ~torqueCurveExact
    error('runCarConfigurationIntegrationChecks:torqueCurveMismatch', ...
        'Current fixed torque data differs from the accepted and construction references.');
end
ggFields = ["ggPoints","ss_info","accel_info","decel_info", ...
    "longAccelLookup","longDecelLookup","ggMask"];
for f = ggFields
    if ~isprop(solved,char(f))
        error('runCarConfigurationIntegrationChecks:missingGGField', ...
            'Accepted solved car does not expose %s.',f);
    end
end
study = DOEStudyConfig();
eventNames = string(study.events(:).');
if ~all(ismember(["skidpad","accel","autocross","endurance"],eventNames))
    error('runCarConfigurationIntegrationChecks:missingScoreEvent', ...
        'All four events are required for the configured score.');
end
objects = cell(2,1);
caseEvidence = cell(2,1);
selectedReferenceEvents = ["skidpad","accel","autocross"];
for i = 1:2
    car = cases(i).lapCar;
    for f = ggFields, car.(char(f)) = solved.(char(f)); end
    cacheCopyExact = isequaln(car.ggPoints,solved.ggPoints) && ...
        isequaln(car.ggMask,solved.ggMask) && ...
        isequaln(car.longAccelLookup,solved.longAccelLookup) && ...
        isequaln(car.longDecelLookup,solved.longDecelLookup);
    if ~cacheCopyExact
        error('runCarConfigurationIntegrationChecks:ggCacheCopy', ...
            'Validated G-G payload did not copy exactly to case %d.',i);
    end
    eventClock = tic;
    ev = Events2(car,cases(i).accelCar,eventParams2);
    for eventName = eventNames
        switch eventName
            case "skidpad", ev.Skidpad();
            case "accel", ev.Accel();
            case "autocross", ev.Autocross();
            case "endurance", ev.Endurance();
        end
    end
    points = ev.computePoints();
    car.comp = ev;
    [metrics,~] = doeCaseMetrics(car,i,[]);
    [score,breakdown] = doeScoreCase(points,metrics,study.objective);
    objects{i} = ev;
    timeNames = ["skidpad","accel","autocross","endurance"];
    eventTimeValues = zeros(1,numel(timeNames));
    for j = 1:numel(timeNames)
        eventTimeValues(j) = double(ev.times.(char(timeNames(j))));
    end
    eventTimesFinitePositive = all(isfinite(eventTimeValues) & eventTimeValues > 0);
    pointValues = cellfun(@double,struct2cell(points));
    pointsFinite = all(isfinite(pointValues));
    breakdownFinite = finiteNumericScalarStruct(breakdown);
    scoreFinite = isscalar(score) && isfinite(score);
    metricValues = [metrics.t_autox(1),metrics.t_accel(1), ...
        metrics.t_skid(1),metrics.total_work_kJ(1),metrics.gg_coverage(1)];
    metricsFinite = all(isfinite(metricValues));
    metricValid = logical(metrics.valid(1));
    referenceEventDeltas = struct();
    referenceTimesMatch = true;
    if i == 1
        for eventName = selectedReferenceEvents
            name = char(eventName);
            delta = eventTimeValues(timeNames == eventName) - ...
                double(n.qss.eventTimes.(name));
            referenceEventDeltas.(name) = delta;
            referenceTimesMatch = referenceTimesMatch && abs(delta) <= 1e-8;
        end
    end
    caseAccepted = eventTimesFinitePositive && pointsFinite && ...
        breakdownFinite && scoreFinite && metricsFinite && metricValid && ...
        (i ~= 1 || referenceTimesMatch);
    failedChecks = strings(0,1);
    if ~eventTimesFinitePositive, failedChecks(end+1) = "event times nonfinite or nonpositive"; end %#ok<AGROW>
    if ~pointsFinite, failedChecks(end+1) = "event points nonfinite"; end %#ok<AGROW>
    if ~breakdownFinite, failedChecks(end+1) = "score breakdown nonfinite"; end %#ok<AGROW>
    if ~scoreFinite, failedChecks(end+1) = "objective score nonfinite"; end %#ok<AGROW>
    if ~metricsFinite, failedChecks(end+1) = "required DOE metrics nonfinite"; end %#ok<AGROW>
    if ~metricValid, failedChecks(end+1) = "DOE metric row is invalid"; end %#ok<AGROW>
    if i == 1 && ~referenceTimesMatch, failedChecks(end+1) = "baseline selected event times differ from accepted reference"; end %#ok<AGROW>
    caseEvidence{i} = struct('caseIndex',i,'caseId',string(cases(i).id), ...
        'designRow',design(i,:),'lapMassKg',cases(i).lapCar.M, ...
        'accelMassKg',cases(i).accelCar.M,'eventNames',eventNames, ...
        'eventTimes',ev.times,'eventTimeFields',timeNames, ...
        'eventTimesFinitePositive',eventTimesFinitePositive, ...
        'baselineReferenceEventDeltas',referenceEventDeltas, ...
        'baselineReferenceTimesMatch',referenceTimesMatch, ...
        'eventArraySizes',eventArraySizes(ev),'points',points, ...
        'pointsFinite',pointsFinite,'score',score,'scoreFinite',scoreFinite, ...
        'scoreBreakdown',breakdown,'scoreBreakdownFinite',breakdownFinite, ...
        'metricRow',metrics,'requiredMetricsFinite',metricsFinite, ...
        'metricRowValid',metricValid,'eventElapsedSeconds',toc(eventClock), ...
        'freshEvents2',isa(ev,'Events2'),'accepted',caseAccepted, ...
        'failedChecks',failedChecks);
end
freshDistinctEvents2 = isa(objects{1},'Events2') && isa(objects{2},'Events2') && ...
    ~(objects{1} == objects{2});
caseResults = [caseEvidence{:}];
allCasesAccepted = all([caseResults.accepted]);
accepted = freshDistinctEvents2 && allCasesAccepted;
acceptanceIssues = strings(0,1);
if ~freshDistinctEvents2, acceptanceIssues(end+1) = "Events2 handles are not distinct"; end %#ok<AGROW>
if ~allCasesAccepted
    acceptanceIssues(end+1) = "one or more Explicit rows failed finite/status/reference gates"; %#ok<AGROW>
end
if isempty(acceptanceIssues), acceptanceDetails = "Both rows passed finite event, point, score, metric, case-index, fresh-handle, and baseline selected-time gates.";
else, acceptanceDetails = strjoin(acceptanceIssues,"; "); end

evidence = struct('referencePath',string(afterPath), ...
    'constructionReferencePath',string(constructionPath), ...
    'referenceGitCommit',string(n.metadata.gitCommit), ...
    'dependencyHashGate',dependencyGate, ...
    'fixedTorqueCurveExactlyMatchesReferences',torqueCurveExact, ...
    'ggOptions',ggOptions,'ggWorkerCount',n.input.ggWorkerCount, ...
    'ggReuse',"matched validated baseline cache; no G-G per Explicit row", ...
    'ggFieldsReused',ggFields,'ggRawSize',size(n.qss.ggParamSets), ...
    'ggPointSize',size(solved.ggPoints),'ggMask',solved.ggMask, ...
    'designRows',design,'changedDesignFields',changed, ...
    'expectedAccelMassDeltaKg',options.explicitDriverDeltaKg, ...
    'observedAccelMassDeltaKg',massDelta,'caseIdsDistinct',true, ...
    'eventParamsEqualBaseline',true, ...
    'lapConfigurationEqualValidatedBaseline',true, ...
    'freshDistinctEvents2Handles',freshDistinctEvents2, ...
    'accepted',accepted,'acceptanceDetails',acceptanceDetails, ...
    'cases',caseResults,'eventsPath',string(which('Events2')), ...
    'metricsPath',string(which('doeCaseMetrics')), ...
    'scorePath',string(which('doeScoreCase')));
end

function gate = verifyAcceptedPhysicsDependencies(repoRoot)
acceptedPath = fullfile(repoRoot,'artifacts','car_configuration_refactor', ...
    '2026-10-05-after-production-28mps-retry1','dependency-hashes-postrun.csv');
reviewedPath = fullfile(repoRoot,'artifacts','car_configuration_refactor', ...
    '2026-10-05-doe-expensive','post-run-sha256-comparison.tsv');
accepted = readtable(acceptedPath,'TextType','string');
reviewed = readtable(reviewedPath,'FileType','text','Delimiter','\t', ...
    'TextType','string','VariableNamingRule','preserve');
if height(accepted) ~= 177 || height(reviewed) ~= 273
    error('runCarConfigurationIntegrationChecks:dependencyManifestSize', ...
        'Expected accepted/reviewed manifests of 177/273 rows; found %d/%d.', ...
        height(accepted),height(reviewed));
end
allowedDrift = ["Full Car Models/+rampSpeed/releaseManifest.m"; ...
    "Full Car Models/+vehicleConfig/buildCar.m"; ...
    "Full Car Models/+vehicleConfig/fromRampSetup.m"; ...
    "Full Car Models/entrypoints/SteadyStateLapsim.m"; ...
    "Full Car Models/setup_paths.m"];
reviewedDrift = string(reviewed.RelativePath( ...
    lower(string(reviewed.PreSHA256)) ~= lower(string(reviewed.PostSHA256))));
if ~isequal(sort(reviewedDrift),sort(allowedDrift))
    error('runCarConfigurationIntegrationChecks:reviewedDriftSet', ...
        'Reviewed DOE manifest drift set differs from the five explicitly reviewed source paths.');
end
for driftPath = allowedDrift.'
    driftIndex = find(string(reviewed.RelativePath) == driftPath);
    if numel(driftIndex) ~= 1
        error('runCarConfigurationIntegrationChecks:missingDriftReview', ...
            'A known drift path is missing or duplicated in DOE manifest: %s',driftPath);
    end
    driftFile = fullfile(repoRoot,replace(driftPath,'/',filesep));
    if ~isfile(driftFile) || ...
            sha256File(driftFile) ~= lower(string(reviewed.PostSHA256(driftIndex)))
        error('runCarConfigurationIntegrationChecks:driftPostHashMismatch', ...
            'Current known-drift file differs from its reviewed post-run hash: %s',driftPath);
    end
end
observedDrift = strings(0,1);
for i = 1:height(accepted)
    relativePath = string(accepted.repoRelativePath(i));
    reviewIndex = find(string(reviewed.RelativePath) == relativePath);
    if numel(reviewIndex) > 1
        error('runCarConfigurationIntegrationChecks:dependencyReviewCoverage', ...
            'Reviewed dependency manifest has %d rows for %s.', ...
            numel(reviewIndex),relativePath);
    end
    filePath = fullfile(repoRoot,replace(relativePath,'/',filesep));
    if ~isfile(filePath)
        error('runCarConfigurationIntegrationChecks:dependencyMissing', ...
            'Accepted physics dependency is missing: %s',filePath);
    end
    currentHash = sha256File(filePath);
    acceptedHash = lower(string(accepted.sha256(i)));
    if isempty(reviewIndex)
        if currentHash ~= acceptedHash
            error('runCarConfigurationIntegrationChecks:unreviewedDependencyDrift', ...
                'Dependency absent from the DOE comparison differs from accepted hash: %s', ...
                relativePath);
        end
        continue
    end
    reviewedHash = lower(string(reviewed.PostSHA256(reviewIndex)));
    if currentHash ~= reviewedHash
        error('runCarConfigurationIntegrationChecks:dependencyChangedAfterReview', ...
            'Current dependency differs from the reviewed post-run hash: %s',relativePath);
    end
    if currentHash ~= acceptedHash
        if ~ismember(relativePath,allowedDrift)
            error('runCarConfigurationIntegrationChecks:unapprovedDependencyDrift', ...
                'Unreviewed accepted-baseline dependency drift: %s',relativePath);
        end
        observedDrift(end+1,1) = relativePath; %#ok<AGROW>
    end
end
if any(~ismember(observedDrift,allowedDrift))
    error('runCarConfigurationIntegrationChecks:dependencyDriftSet', ...
        'Accepted baseline changed outside reviewed paths: %s.', ...
        strjoin(observedDrift,', '));
end
torquePath = "Full Car Models/carComponents/KTM450.m";
torqueIndex = find(string(accepted.repoRelativePath) == torquePath);
if numel(torqueIndex) ~= 1
    error('runCarConfigurationIntegrationChecks:torqueSourceNotReviewed', ...
        'KTM450 torque source is missing or duplicated in accepted manifest.');
end
gate = struct('acceptedManifestPath',string(acceptedPath), ...
    'reviewedManifestPath',string(reviewedPath), ...
    'acceptedDependencyCount',height(accepted), ...
    'reviewedDependencyCount',height(reviewed), ...
    'currentHashesMatchReviewedPostRun',true, ...
    'reviewedCommentPathDriftPaths',reviewedDrift, ...
    'acceptedBaselineDriftPaths',observedDrift, ...
    'ktm450TorqueSourcePath',torquePath, ...
    'ktm450TorqueSourceSha256',lower(string(accepted.sha256(torqueIndex))));
end

function value = sha256File(filePath)
fid = fopen(char(filePath),'r');
if fid < 0
    error('runCarConfigurationIntegrationChecks:dependencyReadFailed', ...
        'Unable to read dependency file: %s',filePath);
end
cleanup = onCleanup(@()fclose(fid)); %#ok<NASGU>
bytes = fread(fid,Inf,'*uint8');
digest = javaMethod('getInstance','java.security.MessageDigest','SHA-256');
digest.update(bytes);
digestBytes = typecast(digest.digest(),'uint8');
value = lower(string(reshape(dec2hex(digestBytes,2).',1,[])));
end

function assertSolvedCarMatchesConfiguration(candidate,solved)
outputFields = ["ggPoints","ss_info","accel_info","decel_info", ...
    "longAccelLookup","longDecelLookup","ggMask","comp"];
reference = solved;
for f = outputFields
    name = char(f);
    if isprop(candidate,name), reference.(name) = candidate.(name); end
end
if ~isequaln(carConfigurationSnapshot(candidate), ...
        carConfigurationSnapshot(reference))
    error('runCarConfigurationIntegrationChecks:solvedCarMismatch', ...
        'Accepted solved car differs from the current exact lap configuration.');
end
end

function valid = finiteNumericScalarStruct(value)
valid = isstruct(value) && isscalar(value);
if ~valid, return, end
fields = fieldnames(value);
for i = 1:numel(fields)
    item = value.(fields{i});
    if ~isnumeric(item) || ~isscalar(item) || ~isfinite(item)
        valid = false;
        return
    end
end
end

function sizes = eventArraySizes(ev)
sizes = struct();
names = ["skidpad","accel","autocross","endurance"];
for eventName = names
    name = char(eventName);
    value = ev.(name);
    sizes.(name) = struct();
    if isstruct(value) && isscalar(value)
        fields = fieldnames(value);
        for i = 1:numel(fields)
            item = value.(fields{i});
            if isnumeric(item) || islogical(item) || istable(item)
                sizes.(name).(fields{i}) = size(item);
            end
        end
    end
end
end

function evidence = runRampPair(options)
[~,config] = carConfigBaseline();
baseline = config.defaultSetup;
edited = baseline;
edited.id = "ramp-rear-spring-300";
edited.label = "Rear spring 300 lb/in";
allowedRearSprings = double(config.options.rearSpringRate_lb_in(:));
if ~ismember(options.rampEditedRearSpringRateLbIn,allowedRearSprings) || ...
        options.rampEditedRearSpringRateLbIn == baseline.rearSpringRate_lb_in
    error('runCarConfigurationIntegrationChecks:rampSpringChoice', ...
        'Edited rear spring must be a different value from the legal discrete catalog.');
end
edited.rearSpringRate_lb_in = options.rampEditedRearSpringRateLbIn;
[cars,cases] = rampSpeed.buildSetupCatalog(config,[baseline;edited]);
if numel(cars) ~= 2 || numel(cases) ~= 2 || ...
        ~isa(cars{1},'Car') || ~isa(cars{2},'Car') || isequaln(cars{1},cars{2})
    error('runCarConfigurationIntegrationChecks:rampConstruction', ...
        'Expected distinct frozen-baseline and edited setup cars.');
end
modes = ["coast","balanced"];
modeEvidence = cell(numel(modes),1);
for k = 1:numel(modes)
    settings = struct('speeds',options.rampSpeedMps,'nRamp',4,'nBisect',0, ...
        'mode',modes(k),'verbose',false,'solverProfile',"fastPreview", ...
        'speedGrid',struct('mode',"fixed"));
    request = struct('rampType',"lateral",'settings',settings, ...
        'parallelRequested',false,'numWorkers',0,'checkpointPath',"", ...
        'appVersion',"task8-integration",'runCaseFcn',@rampSpeed.runCase);
    stageClock = tic;
    [study,events] = rampSpeed.runStudy(cars,cases,request,struct());
    if ~isfield(study,'runs') || numel(study.runs) ~= 2
        error('runCarConfigurationIntegrationChecks:rampRunCount', ...
            'Expected two ramp runs for mode %s.',modes(k));
    end
    perCase = cell(2,1);
    for i = 1:2
        runResult = study.runs(i);
        T = runResult.perSpeed;
        if height(T) ~= 1 || abs(T.speed_mps(1)-options.rampSpeedMps) > 1e-12
            error('runCarConfigurationIntegrationChecks:rampSpeedRow', ...
                'Unexpected speed row for %s case %d.',modes(k),i);
        end
        required = ["status","solver_status","valid","aLat_free_mps2", ...
            "K_linear_rad_per_mps2","cuo_steer_linear_rad", ...
            "cuo_steer_limit_rad","mechanical_balance_front", ...
            "front_camber_rad","rear_camber_rad","front_ride_height_m", ...
            "rear_ride_height_m","ramp_complete_fraction"];
        if ~all(ismember(required,string(T.Properties.VariableNames)))
            error('runCarConfigurationIntegrationChecks:rampFields', ...
                'Ramp result is missing required status, mask, steering, balance, camber, or ride-height fields.');
        end
        numericFields = ["aLat_free_mps2","K_linear_rad_per_mps2", ...
            "cuo_steer_linear_rad","cuo_steer_limit_rad", ...
            "mechanical_balance_front","front_camber_rad","rear_camber_rad", ...
            "front_ride_height_m","rear_ride_height_m", ...
            "ramp_complete_fraction"];
        finiteMasks = struct();
        finiteMetrics = true;
        for field = numericFields
            name = char(field);
            value = T.(name)(1);
            finiteMasks.(name) = isfinite(value);
            finiteMetrics = finiteMetrics && isfinite(value);
        end
        validMask = logical(T.valid(1));
        solverRowStatus = lower(string(T.status(1)));
        solverDetailStatus = lower(string(T.solver_status(1)));
        runStatus = lower(string(runResult.status));
        studyStatus = lower(string(study.status));
        acceptedSolverStatuses = ["converged","near_feasible","complete","completed"];
        solverStatusAccepted = ismember(solverRowStatus,acceptedSolverStatuses) && ...
            (strlength(solverDetailStatus) == 0 || ...
            ismember(solverDetailStatus,acceptedSolverStatuses));
        schedulerStatusAccepted = runStatus == "complete" && studyStatus == "complete";
        rowAccepted = schedulerStatusAccepted && solverStatusAccepted && ...
            validMask && finiteMetrics;
        issues = strings(0,1);
        if ~schedulerStatusAccepted
            issues(end+1) = "study/run status=" + studyStatus + "/" + runStatus; %#ok<AGROW>
        end
        if ~solverStatusAccepted
            issues(end+1) = "solver status=" + solverRowStatus + "/" + solverDetailStatus; %#ok<AGROW>
        end
        if ~validMask, issues(end+1) = "valid mask is false"; end %#ok<AGROW>
        if ~finiteMetrics
            missingFinite = numericFields(~cellfun(@(name)finiteMasks.(name),cellstr(numericFields)));
            issues(end+1) = "nonfinite fields=" + strjoin(missingFinite,","); %#ok<AGROW>
        end
        perCase{i} = struct('id',string(cases(i).id), ...
            'label',string(cases(i).label), ...
            'rearSpringRateLbIn',cases(i).setupSpec.rearSpringRate_lb_in, ...
            'studyStatus',studyStatus,'runStatus',runStatus, ...
            'solverRowStatus',solverRowStatus, ...
            'solverDetailStatus',solverDetailStatus,'validMask',validMask, ...
            'finiteMasks',finiteMasks,'finiteMetrics',finiteMetrics, ...
            'accepted',rowAccepted,'issues',issues,'perSpeed',T, ...
            'runMeta',runResult.runMeta);
    end
    modeCases = [perCase{:}];
    modeAccepted = all([modeCases.accepted]);
    if modeAccepted, modeDetails = "both cases complete with valid finite steering, balance, camber, ride-height, and capability outputs";
    else
        detailParts = strings(2,1);
        for i = 1:2
            detailParts(i) = modeCases(i).id + ": " + strjoin(modeCases(i).issues,"; ");
        end
        modeDetails = strjoin(detailParts," | ");
    end
    modeEvidence{k} = struct('mode',modes(k), ...
        'studyStatus',string(study.status),'progressEventCount',numel(events), ...
        'elapsedSeconds',toc(stageClock),'accepted',modeAccepted, ...
        'acceptanceDetails',modeDetails,'cases',modeCases);
end
modeResults = [modeEvidence{:}];
accepted = all([modeResults.accepted]);
if accepted, acceptanceDetails = "Coast and balanced each completed both cases with valid finite steering, balance, camber, ride-height, and capability outputs.";
else
    detailParts = strings(numel(modeResults),1);
    for i = 1:numel(modeResults)
        detailParts(i) = modeResults(i).mode + ": " + modeResults(i).acceptanceDetails;
    end
    acceptanceDetails = strjoin(detailParts," | ");
end
evidence = struct('rampType',"lateral",'speedMps',options.rampSpeedMps, ...
    'gridMode',"fixed",'nRamp',4,'nBisect',0, ...
    'solverProfile',"fastPreview",'parallel',false, ...
    'baseSetup',baseline,'editedSetup',edited, ...
    'editedField',"rearSpringRate_lb_in", ...
    'editedBaselineValueLbIn',baseline.rearSpringRate_lb_in, ...
    'editedValueLbIn',edited.rearSpringRate_lb_in, ...
    'capabilityMetricMeaning', ...
    "ramp_complete_fraction is numeric gLat_top/max(gLat_max,eps), not a boolean or completed-row count", ...
    'accepted',accepted,'acceptanceDetails',acceptanceDetails, ...
    'modes',modeResults,'carPaths',string(which('rampSpeed.buildCarFromSetup')), ...
    'runCasePath',string(which('rampSpeed.runCase')));
end

function report = appendPhase(report,name,status,elapsed,details)
report.phases(end+1,1) = string(name);
report.phaseStatus(end+1,1) = string(status);
report.phaseElapsedSeconds(end+1,1) = double(elapsed);
report.phaseDetails(end+1,1) = string(details);
end

function report = persistIntegrationReport(report,outputDirectory)
report.phaseTable = table(report.phases,report.phaseStatus, ...
    report.phaseElapsedSeconds,report.phaseDetails, ...
    'VariableNames',{'phase','status','elapsed_s','details'});
report.reportPath = string(fullfile(outputDirectory,'integration-report.mat'));
report.summaryPath = string(fullfile(outputDirectory,'integration-summary.csv'));
save(char(report.reportPath),'report','-v7.3');
writetable(report.phaseTable,char(report.summaryPath));
end

function value = exceptionRecord(exception)
value = struct('identifier',string(exception.identifier), ...
    'message',string(exception.message), ...
    'report',string(getReport(exception,'extended','hyperlinks','off')));
end

function output = uniqueDirectory(requested)
output = char(requested);
if ~isfolder(output) && ~isfile(output), return; end
for i = 1:1000
    candidate = sprintf('%s-%02d',requested,i);
    if ~isfolder(candidate) && ~isfile(candidate)
        output = candidate;
        return
    end
end
error('runCarConfigurationIntegrationChecks:outputCollision', ...
    'Could not choose a unique evidence directory.');
end

function value = logicalScalar(value)
if ~isscalar(value)
    error('runCarConfigurationIntegrationChecks:invalidOption', ...
        'Boolean options must be scalar.');
end
value = logical(value);
end

function [excluded,origins] = filterTask8ArtifactPaths()
entries = string(strsplit(path,pathsep));
normalized = lower(replace(entries,"\","/"));
mask = contains(normalized,"/.superpowers/") | endsWith(normalized,"/.superpowers") | ...
    contains(normalized,"/artifacts/") | endsWith(normalized,"/artifacts");
excluded = entries(mask);
for i = 1:numel(excluded)
    rmpath(char(excluded(i)));
end
remaining = string(strsplit(path,pathsep));
remainingNormalized = lower(replace(remaining,"\","/"));
if any(contains(remainingNormalized,"/.superpowers/") | ...
        endsWith(remainingNormalized,"/.superpowers") | ...
        contains(remainingNormalized,"/artifacts/") | endsWith(remainingNormalized,"/artifacts"))
    error('runCarConfigurationIntegrationChecks:artifactPathRemains', ...
        'Artifact or .superpowers directories remain on the MATLAB path.');
end
origins = strings(0,2);
end

function origins = assertTask8ModelOrigins(fullCarModelsDir)
names = ["Car";"Tire2";"Aero";"Powertrain";"Events2";"gg2";"makeGG"; ...
    "carConfig";"parameters_loop";"doeRunCase";"doeCaseMetrics";"doeScoreCase"; ...
    "vehicleConfig.defaultSetup";"vehicleConfig.sampleSetups"; ...
    "vehicleConfig.buildCases";"vehicleConfig.buildCar"; ...
    "vehicleConfig.fromRampSetup";"vehicleConfig.loadAssets"; ...
    "rampSpeed.buildCarFromSetup";"rampSpeed.executeSpeedPlan"; ...
    "rampSpeed.loadRampAssets";"rampSpeed.runCase"];
expectedPrefix = lower(replace(string(fullCarModelsDir),"\","/")) + "/";
resolved = strings(size(names));
for i = 1:numel(names)
    resolved(i) = string(which(char(names(i))));
    normalized = lower(replace(resolved(i),"\","/"));
    if strlength(resolved(i)) == 0 || ~startsWith(normalized,expectedPrefix) || ...
            contains(normalized,"/.superpowers/") || contains(normalized,"/artifacts/")
        error('runCarConfigurationIntegrationChecks:unexpectedModelOrigin', ...
            'Expected %s under %s; resolved %s.',names(i),expectedPrefix,resolved(i));
    end
end
origins = table(names,resolved,'VariableNames',{'name','resolvedPath'});
fprintf('TASK8_MODEL_ORIGIN_ASSERTIONS=%d\n',height(origins));
end
function restoreSessionState(oldPath,oldFolder)
path(oldPath);
cd(oldFolder);
end
