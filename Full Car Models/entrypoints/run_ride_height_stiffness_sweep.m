function result = run_ride_height_stiffness_sweep(options)
%RUN_RIDE_HEIGHT_STIFFNESS_SWEEP Sweep ride height and spring rates.
%
%   RESULT = RUN_RIDE_HEIGHT_STIFFNESS_SWEEP() runs the default Cartesian
%   sweep for the ride-height aero map and saves a compact MAT result, CSV
%   case table, and summary figure under Full Car Models/sweeps.
%
%   RESULT = RUN_RIDE_HEIGHT_STIFFNESS_SWEEP(OPTIONS) accepts:
%       preview                 return the design without solving (false)
%       bootstrap               configure model paths and repository cwd (true)
%       aeroMapPath             active aero-map CSV path
%       numRideLevels           levels per ride-height axis (5)
%       frontRideHeights_in     explicit front ride-height levels
%       rearRideHeights_in      explicit rear ride-height levels
%       frontSpringRates_lb_in  front damper spring rates (250,300,350)
%       rearSpringRates_lb_in   rear damper spring rates (200,250,300)
%       events                  skidpad, accel, autocross by default
%       numWorkers              DOE case workers (8)
%       allowSerialFallback     use serial execution without PCT (true)
%       requireMapCoverage      reject points outside scattered coverage (false)
%       outputDirectory         result directory
%       outputPath              requested MAT result path
%       overwrite               replace existing output files (false)
%       saveRawResults          retain full per-case car/event objects (false)
%       plot                    save the summary figure (true)
%       verbose                 print progress (true)
%
% Ride heights are static offsets in inches relative to the aero-map
% reference. Spring rates are damper rates in lb/in. The default ride axes
% are generated from the rectangular bounds exposed by AeroMap, so this
% entry point never creates a requested static ride height outside that
% envelope.

if nargin < 1 || isempty(options)
    options = struct();
end

entrypointDir = fileparts(mfilename('fullpath'));
fullCarModelsDir = fileparts(entrypointDir);
repoRoot = fileparts(fullCarModelsDir);
options = normalizeOptions(options,fullCarModelsDir,repoRoot);

if options.bootstrap
    run(fullfile(entrypointDir,'bootstrap.m'));
end

map = AeroMap(options.aeroMapPath);
mapInfo = makeMapInfo(map);
[designTable,settings] = makeDesign(map,options);
mapPreflight = preflightMap(map,designTable,options);

if options.preview
    result = struct( ...
        'designTable',designTable, ...
        'caseCount',height(designTable), ...
        'mapBounds',struct('front',map.frontRangeIn,'rear',map.rearRangeIn), ...
        'mapInfo',mapInfo, ...
        'mapPreflight',mapPreflight, ...
        'settings',settings, ...
        'output',emptyOutput());
    return
end

effectiveWorkers = resolveWorkers(options.numWorkers,options.allowSerialFallback);
settings.requestedWorkers = options.numWorkers;
settings.effectiveWorkers = effectiveWorkers;

if options.verbose
    fprintf('Ride-height/spring-rate sweep: %d cases, %d worker(s)\n', ...
        height(designTable),effectiveWorkers);
    if mapPreflight.unsupportedCoverageCount > 0
        fprintf(['  %d static points are inside the rectangular map bounds but ' ...
            'outside scattered coverage.\n'],mapPreflight.unsupportedCoverageCount);
    end
end

study = DOEStudyConfig();
study.events = options.events;
study.numWorkers = effectiveWorkers;
study.ramps.enabled = false;
study.ramps.saveFullPoints = false;

job = simLog.start('ride height stiffness sweep');
totalTimer = tic;
try
    [carCell,eventParams,returnedDesign] = carConfig('Explicit',designTable);
    sweepNames = designTable.Properties.VariableNames;
    returnedHasSweep = all(ismember(sweepNames, ...
        returnedDesign.Properties.VariableNames));
    returnedMatchesSweep = returnedHasSweep && ...
        isequaln(returnedDesign{:,sweepNames},designTable{:,:});
    if ~returnedMatchesSweep
        error('run_ride_height_stiffness_sweep:designMismatch', ...
            'carConfig did not retain the requested sweep override values.');
    end

    caseResults = doeRunBatch(carCell,eventParams,study, ...
        (1:height(designTable)).');
    result = assembleResult(designTable,caseResults,eventParams, ...
        mapInfo,mapPreflight,settings,options,totalTimer);

    result.output = writeOutputs(result,options);

    logDetails = sprintf(['status=complete; mat=%s; csv=%s; ' ...
        'figure=%s'],result.output.matPath,result.output.csvPath, ...
        result.output.figurePath);
    simLog.finish(job,'events',cellstr(options.events), ...
        'workers',effectiveWorkers,'nCases',height(designTable), ...
        'details',logDetails);
catch ME
    try
        simLog.finish(job,'events',cellstr(options.events), ...
            'workers',effectiveWorkers,'nCases',height(designTable), ...
            'details',"status=failed; " + string(ME.identifier));
    catch
    end
    rethrow(ME)
end
end

function options = normalizeOptions(options,fullCarModelsDir,repoRoot)
if ~isstruct(options) || numel(options) ~= 1
    error('run_ride_height_stiffness_sweep:badOptions', ...
        'options must be a scalar struct.');
end

options.preview = getOption(options,'preview',false);
options.bootstrap = getOption(options,'bootstrap',true);
options.aeroMapPath = getOption(options,'aeroMapPath', ...
    fullfile(fullCarModelsDir,'aeromap_b26.csv'));
options.numRideLevels = getOption(options,'numRideLevels',5);
options.frontRideHeights_in = getOption(options,'frontRideHeights_in',[]);
options.rearRideHeights_in = getOption(options,'rearRideHeights_in',[]);
options.frontSpringRates_lb_in = getOption(options, ...
    'frontSpringRates_lb_in',[250 300 350]);
options.rearSpringRates_lb_in = getOption(options, ...
    'rearSpringRates_lb_in',[200 250 300]);
options.events = string(getOption(options,'events', ...
    ["skidpad","accel","autocross"]));
options.numWorkers = getOption(options,'numWorkers',8);
options.allowSerialFallback = getOption(options,'allowSerialFallback',true);
options.requireMapCoverage = getOption(options,'requireMapCoverage',false);
options.outputDirectory = getOption(options,'outputDirectory', ...
    fullfile(fullCarModelsDir,'sweeps'));
hasOutputPath = isfield(options,'outputPath') && ~isempty(options.outputPath);
if ~hasOutputPath
    options.outputPath = fullfile(options.outputDirectory, ...
        'ride_height_stiffness_sweep.mat');
end
options.overwrite = getOption(options,'overwrite',false);
options.saveRawResults = getOption(options,'saveRawResults',false);
options.plot = getOption(options,'plot',true);
options.verbose = getOption(options,'verbose',true);

options.aeroMapPath = char(string(options.aeroMapPath));
options.outputDirectory = char(string(options.outputDirectory));
options.outputPath = char(string(options.outputPath));
if ~isfile(options.aeroMapPath)
    error('run_ride_height_stiffness_sweep:aeroMapNotFound', ...
        'Aero-map file not found: %s',options.aeroMapPath);
end

validateattributes(options.numRideLevels,{'numeric'}, ...
    {'scalar','integer','>=',2},mfilename,'options.numRideLevels');
validateattributes(options.numWorkers,{'numeric'}, ...
    {'scalar','integer','nonnegative'},mfilename,'options.numWorkers');
options.allowSerialFallback = scalarLogical(options.allowSerialFallback, ...
    'allowSerialFallback');
options.requireMapCoverage = scalarLogical(options.requireMapCoverage, ...
    'requireMapCoverage');
options.overwrite = scalarLogical(options.overwrite,'overwrite');
options.saveRawResults = scalarLogical(options.saveRawResults,'saveRawResults');
options.plot = scalarLogical(options.plot,'plot');
options.verbose = scalarLogical(options.verbose,'verbose');

allowedEvents = ["skidpad","accel","autocross","endurance"];
if isempty(options.events) || any(~ismember(options.events,allowedEvents)) || ...
        numel(unique(options.events)) ~= numel(options.events)
    error('run_ride_height_stiffness_sweep:badEvents', ...
        'events must contain unique names from skidpad, accel, autocross, endurance.');
end
requiredEvents = ["skidpad","accel","autocross"];
if ~all(ismember(requiredEvents,options.events))
    error('run_ride_height_stiffness_sweep:missingMetricsEvents', ...
        'events must include skidpad, accel, and autocross.');
end

if ~isfolder(options.outputDirectory) && ~options.preview
    mkdir(options.outputDirectory);
end
if isempty(fileparts(options.outputPath))
    options.outputPath = fullfile(options.outputDirectory,options.outputPath);
end

% Keep the repository root visible in the saved settings without making the
% runtime depend on a caller's current folder.
options.repoRoot = repoRoot;
end

function [designTable,settings] = makeDesign(map,options)
front = axisFromOptions(options.frontRideHeights_in,options.numRideLevels, ...
    map.frontRangeIn,'frontRideHeights_in');
rear = axisFromOptions(options.rearRideHeights_in,options.numRideLevels, ...
    map.rearRangeIn,'rearRideHeights_in');
frontSpring = positiveAxis(options.frontSpringRates_lb_in, ...
    'frontSpringRates_lb_in');
rearSpring = positiveAxis(options.rearSpringRates_lb_in, ...
    'rearSpringRates_lb_in');

validateWithinBounds(front,map.frontRangeIn,'frontRideHeights_in');
validateWithinBounds(rear,map.rearRangeIn,'rearRideHeights_in');

[frontGrid,rearGrid,frontSpringGrid,rearSpringGrid] = ndgrid( ...
    front,rear,frontSpring,rearSpring);
designTable = table( ...
    frontGrid(:),rearGrid(:),frontSpringGrid(:),rearSpringGrid(:), ...
    'VariableNames',{ ...
        'static_front_ride_height_in', ...
        'static_rear_ride_height_in', ...
        'spring_rate_front_lb_in', ...
        'spring_rate_rear_lb_in'});

settings = struct( ...
    'frontRideHeights_in',front, ...
    'rearRideHeights_in',rear, ...
    'frontSpringRates_lb_in',frontSpring, ...
    'rearSpringRates_lb_in',rearSpring, ...
    'events',options.events, ...
    'requireMapCoverage',options.requireMapCoverage, ...
    'saveRawResults',options.saveRawResults, ...
    'aeroMapPath',options.aeroMapPath, ...
    'repoRoot',options.repoRoot, ...
    'outputDirectory',options.outputDirectory, ...
    'outputPath',options.outputPath);
end

function axis = axisFromOptions(explicitValues,numLevels,bounds,name)
if isempty(explicitValues)
    axis = linspace(bounds(1),bounds(2),numLevels).';
else
    axis = numericAxis(explicitValues,name);
end
axis = unique(sort(axis(:)));
end

function axis = positiveAxis(values,name)
axis = numericAxis(values,name);
if any(axis <= 0)
    error('run_ride_height_stiffness_sweep:invalidSpringRate', ...
        '%s must contain positive spring rates.',name);
end
end

function axis = numericAxis(values,name)
if ~isnumeric(values) || isempty(values) || any(~isfinite(values(:)))
    error('run_ride_height_stiffness_sweep:badAxis', ...
        '%s must contain finite numeric values.',name);
end
axis = double(values(:));
end

function validateWithinBounds(axis,bounds,name)
tolerance = 1e-10*max(1,max(abs(bounds)));
if any(axis < bounds(1)-tolerance) || any(axis > bounds(2)+tolerance)
    error('run_ride_height_stiffness_sweep:rideHeightOutsideAeroMap', ...
        '%s must stay within [%g, %g] in; requested values include an out-of-map height.', ...
        name,bounds(1),bounds(2));
end
end

function preflight = preflightMap(map,designTable,options)
[~,~,~,~,outsideMap,coverageValid] = map.evaluateNumeric( ...
    designTable.static_front_ride_height_in, ...
    designTable.static_rear_ride_height_in);
outsideMap = logical(outsideMap(:));
coverageValid = logical(coverageValid(:));
if any(outsideMap)
    error('run_ride_height_stiffness_sweep:rideHeightOutsideAeroMap', ...
        'The generated design contains static ride heights outside the aero-map envelope.');
end
if options.requireMapCoverage && any(~coverageValid)
    error('run_ride_height_stiffness_sweep:rideHeightOutsideAeroCoverage', ...
        ['The generated design contains points outside the aero-map scattered ' ...
         'coverage. Set requireMapCoverage=false to retain them as diagnostics.']);
end
preflight = struct( ...
    'outsideMap',outsideMap, ...
    'coverageValid',coverageValid, ...
    'outsideMapCount',nnz(outsideMap), ...
    'unsupportedCoverageCount',nnz(~coverageValid));
end

function mapInfo = makeMapInfo(map)
mapInfo = struct( ...
    'sourcePath',char(map.sourcePath), ...
    'frontRangeIn',double(map.frontRangeIn), ...
    'rearRangeIn',double(map.rearRangeIn), ...
    'sampleCount',numel(map.sampleFrontOffsetIn), ...
    'interpolationMode',char(map.interpolationMode));
end

function workers = resolveWorkers(requested,allowSerialFallback)
workers = requested;
if workers == 0
    return
end
hasParallel = false;
try
    hasParallel = license('test','Distrib_Computing_Toolbox');
catch
end
if hasParallel
    return
end
if allowSerialFallback
    warning('run_ride_height_stiffness_sweep:serialFallback', ...
        'Parallel Computing Toolbox is unavailable; running serially.');
    workers = 0;
else
    error('run_ride_height_stiffness_sweep:noParallelToolbox', ...
        'Parallel Computing Toolbox is required for numWorkers > 0.');
end
end

function result = assembleResult(designTable,caseResults,eventParams, ...
        mapInfo,mapPreflight,settings,options,totalTimer)
n = numel(caseResults);
metricCells = cell(n,1);
points = cell(n,1);
status = strings(n,1);
elapsed = nan(n,1);
errorIdentifier = strings(n,1);
errorMessage = strings(n,1);
for i = 1:n
    one = caseResults{i};
    metricCells{i} = one.metricRow;
    points{i} = one.points;
    status(i) = string(one.status);
    elapsed(i) = double(one.elapsed);
    errorIdentifier(i) = string(one.errorIdentifier);
    errorMessage(i) = string(one.errorMessage);
end
metricTable = vertcat(metricCells{:});
caseTable = designTable;
caseTable.case_index = (1:n).';
caseTable.status = status;
caseTable.elapsed_s = elapsed;
caseTable.error_identifier = errorIdentifier;
caseTable.error_message = errorMessage;
caseTable.static_map_coverage_valid = mapPreflight.coverageValid;

result = struct( ...
    'designTable',designTable, ...
    'caseTable',caseTable, ...
    'metricTable',metricTable, ...
    'points',{points}, ...
    'status',status, ...
    'total_elapsed_s',toc(totalTimer), ...
    'eventParams',eventParams, ...
    'mapInfo',mapInfo, ...
    'mapPreflight',mapPreflight, ...
    'settings',settings, ...
    'output',emptyOutput(), ...
    'rawResults',[]);
if options.saveRawResults
    result.rawResults = caseResults;
end
end

function output = writeOutputs(result,options)
output = emptyOutput();
if ~options.overwrite
    matPath = uniqueOutputPath(options.outputPath,options.outputDirectory, ...
        'ride_height_stiffness_sweep');
else
    matPath = normalizedMatPath(options.outputPath,options.outputDirectory, ...
        'ride_height_stiffness_sweep');
end
[folder,base,~] = fileparts(matPath);
csvPath = fullfile(folder,[base '_summary.csv']);
figurePath = fullfile(folder,[base '_summary.png']);
if ~options.overwrite
    [csvPath,figurePath] = uniqueCompanionPaths(csvPath,figurePath,folder,base);
end
if ~isfolder(folder), mkdir(folder); end

if options.plot
    writeSummaryFigure(caseTable,metricTable,char(figurePath));
else
    figurePath = '';
end

output.matPath = char(matPath);
output.csvPath = char(csvPath);
output.figurePath = char(figurePath);

designTable = result.designTable; %#ok<NASGU>
caseTable = result.caseTable; %#ok<NASGU>
metricTable = result.metricTable; %#ok<NASGU>
points = result.points; %#ok<NASGU>
settings = result.settings; %#ok<NASGU>
mapInfo = result.mapInfo; %#ok<NASGU>
mapPreflight = result.mapPreflight; %#ok<NASGU>
eventParams = result.eventParams; %#ok<NASGU>
status = result.status; %#ok<NASGU>
rawResults = result.rawResults; %#ok<NASGU>
summaryTable = mergeSummaryTable(caseTable,metricTable); %#ok<NASGU>
sweepResult = result;
sweepResult.output = output; %#ok<NASGU>
save(char(matPath),'sweepResult','designTable','caseTable','metricTable', ...
    'summaryTable','points','settings','mapInfo','mapPreflight','eventParams', ...
    'status','rawResults','-v7.3');
writetable(summaryTable,char(csvPath));
end

function summaryTable = mergeSummaryTable(caseTable,metricTable)
summaryTable = caseTable;
metricNames = metricTable.Properties.VariableNames;
summaryNames = summaryTable.Properties.VariableNames;
for i = 1:numel(metricNames)
    name = metricNames{i};
    if ~ismember(name,summaryNames)
        summaryTable.(name) = metricTable.(name);
    end
end
end

function path = normalizedMatPath(requested,outputDirectory,defaultBase)
[folder,base,ext] = fileparts(char(requested));
if isempty(folder), folder = outputDirectory; end
if isempty(base), base = defaultBase; end
if isempty(ext), ext = '.mat'; end
path = fullfile(folder,[base ext]);
end

function path = uniqueOutputPath(requested,outputDirectory,defaultBase)
path = normalizedMatPath(requested,outputDirectory,defaultBase);
if ~any(isfile(path) | isfolder(path))
    return
end
[folder,base,ext] = fileparts(path);
for k = 1:100000
    candidate = fullfile(folder,sprintf('%s_%03d%s',base,k,ext));
    [~,candidateBase,~] = fileparts(candidate);
    companions = string({ ...
        fullfile(folder,[candidateBase '_summary.csv']), ...
        fullfile(folder,[candidateBase '_summary.png'])});
    if ~isfile(candidate) && ~any(isfile(companions))
        path = candidate;
        return
    end
end
error('run_ride_height_stiffness_sweep:noUniqueOutputPath', ...
    'Could not find a unique output path near %s.',path);
end

function [csvPath,figurePath] = uniqueCompanionPaths(csvPath,figurePath,folder,base)
if ~any(isfile(string({csvPath,figurePath})))
    return
end
for k = 1:100000
    candidateBase = sprintf('%s_%03d',base,k);
    candidateCsv = fullfile(folder,[candidateBase '_summary.csv']);
    candidateFigure = fullfile(folder,[candidateBase '_summary.png']);
    if ~any(isfile(string({candidateCsv,candidateFigure})))
        csvPath = candidateCsv;
        figurePath = candidateFigure;
        return
    end
end
error('run_ride_height_stiffness_sweep:noUniqueOutputPath', ...
    'Could not find unique companion output paths near %s.',base);
end

function writeSummaryFigure(caseTable,metricTable,path)
fig = figure('Visible','off','Color','w');
cleanup = onCleanup(@() close(fig)); %#ok<NASGU>
tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');
plotMetric(nexttile,'Autocross time (s)',caseTable,metricTable,'t_autox');
plotMetric(nexttile,'Peak lateral g',caseTable,metricTable,'gLat_peak_g');
plotMetric(nexttile,'g-g coverage',caseTable,metricTable,'gg_coverage');
exportgraphics(fig,path,'Resolution',150);
end

function plotMetric(ax,titleText,caseTable,metricTable,fieldName)
values = metricTable.(fieldName);
front = caseTable.static_front_ride_height_in;
rear = caseTable.static_rear_ride_height_in;
valid = isfinite(values) & isfinite(front) & isfinite(rear);
if any(valid)
    scatter(ax,front(valid),rear(valid),32,values(valid),'filled');
    colorbar(ax);
else
    text(ax,0.5,0.5,'No finite data','HorizontalAlignment','center');
end
title(ax,titleText,'Interpreter','none');
xlabel(ax,'Front ride height offset (in)');
ylabel(ax,'Rear ride height offset (in)');
grid(ax,'on');
end

function output = emptyOutput()
output = struct('matPath','','csvPath','','figurePath','');
end

function value = getOption(options,name,defaultValue)
if isfield(options,name) && ~isempty(options.(name))
    value = options.(name);
else
    value = defaultValue;
end
end

function value = scalarLogical(value,name)
if ~(islogical(value) || isnumeric(value)) || ~isscalar(value) || ...
        ~isfinite(double(value))
    error('run_ride_height_stiffness_sweep:badLogical', ...
        '%s must be a scalar logical value.',name);
end
value = logical(value);
end

