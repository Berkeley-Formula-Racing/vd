function result = exportStudy(study,outputDirectory,options)
%EXPORTSTUDY Export a validated schema-v1 study and its report tables.

if nargin < 2 || isempty(outputDirectory)
    outputDirectory = pwd;
end
if nargin < 3 || isempty(options)
    options = struct();
end
if ~isstruct(options) || ~isscalar(options)
    error("rampSpeed:invalidExportOptions", ...
        "options must be a scalar struct.");
end

[ok,issues] = rampSpeed.validateStudy(study);
if ~ok
    error("rampSpeed:invalidStudy", ...
        "Invalid ramp-speed study: %s",joinIssues(issues));
end

outputDirectory = ensureOutputDirectory(outputDirectory);
baseName = exportBaseName(options);
[perSpeed,points] = makeExportTables(study);

result = struct();
result.outputDirectory = outputDirectory;
result.baseName = baseName;
result.perSpeed = perSpeed;
result.points = points;
result.perSpeedCsv = absolutePath(fullfile(outputDirectory, ...
    baseName + "_per_speed.csv"));
result.pointsCsv = absolutePath(fullfile(outputDirectory, ...
    baseName + "_points.csv"));
writetable(perSpeed,result.perSpeedCsv);
writetable(points,result.pointsCsv);

result.studyMat = absolutePath(fullfile(outputDirectory,baseName + "_study.mat"));
rampSpeed.saveStudy(result.studyMat,study);

result.terminalLog = makeTerminalLog(study);
result.terminalLogCsv = absolutePath(fullfile(outputDirectory, ...
    baseName + "_terminal_log.csv"));
writetable(result.terminalLog,result.terminalLogCsv);
emitTerminalLog(result.terminalLog);

metadata = makeMetadata(study,baseName,perSpeed,points);
result.metadata = metadata;
result.metadataJson = absolutePath(fullfile(outputDirectory, ...
    baseName + "_metadata.json"));
writeJson(result.metadataJson,metadata);

result.figureFiles = exportRequestedFigures(options,outputDirectory,baseName);
result.files = [result.perSpeedCsv; result.pointsCsv; result.studyMat; ...
    result.terminalLogCsv; result.metadataJson; result.figureFiles(:)];
end

function outputDirectory = ensureOutputDirectory(outputDirectory)
outputDirectory = string(outputDirectory);
if ~isscalar(outputDirectory) || strlength(strtrim(outputDirectory)) == 0
    error("rampSpeed:invalidOutputDirectory", ...
        "outputDirectory must be a non-empty path.");
end
if ~isfolder(outputDirectory)
    [created,message,messageId] = mkdir(outputDirectory);
    if ~created
        if isempty(messageId)
            error("rampSpeed:invalidOutputDirectory", ...
                "Could not create output directory: %s",message);
        end
        error("rampSpeed:invalidOutputDirectory", ...
            "Could not create output directory (%s): %s",messageId,message);
    end
end
outputDirectory = absolutePath(outputDirectory);
end

function baseName = exportBaseName(options)
baseName = "ramp_speed_study_" + string(datetime("now", ...
    "Format","yyyyMMdd_HHmmss"));
if isfield(options,"baseName") && ~isempty(options.baseName)
    baseName = string(options.baseName);
end
if ~isscalar(baseName) || strlength(strtrim(baseName)) == 0
    error("rampSpeed:invalidExportOptions", ...
        "options.baseName must be non-empty scalar text.");
end
baseName = strtrim(baseName);
[~,name,~] = fileparts(char(baseName));
if isempty(name)
    name = char(baseName);
end
baseName = string(name);
baseName = string(regexprep(char(baseName),"[^A-Za-z0-9_.-]","_"));
if strlength(baseName) == 0
    error("rampSpeed:invalidExportOptions", ...
        "options.baseName does not contain a usable file name.");
end
end

function [perSpeed,points] = makeExportTables(study)
if isempty(study.runs)
    template = rampSpeed.makeRun("lateral","coast",struct(),struct());
    emptyRun = struct("caseId","","type","","mode","","status","");
    perSpeed = addRunContext(template.perSpeed,0,emptyRun);
    points = addRunContext(template.points,0,emptyRun);
    points = addvars(points,strings(0,1),'After','status', ...
        'NewVariableNames','point_reason');
    return
end

perSpeed = table();
points = table();
for runIndex = 1:numel(study.runs)
    run = study.runs(runIndex);
    currentPerSpeed = addRunContext(run.perSpeed,runIndex,run);
    currentPoints = addRunContext(run.points,runIndex,run);
    reasons = pointReasons(run);
    currentPoints = addvars(currentPoints,reasons,'After','status', ...
        'NewVariableNames','point_reason');
    if runIndex == 1
        perSpeed = currentPerSpeed;
        points = currentPoints;
    else
        perSpeed = [perSpeed; currentPerSpeed]; %#ok<AGROW>
        points = [points; currentPoints]; %#ok<AGROW>
    end
end
end

function T = addRunContext(T,runIndex,run)
n = height(T);
context = table(repmat(runIndex,n,1), ...
    repmat(string(run.caseId),n,1), ...
    repmat(string(run.type),n,1), ...
    repmat(string(run.mode),n,1), ...
    repmat(string(run.status),n,1), ...
    'VariableNames',{'run_index','case_id','run_type', ...
    'run_mode','run_status'});
T = [context T];
end

function reasons = pointReasons(run)
n = height(run.points);
reasons = strings(n,1);
if n == 0 || ~ismember("reason",string(run.perSpeed.Properties.VariableNames))
    return
end

if ismember("speed_index",string(run.points.Properties.VariableNames))
    indices = run.points.speed_index;
    valid = isfinite(indices) & indices >= 1 & ...
        indices <= height(run.perSpeed) & indices == floor(indices);
    indices = round(indices);
    reasons(valid) = run.perSpeed.reason(indices(valid));
    return
end

if ~ismember("speed_mps",string(run.points.Properties.VariableNames)) || ...
        ~ismember("speed_mps",string(run.perSpeed.Properties.VariableNames))
    return
end
for row = 1:n
    match = find(run.perSpeed.speed_mps == run.points.speed_mps(row),1,"first");
    if ~isempty(match)
        reasons(row) = run.perSpeed.reason(match);
    end
end
end

function terminalLog = makeTerminalLog(study)
indices = zeros(0,1);
caseIds = strings(0,1);
types = strings(0,1);
modes = strings(0,1);
statuses = strings(0,1);
messages = strings(0,1);
terminalStatuses = ["complete","completed","failed","cancelled"];

for runIndex = 1:numel(study.runs)
    run = study.runs(runIndex);
    status = lower(string(run.status));
    if ~any(status == terminalStatuses)
        continue
    end
    indices(end+1,1) = runIndex; %#ok<AGROW>
    caseIds(end+1,1) = string(run.caseId); %#ok<AGROW>
    types(end+1,1) = string(run.type); %#ok<AGROW>
    modes(end+1,1) = string(run.mode); %#ok<AGROW>
    statuses(end+1,1) = string(run.status); %#ok<AGROW>
    messages(end+1,1) = runMessage(run); %#ok<AGROW>
end

terminalLog = table(indices,caseIds,types,modes,statuses,messages, ...
    'VariableNames',{'run_index','case_id','run_type','run_mode', ...
    'status','message'});
end

function message = runMessage(run)
message = "";
if isfield(run,"runMeta") && isstruct(run.runMeta)
    message = firstTextField(run.runMeta,"errors");
    if strlength(message) == 0
        message = firstTextField(run.runMeta,"warnings");
    end
end
if strlength(message) == 0 && ismember("reason", ...
        string(run.perSpeed.Properties.VariableNames))
    message = firstNonEmpty(string(run.perSpeed.reason));
end
end

function value = firstTextField(record,name)
value = "";
if ~isfield(record,name) || isempty(record.(name))
    return
end
value = firstNonEmpty(string(record.(name)));
end

function value = firstNonEmpty(values)
values = string(values(:));
values = values(strlength(strtrim(values)) > 0);
if isempty(values)
    value = "";
else
    value = values(1);
end
end

function emitTerminalLog(terminalLog)
for row = 1:height(terminalLog)
    fprintf(1,"[rampSpeed] run %d case=%s status=%s: %s\n", ...
        terminalLog.run_index(row),terminalLog.case_id(row), ...
        terminalLog.status(row),terminalLog.message(row));
end
end

function metadata = makeMetadata(study,baseName,perSpeed,points)
metadata = struct();
metadata.schemaVersion = double(study.schemaVersion);
metadata.appVersion = char(string(study.appVersion));
metadata.created = datetimeText(study.created);
metadata.exportedAt = datetimeText(datetime("now"));
metadata.baseName = char(baseName);
metadata.canonicalUnits = canonicalUnitText();
metadata.displayUnits = unitText(study.displayUnits);
metadata.tables = struct( ...
    "perSpeed",struct("rows",height(perSpeed),"coordinate","speed_mps"), ...
    "points",struct("rows",height(points),"coordinate","speed_mps"));
metadata.validityFields = {"valid","status","reason"};
end

function units = unitText(displayUnits)
units = struct();
names = string(fieldnames(displayUnits));
for i = 1:numel(names)
    value = displayUnits.(char(names(i)));
    units.(char(names(i))) = char(string(value));
end
end

function units = canonicalUnitText()
units = struct();
units.speed = 'm/s';
units.acceleration = 'm/s^2';
units.force = 'N';
units.length = 'm';
units.angle = 'rad';
units.angularRate = 'rad/s';
end
function textValue = datetimeText(value)
if isempty(value)
    textValue = "";
    return
end
textValue = string(value);
if ~isscalar(textValue)
    textValue = textValue(1);
end
textValue = char(textValue);
end

function writeJson(fileName,metadata)
jsonText = jsonencode(metadata,"PrettyPrint",true);
fid = fopen(fileName,"w","n","UTF-8");
if fid < 0
    error("rampSpeed:exportFailed", ...
        "Could not create metadata file: %s",fileName);
end
cleanup = onCleanup(@()fclose(fid));
fprintf(fid,"%s",jsonText);
end

function figureFiles = exportRequestedFigures(options,outputDirectory,baseName)
figureFiles = strings(0,1);
requests = [];
if isfield(options,"visibleFigures") && ~isempty(options.visibleFigures)
    requests = options.visibleFigures;
elseif isfield(options,"figures") && ~isempty(options.figures)
    requests = options.figures;
end
if isempty(requests)
    return
end
if isempty(which("exportgraphics"))
    error("rampSpeed:figureExportUnavailable", ...
        "MATLAB exportgraphics is unavailable.");
end

if iscell(requests)
    requestList = requests;
elseif isstruct(requests)
    requestList = num2cell(requests);
else
    requestList = num2cell(requests);
end
for i = 1:numel(requestList)
    request = requestList{i};
    [handle,name] = figureRequest(request,i);
    if ~isgraphics(handle)
        error("rampSpeed:invalidFigure", ...
            "Figure request %d does not contain a graphics handle.",i);
    end
    fileName = fullfile(outputDirectory, ...
        baseName + "_" + name + ".png");
    fileName = absolutePath(fileName);
    exportgraphics(handle,fileName);
    figureFiles(end+1,1) = fileName; %#ok<AGROW>
end
end

function [handle,name] = figureRequest(request,index)
handle = [];
name = "figure_" + sprintf("%02d",index);
if isstruct(request)
    if isfield(request,"figure")
        handle = request.figure;
    elseif isfield(request,"handle")
        handle = request.handle;
    end
    if isfield(request,"name") && ~isempty(request.name)
        name = string(request.name);
    end
else
    handle = request;
end
name = string(regexprep(char(name),"[^A-Za-z0-9_.-]","_"));
if strlength(name) == 0
    name = "figure_" + sprintf("%02d",index);
end
end

function path = absolutePath(fileName)
path = string(char(java.io.File(char(fileName)).getAbsolutePath()));
end

function message = joinIssues(issues)
message = char(strjoin(string(issues(:)),"; "));
end
