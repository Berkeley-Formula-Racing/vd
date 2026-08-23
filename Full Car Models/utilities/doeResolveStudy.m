function resolved = doeResolveStudy(study,baselineTable)
%DOERESOLVESTUDY Validate a DOE study and resolve its physical design bounds.

validateStudy(study);
P = validateParameters(study.parameters,baselineTable);

payload = struct('name',{cellstr(P.name)},'baseline',P.baseline, ...
    'rangeType',{cellstr(P.rangeType)}, ...
    'lower',P.lowerPhysical,'upper',P.upperPhysical);

resolved = struct();
resolved.parameters = P;
resolved.toPhysical = @(U) toPhysical(U,P);
resolved.signature = string(jsonencode(payload));
end

function P = validateParameters(parameters,baselineTable)
required = {'name','lower','upper','rangeType'};
if ~istable(parameters) || ~all(ismember(required,parameters.Properties.VariableNames))
    error('doeResolveStudy:badParameters', ...
        'study.parameters must be a table with name, lower, upper, and rangeType columns.')
end
if height(parameters) == 0
    error('doeResolveStudy:badParameters','study.parameters must not be empty.')
end

name = string(parameters.name(:));
rangeType = lower(string(parameters.rangeType(:)));
lowerInput = parameters.lower(:);
upperInput = parameters.upper(:);
if any(ismissing(name) | strlength(name) == 0) || numel(unique(name)) ~= numel(name)
    error('doeResolveStudy:duplicateParameter', ...
        'Each study parameter name must be present exactly once.')
end
if ~isnumeric(lowerInput) || ~isnumeric(upperInput) || ...
        ~isreal(lowerInput) || ~isreal(upperInput) || ...
        any(~isfinite(lowerInput)) || any(~isfinite(upperInput))
    error('doeResolveStudy:badBounds','Parameter bounds must be finite real numbers.')
end
if any(~ismember(rangeType,["percent","absolute"]))
    error('doeResolveStudy:badRangeType', ...
        'Parameter rangeType values must be percent or absolute.')
end
if any(lowerInput >= upperInput)
    error('doeResolveStudy:badBounds', ...
        'Each parameter lower bound must be less than its upper bound.')
end

knownNames = parameterRegistry();
if any(~ismember(name,knownNames))
    error('doeResolveStudy:unknownParameter', ...
        'Every study parameter must be listed by parameters_loop.m.')
end
if ~istable(baselineTable) || height(baselineTable) ~= 1
    error('doeResolveStudy:badBaseline', ...
        'baselineTable must be a one-row table.')
end
if any(~ismember(cellstr(name),baselineTable.Properties.VariableNames))
    error('doeResolveStudy:missingBaselineParameter', ...
        'baselineTable must contain every requested study parameter.')
end

baseline = zeros(numel(name),1);
for i = 1:numel(name)
    value = baselineTable.(char(name(i)));
    if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ~isfinite(value)
        error('doeResolveStudy:badBaseline', ...
            'Baseline values must be finite real numeric scalars.')
    end
    baseline(i) = value;
end
isPercent = rangeType == "percent";
if any(isPercent & baseline <= 0)
    error('doeResolveStudy:badPercentBaseline', ...
        'Percentage ranges require positive, nonzero baseline values.')
end

lowerPhysical = lowerInput;
upperPhysical = upperInput;
lowerPhysical(isPercent) = baseline(isPercent) .* (1 + lowerInput(isPercent)./100);
upperPhysical(isPercent) = baseline(isPercent) .* (1 + upperInput(isPercent)./100);
if any(isPercent & (lowerPhysical <= 0 | upperPhysical <= 0))
    error('doeResolveStudy:badPercentPhysicalBounds', ...
        'Percentage ranges must resolve to strictly positive physical bounds.')
end
if any(~isfinite(lowerPhysical) | ~isfinite(upperPhysical) | ...
        lowerPhysical >= upperPhysical)
    error('doeResolveStudy:badBounds', ...
        'Resolved physical bounds must be finite and ordered.')
end

P = table(name,rangeType,lowerInput,upperInput,baseline,lowerPhysical,upperPhysical, ...
    'VariableNames',{'name','rangeType','lowerInput','upperInput','baseline', ...
    'lowerPhysical','upperPhysical'});
end

function knownNames = parameterRegistry()
registryFile = which('parameters_loop');
if isempty(registryFile)
    registryFile = fullfile(fileparts(mfilename('fullpath')),'parameters_loop.m');
end
source = fileread(registryFile);
match = regexp(source,'names\s*=\s*\{([\s\S]*?)\};','tokens','once');
if isempty(match)
    error('doeResolveStudy:registryUnavailable', ...
        'Could not read the parameter registry from parameters_loop.m.')
end
tokens = regexp(match{1},'''([^'']+)''','tokens');
knownNames = string(cellfun(@(token) token{1},tokens,'UniformOutput',false))';
end

function physical = toPhysical(U,P)
if ~isnumeric(U) || ~isreal(U) || ndims(U) ~= 2 || size(U,2) ~= height(P) || ...
        any(~isfinite(U),'all') || any(U < 0 | U > 1,'all')
    error('doeResolveStudy:badNormalizedInput', ...
        'U must be a finite N-by-P numeric matrix with values in [0,1].')
end
X = P.lowerPhysical' + U .* (P.upperPhysical' - P.lowerPhysical');
physical = array2table(X,'VariableNames',cellstr(P.name)');
end

function validateStudy(study)
if ~isstruct(study) || ~isscalar(study) || ~isfield(study,'parameters')
    error('doeResolveStudy:badStudy', ...
        'study must be a scalar struct with a parameters table.')
end

validateScalarField(study,'name',@(x) (ischar(x) || isstring(x)) && ...
    isscalar(x) && strlength(string(x)) > 0,'badName');
validateScalarField(study,'mode',@(x) ...
    (isstring(x) && isscalar(x)) || (ischar(x) && isrow(x)),'badMode');
if isfield(study,'mode') && ~ismember(lower(string(study.mode)), ...
        ["sensitivity","optimization","hybrid"])
    error('doeResolveStudy:badMode', ...
        'study.mode must be sensitivity, optimization, or hybrid.')
end
validateIntegerField(study,'randomSeed',0);
validateIntegerField(study,'initialCases',1);
validateIntegerField(study,'batchSize',1);
validateIntegerField(study,'maxCases',1);
validateIntegerField(study,'numWorkers',1);
validateLogicalField(study,'allowSerialFallback');
validateLogicalField(study,'resume');

if all(isfield(study,{'initialCases','batchSize','maxCases'})) && ...
        (study.initialCases > study.maxCases || study.batchSize > study.maxCases)
    error('doeResolveStudy:badCaseCounts', ...
        'initialCases and batchSize must not exceed maxCases.')
end
validateAdaptive(study);
validateEvents(study);
validateRamps(study);
end

function validateScalarField(study,fieldName,predicate,errorSuffix)
if isfield(study,fieldName) && ~predicate(study.(fieldName))
    error(['doeResolveStudy:' errorSuffix], ...
        'study.%s must be a valid scalar value.',fieldName)
end
end

function validateIntegerField(study,fieldName,minimum)
if isfield(study,fieldName)
    value = study.(fieldName);
    if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value) || value < minimum || value ~= floor(value)
        error('doeResolveStudy:badScalarSetting', ...
            'study.%s must be an integer greater than or equal to %d.',fieldName,minimum)
    end
end
end

function validateLogicalField(study,fieldName)
if isfield(study,fieldName) && (~islogical(study.(fieldName)) || ~isscalar(study.(fieldName)))
    error('doeResolveStudy:badScalarSetting', ...
        'study.%s must be a logical scalar.',fieldName)
end
end

function validateAdaptive(study)
if ~isfield(study,'adaptive'), return, end
adaptive = study.adaptive;
if ~isstruct(adaptive) || ~isscalar(adaptive)
    error('doeResolveStudy:badAdaptive','study.adaptive must be a scalar struct.')
end
if isfield(adaptive,'responses')
    responses = string(adaptive.responses);
    valid = ["modeled_dynamic_points","t_autox","t_accel","t_skid", ...
        "total_work_kJ","understeer_gradient_10_deg_per_g", ...
        "understeer_gradient_25_deg_per_g"];
    if ~isvector(responses) || any(~ismember(responses,valid))
        error('doeResolveStudy:badAdaptiveResponse', ...
            'study.adaptive.responses contains an unsupported response name.')
    end
end
validateAdaptiveNumber(adaptive,'hybridSensitivityFraction',0,1);
validateAdaptiveNumber(adaptive,'candidatePoolSize',1,Inf,true);
validateAdaptiveNumber(adaptive,'minimumDistance',0,Inf);
end

function validateAdaptiveNumber(adaptive,fieldName,lower,upper,integerOnly)
if nargin < 5, integerOnly = false; end
if isfield(adaptive,fieldName)
    value = adaptive.(fieldName);
    if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value) || ...
            value < lower || value > upper || (integerOnly && value ~= floor(value))
        error('doeResolveStudy:badAdaptive', ...
            'study.adaptive.%s is outside its valid range.',fieldName)
    end
end
end

function validateEvents(study)
if ~isfield(study,'events'), return, end
events = string(study.events(:));
scored = ["skidpad","accel","autocross","endurance"];
if any(~ismember(events,scored)) || numel(unique(events)) ~= numel(events)
    error('doeResolveStudy:badEvents','study.events contains an invalid or duplicate event.')
end
if isfield(study,'mode') && ismember(lower(string(study.mode)),["optimization","hybrid"]) && ...
        ~all(ismember(scored,events))
    error('doeResolveStudy:missingScoredEvents', ...
        'Optimization and hybrid studies require all four scored events.')
end
if isfield(study,'adaptive') && isfield(study.adaptive,'responses') && ...
        any(string(study.adaptive.responses) == "modeled_dynamic_points") && ...
        ~all(ismember(scored,events))
    error('doeResolveStudy:missingPointEvents', ...
        'modeled_dynamic_points requires all four scored events.')
end
end

function validateRamps(study)
if ~isfield(study,'ramps'), return, end
ramps = study.ramps;
if ~isstruct(ramps) || ~isscalar(ramps)
    error('doeResolveStudy:badRamps','study.ramps must be a scalar struct.')
end
if isfield(ramps,'enabled') && (~islogical(ramps.enabled) || ~isscalar(ramps.enabled))
    error('doeResolveStudy:badRamps','study.ramps.enabled must be a logical scalar.')
end
for fieldName = ["nRamp","nBisect"]
    if isfield(ramps,fieldName)
        value = ramps.(fieldName);
        if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value) || value < 0 || value ~= floor(value)
            error('doeResolveStudy:badRamps','study.ramps.%s must be a nonnegative integer.',fieldName)
        end
    end
end
end
