function [Unew,selection] = doeSelectAdaptiveBatch(state,study,candidateU)
%DOESELECTADAPTIVEBATCH Select a diverse normalized batch by DOE acquisition.

validateState(state);
mode = lower(string(study.mode));
if ~ismember(mode,["sensitivity","optimization","hybrid"])
    error('doeSelectAdaptiveBatch:badMode', ...
        'study.mode must be sensitivity, optimization, or hybrid.')
end
validateattributes(study.batchSize,{'numeric'},{'scalar','integer','positive'})
if ~isfield(study,'adaptive') || ~isstruct(study.adaptive)
    error('doeSelectAdaptiveBatch:badStudy','study.adaptive is required.')
end
adaptive = study.adaptive;
for required = ["responses","minimumDistance"]
    if ~isfield(adaptive,required)
        error('doeSelectAdaptiveBatch:badStudy', ...
            'study.adaptive.%s is required.',required)
    end
end
validateattributes(adaptive.minimumDistance,{'numeric'}, ...
    {'scalar','real','finite','nonnegative'})

if nargin < 3 || isempty(candidateU)
    [candidateU,nextRandomState] = generateCandidatePool(state,study);
else
    candidateU = validateCandidates(candidateU,size(state.U,2));
    nextRandomState = [];
end
candidateU = removeCompletedAndDuplicateRows(candidateU,state.U);
if size(candidateU,1) < study.batchSize
    error('doeSelectAdaptiveBatch:insufficientCandidates', ...
        'Fewer than study.batchSize unused candidates are available.')
end

models = doeFitSurrogates(state.U,state.metricTable,adaptive.responses);
sensitivityScore = sensitivityAcquisition(models,candidateU);
optimizationScore = optimizationAcquisition(models,state.metricTable,candidateU);

nSensitivity = 0;
nOptimization = 0;
switch mode
    case "sensitivity"
        nSensitivity = study.batchSize;
    case "optimization"
        nOptimization = study.batchSize;
    case "hybrid"
        if ~isfield(adaptive,'hybridSensitivityFraction')
            error('doeSelectAdaptiveBatch:badStudy', ...
                'study.adaptive.hybridSensitivityFraction is required for hybrid mode.')
        end
        validateattributes(adaptive.hybridSensitivityFraction,{'numeric'}, ...
            {'scalar','real','finite','>=',0,'<=',1})
        nSensitivity = round(study.batchSize * adaptive.hybridSensitivityFraction);
        nOptimization = study.batchSize - nSensitivity;
end

[picked,source] = diverseSelection(candidateU,state.U,sensitivityScore, ...
    optimizationScore,nSensitivity,nOptimization,adaptive.minimumDistance);
Unew = candidateU(picked,:);
nearest = nearestDistances(Unew,state.U);
acquisition = zeros(numel(picked),1);
isSensitivity = source == "sensitivity";
acquisition(isSensitivity) = sensitivityScore(picked(isSensitivity));
acquisition(~isSensitivity) = optimizationScore(picked(~isSensitivity));
selection = table(repmat(mode,numel(picked),1),source,acquisition,nearest,picked, ...
    'VariableNames',{'mode','source','acquisition_value', ...
    'nearest_existing_distance','candidate_index'});
selection.Properties.UserData = struct('randomState',nextRandomState);
end

function validateState(state)
if ~isstruct(state) || ~isfield(state,'U') || ~isfield(state,'metricTable')
    error('doeSelectAdaptiveBatch:badState', ...
        'state must contain U and metricTable.')
end
validateattributes(state.U,{'numeric'},{'2d','real','finite','nonempty','>=',0,'<=',1})
if ~istable(state.metricTable) || height(state.metricTable) ~= size(state.U,1)
    error('doeSelectAdaptiveBatch:badState', ...
        'state.metricTable must align with state.U.')
end
end

function [candidateU,nextRandomState] = generateCandidatePool(state,study)
if ~isfield(study.adaptive,'candidatePoolSize')
    error('doeSelectAdaptiveBatch:badStudy', ...
        'study.adaptive.candidatePoolSize is required when candidateU is omitted.')
end
validateattributes(study.adaptive.candidatePoolSize,{'numeric'}, ...
    {'scalar','integer','positive'})
previousRng = rng;
restoreRng = onCleanup(@() rng(previousRng)); %#ok<NASGU>
if isfield(state,'randomState') && isstruct(state.randomState)
    rng(state.randomState);
elseif isfield(study,'randomSeed')
    rng(study.randomSeed + size(state.U,1),'twister');
end
candidateU = rand(study.adaptive.candidatePoolSize,size(state.U,2));
nextRandomState = rng;
end

function candidateU = validateCandidates(candidateU,nVar)
validateattributes(candidateU,{'numeric'}, ...
    {'2d','real','finite','ncols',nVar,'>=',0,'<=',1})
end

function C = removeCompletedAndDuplicateRows(C,U)
if isempty(C), return, end
isCompleted = false(size(C,1),1);
for i = 1:size(C,1)
    isCompleted(i) = any(all(abs(U-C(i,:)) <= 1e-12,2));
end
C = C(~isCompleted,:);
[~,first] = unique(C,'rows','stable');
C = C(first,:);
end

function score = sensitivityAcquisition(models,C)
names = string(models.responseNames(:));
scaledVariance = zeros(size(C,1),0);
for i = 1:numel(names)
    fieldName = matlab.lang.makeValidName(char(names(i)));
    if ~isfield(models.responses,fieldName), continue, end
    [~,sd] = predict(models.responses.(fieldName),C);
    scaledVariance(:,end+1) = (sd ./ models.responseSpreads.(fieldName)).^2; %#ok<AGROW>
end
score = sqrt(mean(scaledVariance,2));
end

function score = optimizationAcquisition(models,metricTable,C)
if isempty(models.objective)
    error('doeSelectAdaptiveBatch:noObjective', ...
        'Optimization requires a finite, nonconstant objective_score.')
end
objective = metricTable.objective_score;
valid = metricTable.valid;
ok = valid & isnumeric(objective) & isfinite(objective);
if ~any(ok)
    error('doeSelectAdaptiveBatch:noObjective', ...
        'Optimization requires at least one valid finite objective score.')
end
[mu,sd] = predict(models.objective,C);
score = doeExpectedImprovement(mu,sd,max(objective(ok))) .* ...
    feasibilityProbability(models.feasibility,C);
end

function probability = feasibilityProbability(feasibility,C)
if isnumeric(feasibility)
    probability = repmat(feasibility,size(C,1),1);
    return
end
[~,score] = predict(feasibility,C);
classNames = string(feasibility.ClassNames);
positive = find(classNames == "true" | classNames == "1",1);
if isempty(positive)
    error('doeSelectAdaptiveBatch:badFeasibility', ...
        'The feasibility classifier does not have a true class.')
end
if any(score(:) < 0) || any(abs(sum(score,2)-1) > 1e-8)
    shifted = score - max(score,[],2);
    score = exp(shifted);
    score = score ./ sum(score,2);
end
probability = score(:,positive);
end

function [picked,source] = diverseSelection(C,U,sensitivity,optimization,nSensitivity,nOptimization,minimumDistance)
threshold = minimumDistance;
while true
    picked = zeros(0,1);
    source = strings(0,1);
    [picked,source] = appendRanked(C,U,picked,source,sensitivity, ...
        "sensitivity",nSensitivity,threshold);
    [picked,source] = appendRanked(C,U,picked,source,optimization, ...
        "optimization",nOptimization,threshold);
    if numel(picked) == nSensitivity + nOptimization
        return
    end
    if threshold <= eps
        error('doeSelectAdaptiveBatch:insufficientDiverseCandidates', ...
            'The candidate pool cannot fill the requested batch.')
    end
    threshold = threshold / 2;
end
end

function [picked,source] = appendRanked(C,U,picked,source,score,label,count,threshold)
if count == 0, return, end
ranked = rankCandidates(score);
for i = 1:numel(ranked)
    candidate = ranked(i);
    if ismember(candidate,picked), continue, end
    reference = [U; C(picked,:)];
    if all(nearestDistances(C(candidate,:),reference) >= threshold)
        picked(end+1,1) = candidate; %#ok<AGROW>
        source(end+1,1) = label; %#ok<AGROW>
        if sum(source == label) == count, return, end
    end
end
end

function ranked = rankCandidates(score)
score = score(:);
score(~isfinite(score)) = -Inf;
[~,ranked] = sortrows([-score,(1:numel(score))'],[1 2]);
end

function distances = nearestDistances(points,reference)
if isempty(reference)
    distances = inf(size(points,1),1);
    return
end
distances = inf(size(points,1),1);
for i = 1:size(points,1)
    distances(i) = min(sqrt(mean((reference-points(i,:)).^2,2)));
end
end
