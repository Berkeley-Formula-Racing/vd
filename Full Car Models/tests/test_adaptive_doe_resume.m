%% test_adaptive_doe_resume
% Removing persisted state reuse must rerun or change initial designs, while
% removing adaptive selection metadata must fail to mark resumed rows as optimization.
testDir = fileparts(mfilename('fullpath'));
modelRoot = fileparts(testDir);
cd(modelRoot); setup_paths

d = tempname;
mkdir(d)
cleaner = onCleanup(@() rmdir(d,'s')); %#ok<NASGU>

study = DOEStudyConfig();
study.name = "resume_smoke";
study.parameters = table("mass",-2,2,"percent", ...
    'VariableNames',{'name','lower','upper','rangeType'});
study.initialCases = 4;
study.batchSize = 2;
study.maxCases = 4;
study.numWorkers = 0;
study.mode = "sensitivity";
study.events = ["skidpad","accel","autocross","endurance"];
study.ramps.enabled = false;
study.output.directory = d;

state1 = runAdaptiveDOE(study);
assert(height(state1.designTable) == 4)
assert(size(unique(state1.U,'rows'),1) == 4)

changedEvents = study;
changedEvents.resume = true;
changedEvents.events = ["endurance","autocross","accel","skidpad"];
assertError(@() runAdaptiveDOE(changedEvents), ...
    'runAdaptiveDOE:resumeStudyMismatch')
changedAdaptive = study;
changedAdaptive.resume = true;
changedAdaptive.adaptive.minimumDistance = 0.04;
assertError(@() runAdaptiveDOE(changedAdaptive), ...
    'runAdaptiveDOE:resumeStudyMismatch')

study.mode = "optimization";
study.maxCases = 6;
study.resume = true;
state2 = runAdaptiveDOE(study);
assert(height(state2.designTable) == 6)
assert(size(unique(state2.U,'rows'),1) == 6)
assert(isequal(state2.U(1:4,:),state1.U))
assert(isequal(state2.designTable(1:4,:),state1.designTable))
assert(all(state2.selectionHistory.source(end-1:end) == "optimization"))
resultsPath = fullfile(d,study.output.results);
assert(isfile(resultsPath))
cache = load(resultsPath,'metricTable','rampData','selectionHistory');
assert(height(cache.metricTable) == 6 && numel(cache.rampData) == 6)
assert(all(cache.selectionHistory.source(end-1:end) == "optimization"))

% Changing only the objective must refresh cached scores without new cases.
objectiveStudy = study;
objectiveStudy.objective.penalties.energy.enabled = true;
objectiveStudy.objective.penalties.energy.threshold_kJ = 0;
objectiveStudy.objective.penalties.energy.pointsPerKJ = 0.1;
saveDOECheckpoint(fullfile(d,study.output.checkpoint),state2)
rescored = runAdaptiveDOE(objectiveStudy);
for i = 1:height(state2.metricTable)
    [~,expected] = doeScoreCase(state2.pointData{i},state2.metricTable(i,:), ...
        objectiveStudy.objective);
    expected = struct2table(expected);
    scoreNames = expected.Properties.VariableNames;
    assert(isequaln(rescored.metricTable{i,scoreNames},expected{1,scoreNames}))
end

% A nonempty cache from saveFullPoints=false is incompatible with true.
oldRamps = study.ramps;
oldRamps.enabled = true;
oldRamps.saveFullPoints = false;
oldSignature = string(jsonencode(orderfields(oldRamps)));
incompatible = state2;
incompatible.rampData = cell(height(state2.metricTable),1);
for i = 1:height(state2.metricTable)
    incompatible.rampData{i} = struct('summary',table(10, ...
        'VariableNames',{'vCar'}),'points',[], ...
        'configSignature',oldSignature);
    car = incompatible.carCell{i,1};
    car.comp = [];
    incompatible.carCell{i,1} = car;
end
saveDOECheckpoint(fullfile(d,study.output.checkpoint),incompatible)
rampConfigStudy = study;
rampConfigStudy.ramps.enabled = true;
rampConfigStudy.ramps.saveFullPoints = true;
incompatibleBackfill = runAdaptiveDOE(rampConfigStudy);
assertMetricTableUnchanged(incompatibleBackfill.metricTable,state2.metricTable)
assert(isequal(incompatibleBackfill.pointData,state2.pointData))
assert(isequal(incompatibleBackfill.caseStatus,state2.caseStatus))
assert(height(incompatibleBackfill.rampBackfillDiagnostics) == 6)
assert(all(incompatibleBackfill.rampBackfillDiagnostics.error_identifier == ...
    "doeRunCase:unsolvedCar"))
assert(isequaln(incompatibleBackfill.rampData,incompatible.rampData))

% A failed ramp-only pass must not discard complete full-run cache data.
metricsBefore = state2.metricTable;
pointsBefore = state2.pointData;
statusBefore = state2.caseStatus;
massBefore = cellfun(@(car) car.M,state2.carCell(:,1));
ggBefore = cellfun(@(car) car.ggPoints,state2.carCell(:,1), ...
    'UniformOutput',false);
for i = 1:size(state2.carCell,1)
    car = state2.carCell{i,1};
    car.comp = [];
    state2.carCell{i,1} = car;
end
checkpointPath = fullfile(d,study.output.checkpoint);
saveDOECheckpoint(checkpointPath,state2)
[checkpointState,~] = loadDOECheckpoint(checkpointPath,state2.resolvedStudy);
assertMetricTableUnchanged(checkpointState.metricTable,metricsBefore)

rampStudy = study;
rampStudy.ramps.enabled = true;
backfilled = runAdaptiveDOE(rampStudy);
assertMetricTableUnchanged(backfilled.metricTable,metricsBefore)
assert(isequal(backfilled.pointData,pointsBefore))
assert(isequal(backfilled.caseStatus,statusBefore))
assert(isequal(cellfun(@(car) car.M,backfilled.carCell(:,1)),massBefore))
assert(isequal(cellfun(@(car) car.ggPoints,backfilled.carCell(:,1), ...
    'UniformOutput',false),ggBefore))
assert(height(backfilled.rampBackfillDiagnostics) == 6)
assert(all(backfilled.rampBackfillDiagnostics.status == "failed"))
assert(all(backfilled.rampBackfillDiagnostics.error_identifier == ...
    "doeRunCase:unsolvedCar"))

function assertError(f,id)
try
    f();
    error('test:missingError','Expected %s',id)
catch ME
    assert(strcmp(ME.identifier,id),ME.message)
end
end

function assertMetricTableUnchanged(actual,expected)
assert(isequal(actual.Properties.VariableNames,expected.Properties.VariableNames))
changed = strings(0,1);
for j = 1:width(expected)
    name = expected.Properties.VariableNames{j};
    if isequaln(actual.(name),expected.(name)), continue, end
    rows = zeros(0,1);
    for i = 1:height(expected)
        if ~isequaln(actual.(name)(i,:),expected.(name)(i,:))
            rows(end+1,1) = i; %#ok<AGROW>
        end
    end
    changed(end+1,1) = string(name) + " rows " + mat2str(rows'); %#ok<AGROW>
end
assert(isempty(changed),"Ramp backfill changed metricTable: " + strjoin(changed,"; "))
assert(isequaln(actual,expected), ...
    'Ramp backfill changed metricTable values or table metadata.')
end
