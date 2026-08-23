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
