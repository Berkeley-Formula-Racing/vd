%% test_doe_checkpoint
fullCarModels = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(fullCarModels,'setup_paths.m'))
d = tempname;
mkdir(d)
cleaner = onCleanup(@() rmdir(d,'s'));
path = fullfile(d,'DOE_checkpoint.mat');

R = struct('signature',"fixed-space",'mode',"sensitivity", ...
    'maxCases',16,'numWorkers',2,'batchSize',4, ...
    'objective',struct('penalty',1),'ramps',struct('enabled',false));
state = struct('resolvedStudy',R,'U',[0.2;0.8], ...
    'metricTable',table([1;2],'VariableNames',{'case_index'}), ...
    'batchNumber',1);

saveDOECheckpoint(path,state)
assert(isfile(path))
assert(~isfile(path + ".tmp.mat"))
[L,info] = loadDOECheckpoint(path,R);
assert(isequal(L.U,state.U) && L.batchNumber == 1)
assert(isempty(info.changed))

R2 = R;
R2.mode = "optimization";
R2.maxCases = 32;
R2.numWorkers = 8;
R2.batchSize = 6;
R2.objective = struct('penalty',2);
R2.ramps = struct('enabled',true);
[L2,info] = loadDOECheckpoint(path,R2);
expectedChanged = ["mode","maxCases","numWorkers","batchSize", ...
    "objective","ramps"];
assert(all(ismember(expectedChanged,info.changed)))
assert(isequal(L2.resolvedStudy.mode,R2.mode))
assert(isequal(L2.resolvedStudy.maxCases,R2.maxCases))
assert(isequal(L2.resolvedStudy.numWorkers,R2.numWorkers))
assert(isequal(L2.resolvedStudy.batchSize,R2.batchSize))
assert(isequal(L2.resolvedStudy.objective,R2.objective))
assert(isequal(L2.resolvedStudy.ramps,R2.ramps))

Rbad = R;
Rbad.signature = "different-space";
assertError(@() loadDOECheckpoint(path,Rbad),'loadDOECheckpoint:designMismatch')

function assertError(f,id)
try
    f();
    error('test:missingError','Expected %s',id)
catch ME
    assert(strcmp(ME.identifier,id),ME.message)
end
end
