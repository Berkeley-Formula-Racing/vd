%% test_doe_parallel_batch
testDir = fileparts(mfilename('fullpath'));
modelRoot = fileparts(testDir);
cd(modelRoot); setup_paths

[cars,eventParams] = carConfig();
study = DOEStudyConfig();
study.events = ["skidpad","accel","autocross","endurance"];
study.ramps.enabled = true;
study.ramps.speeds = 10;
study.ramps.nRamp = 4;
study.ramps.nBisect = 1;
study.ramps.saveFullPoints = false;
study.numWorkers = 0;
study.objective.penalties.energy.enabled = true;
study.objective.penalties.energy.threshold_kJ = 0;
study.objective.penalties.energy.pointsPerKJ = 0.01;

tic
serial = doeRunBatch(cars,eventParams,study,1);
serialElapsed = toc;
assert(serial{1}.status == "complete")
assert(isfinite(serial{1}.scoreBreakdown.objective_score))
assert(~isempty(serial{1}.rampSummary))
assert(serial{1}.scoreBreakdown.penalty_energy > 0)
publicPointFields = {'skidpad';'accel';'autocross';'endurance';'total'};
assert(isequal(fieldnames(serial{1}.car.comp.points),publicPointFields))
assert(all(structfun(@(x) isnumeric(x) && isscalar(x), ...
    serial{1}.car.comp.points)))
assert(isequal(serial{1}.car.comp.doeScoreBreakdown, ...
    serial{1}.scoreBreakdown))
fprintf('DOE ramp batch serial elapsed: %.3f s\n',serialElapsed)

rampStudy = study;
rampStudy.objective.penalties.energy.pointsPerKJ = 1;
rampOnly = doeRunBatch({serial{1}.car,serial{1}.accelCar}, ...
    eventParams,rampStudy,1,"rampOnly");
assert(rampOnly{1}.status == "complete")
assert(isequal(rampOnly{1}.points,serial{1}.points))
assert(isequal(rampOnly{1}.scoreBreakdown,serial{1}.scoreBreakdown))

mixed = doeRunBatch({serial{1}.car,serial{1}.accelCar;[] ,[]}, ...
    eventParams,rampStudy,[1;99],"rampOnly");
assert(mixed{1}.status == "complete")
assert(mixed{2}.status == "failed")
assert(strlength(mixed{2}.errorIdentifier) > 0)

subsetStudy = study;
subsetStudy.events = "skidpad";
subsetStudy.ramps.enabled = false;
subset = doeRunCase(cars{1,1},cars{1,2},eventParams,subsetStudy,3);
assert(subset.status == "complete")
assert(isempty(subset.car.comp.points))
assert(all(structfun(@isnan,subset.points)))
assert(isequaln(subset.car.comp.doeScoreBreakdown,subset.scoreBreakdown))

subsetStudy.ramps.enabled = true;
subsetBackfill = doeRunCase(subset.car,subset.accelCar,eventParams, ...
    subsetStudy,3,"rampOnly");
assert(subsetBackfill.status == "complete")
assert(isequaln(subsetBackfill.points,subset.points))
assert(isequaln(subsetBackfill.scoreBreakdown,subset.scoreBreakdown))

if license('test','Distrib_Computing_Toolbox')
    study.numWorkers = min(2,feature('numcores'));
    tic
    parallel = doeRunBatch([cars;cars],eventParams,study,[1;2]);
    parallelElapsed = toc;
    assert(all(cellfun(@(x) x.status == "complete",parallel)))
    a = serial{1}.rampSummary;
    b = parallel{1}.rampSummary;
    assertPointsEqual(serial{1}.points,parallel{1}.points)
    assertRampSummaryEqual(a,b)
    assert(abs(serial{1}.scoreBreakdown.objective_score - ...
        parallel{1}.scoreBreakdown.objective_score) < 1e-8)
    fprintf('DOE ramp batch parallel elapsed: %.3f s\n',parallelElapsed)
else
    fprintf('DOE ramp batch parallel elapsed: skipped (no Parallel Computing Toolbox)\n')
end

function assertPointsEqual(a,b)
names = {'skidpad','accel','autocross','endurance','total'};
for i = 1:numel(names)
    assert(abs(a.(names{i})-b.(names{i})) < 1e-8)
end
end

function assertRampSummaryEqual(a,b)
assert(isequal(a.Properties.VariableNames,b.Properties.VariableNames))
assert(height(a) == height(b))
isNumeric = varfun(@isnumeric,a,'OutputFormat','uniform');
names = a.Properties.VariableNames(isNumeric);
for i = 1:numel(names)
    x = a.(names{i});
    y = b.(names{i});
    finite = isfinite(x);
    assert(isequal(finite,isfinite(y)))
    assert(all(abs(x(finite)-y(finite)) < 1e-8,'all'))
end
end
