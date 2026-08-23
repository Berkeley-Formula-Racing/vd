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

tic
serial = doeRunBatch(cars,eventParams,study,1);
serialElapsed = toc;
assert(serial{1}.status == "complete")
assert(isfinite(serial{1}.scoreBreakdown.objective_score))
assert(~isempty(serial{1}.rampSummary))
fprintf('DOE ramp batch serial elapsed: %.3f s\n',serialElapsed)

rampOnly = doeRunBatch({serial{1}.car,serial{1}.accelCar}, ...
    eventParams,study,1,"rampOnly");
assert(rampOnly{1}.status == "complete")
assert(isequal(rampOnly{1}.points,serial{1}.points))
assert(abs(rampOnly{1}.scoreBreakdown.objective_score - ...
    serial{1}.scoreBreakdown.objective_score) < 1e-8)

failed = doeRunCase([],[],eventParams,study,99);
assert(failed.status == "failed")
assert(strlength(failed.errorIdentifier) > 0)

if license('test','Distrib_Computing_Toolbox')
    study.numWorkers = min(2,feature('numcores'));
    tic
    parallel = doeRunBatch([cars;cars],eventParams,study,[1;2]);
    parallelElapsed = toc;
    assert(all(cellfun(@(x) x.status == "complete",parallel)))
    a = serial{1}.rampSummary;
    b = parallel{1}.rampSummary;
    assert(max(abs(a.K_linear-b.K_linear),[],'omitnan') < 1e-8)
    assert(abs(serial{1}.scoreBreakdown.objective_score - ...
        parallel{1}.scoreBreakdown.objective_score) < 1e-8)
    fprintf('DOE ramp batch parallel elapsed: %.3f s\n',parallelElapsed)
else
    fprintf('DOE ramp batch parallel elapsed: skipped (no Parallel Computing Toolbox)\n')
end
