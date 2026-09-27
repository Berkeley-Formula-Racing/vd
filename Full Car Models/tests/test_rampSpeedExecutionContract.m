function tests = test_rampSpeedExecutionContract
tests = functiontests(localfunctions);
end

function testOneRowPerTaskAndNaNForMissingMetrics(testCase)
tasks = taskTable([5; 10; 15],[1; 2; 3],["requested";"requested";"refined"], ...
    [1; 1; 2],["";"";"gear_transition"]);
request = struct("rampType","longitudinal","setupKey","setup-1", ...
    "solverProfile","accurate");
[results,execution] = rampSpeed.executeSpeedPlan(request,tasks, ...
    @mockPointSolver,struct());
run = rampSpeed.assembleRun(tasks,results,request,execution);

verifyEqual(testCase,height(run.perSpeed),3);
verifyEqual(testCase,run.perSpeed.speedIndex,[1;2;3]);
verifyEqual(testCase,run.perSpeed.status, ...
    ["converged";"infeasible";"solver_failed"]);
verifyTrue(testCase,run.perSpeed.valid(1));
verifyFalse(testCase,any(run.perSpeed.valid(2:3)));
verifyEqual(testCase,run.perSpeed.aLong_mps2(1),5);
verifyTrue(testCase,isnan(run.perSpeed.aLong_mps2(2)));
verifyTrue(testCase,isnan(run.perSpeed.aLong_mps2(3)));
verifyEqual(testCase,numel(run.raw.diagnostics),3);
end

function testSolverFailedAndInfeasibleRemainDistinct(testCase)
tasks = taskTable([5; 10],[1; 2],["requested";"requested"],[1;1],["";""]);
request = struct("rampType","longitudinal","setupKey","setup-1");
[results,execution] = rampSpeed.executeSpeedPlan(request,tasks, ...
    @mockPointSolver,struct());
run = rampSpeed.assembleRun(tasks,results,request,execution);

verifyEqual(testCase,run.perSpeed.status(2),"infeasible");
verifyEqual(testCase,run.perSpeed.reason(2),"mock physical infeasibility");
verifyEqual(testCase,run.perSpeed.status(1),"converged");
end

function testCancellationFillsRemainingTasks(testCase)
tasks = taskTable([5; 10; 15],[1;2;3],repmat("requested",3,1), ...
    ones(3,1),repmat("",3,1));
request = struct("rampType","longitudinal","setupKey","setup-1");
cancelled = false;
control = struct("onProgress",@afterFirst, ...
    "shouldCancel",@isCancelled);
[results,execution] = rampSpeed.executeSpeedPlan(request,tasks, ...
    @mockPointSolver,control);
run = rampSpeed.assembleRun(tasks,results,request,execution);

verifyEqual(testCase,execution.status,"cancelled");
verifyEqual(testCase,run.perSpeed.status, ...
    ["converged";"cancelled";"cancelled"]);
verifyTrue(testCase,all(isnan(run.perSpeed.aLong_mps2(2:3))));

    function afterFirst(event)
        if event.completedTasks >= 1
            cancelled = true;
        end
    end

    function value = isCancelled()
        value = cancelled;
    end
end

function testDuplicateStableSpeedIdsAreRejected(testCase)
tasks = taskTable([5; 10],[1;1],["requested";"refined"],[1;2],["";"retry"]);
request = struct("rampType","longitudinal");
verifyError(testCase,@() rampSpeed.executeSpeedPlan(request,tasks, ...
    @mockPointSolver,struct()),"rampSpeed:duplicateSpeedIndex");
end

function testProgressEventsCarryStableTaskIdentity(testCase)
tasks = taskTable([5; 10],[8;9],repmat("requested",2,1), ...
    ones(2,1),repmat("",2,1));
request = struct("rampType","longitudinal");
events = struct("speedIndex",{}, "speed_mps",{}, ...
    "completedTasks",{}, "plannedTasks",{}, "status",{});
control = struct("onProgress",@capture);
rampSpeed.executeSpeedPlan(request,tasks,@mockPointSolver,control);
verifyEqual(testCase,[events.speedIndex],[8 9]);
verifyEqual(testCase,[events.completedTasks],[1 2]);

    function capture(event)
        events(end+1) = event;
    end
end

function result = mockPointSolver(task,~,~,~)
status = "converged";
metrics = struct("aLong_mps2",task.speed_mps);
diagnostics = struct("reason","");
if task.speed_mps == 10
    status = "infeasible";
    metrics = struct();
    diagnostics.reason = "mock physical infeasibility";
elseif task.speed_mps == 15
    status = "solver_failed";
    metrics = struct();
    diagnostics.reason = "mock optimizer failure";
end
result = rampSpeed.makeSpeedResult(task,status,metrics,diagnostics);
end

function tasks = taskTable(speeds,indices,origins,passes,reasons)
tasks = table(indices,double(speeds),origins,passes,reasons, ...
    repmat("planned",numel(speeds),1), ...
    'VariableNames',{'speedIndex','speed_mps','origin','passIndex', ...
    'refinementReason','status'});
end
