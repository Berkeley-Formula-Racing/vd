# Adaptive DOE and Parallel Ramp Metrics Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a resumable adaptive vehicle DOE that supports sensitivity, points optimization, and hybrid sampling while running optional ramp metrics across the outer parallel case pool.

**Architecture:** `DOEStudyConfig` owns study ranges and settings, while `carConfig` remains the scalar baseline and explicit car builder. `runAdaptiveDOE` executes checkpointed batches; focused utilities resolve parameters, fit GPs, select candidates, score cases, and persist state. Per-case event, metric, and ramp work stays inside one outer `parfor`, and final analysis consumes cached metrics instead of reprocessing all saved cars.

**Tech Stack:** MATLAB, Optimization Toolbox, Statistics and Machine Learning Toolbox (`fitrgp`, `fitcensemble`), Parallel Computing Toolbox, existing QSS `gg2`/`makeGG`/`Events2`/`rampSweep` model chain.

**Spec:** `docs/superpowers/specs/2026-08-23-adaptive-doe-design.md`

## Global Constraints

- Preserve all unrelated and pre-existing working-tree edits, especially current changes in `carConfig.m`, `SteadyStateLapsim.m`, `DOE_Fitting.m`, and DOE analysis utilities.
- `carConfig()` must continue to construct the current nominal `carCell{1,1}` and acceleration car with no external study configuration.
- DOE ranges live only in `DOEStudyConfig.m`; no range arrays remain embedded in baseline parameter assignments.
- Percentage ranges are valid only for positive, nonzero baselines; zero, negative, or sign-crossing parameters use absolute ranges.
- The only active parallel loop is the outer loop over car cases.
- Mode, penalties, worker count, ramp storage, batch size, and maximum cases can change on resume; baseline and resolved design bounds cannot.
- Optimization maximizes existing modeled skidpad, acceleration, autocross, and endurance points minus explicit penalties. Efficiency points remain excluded and are never labeled as modeled.
- Persist every completed batch before fitting or selecting the next batch.
- Follow TDD: every behavioral change starts with a failing test and ends with a focused commit.

## File map

- Create `Full Car Models/DOEStudyConfig.m`: user-editable study definition.
- Create `Full Car Models/runAdaptiveDOE.m`: adaptive orchestration and final result export.
- Modify `Full Car Models/SteadyStateLapsim.m`: thin entry-point wrapper.
- Modify `Full Car Models/carConfig.m`: scalar baseline plus explicit-table construction and deprecated one-shot sampling compatibility.
- Modify `Full Car Models/utilities/parameters_loop.m`: explicit design-table sampling mode.
- Create `Full Car Models/utilities/doeResolveStudy.m`: study validation, absolute bounds, normalized transforms, and signature.
- Create `Full Car Models/utilities/doeInitialDesign.m`: repeatable LHS plus baseline point.
- Create `Full Car Models/utilities/doeFitSurrogates.m`: response GPs, objective GP, and feasibility model.
- Create `Full Car Models/utilities/doeSelectAdaptiveBatch.m`: sensitivity, optimization, and hybrid acquisition with diversity filtering.
- Create `Full Car Models/utilities/doeExpectedImprovement.m`: maximization expected improvement.
- Create `Full Car Models/utilities/doeScoreCase.m`: raw points and penalty accounting.
- Create `Full Car Models/utilities/doeCaseMetrics.m`: one-car metric extraction.
- Modify `Full Car Models/utilities/doeMetrics.m`: compatibility wrapper around the per-case reducer.
- Create `Full Car Models/utilities/doeRunCase.m`: one complete QSS/event/ramp case.
- Create `Full Car Models/utilities/doeRunBatch.m`: serial or outer-`parfor` batch execution.
- Create `Full Car Models/utilities/saveDOECheckpoint.m`: atomic checkpoint writer.
- Create `Full Car Models/utilities/loadDOECheckpoint.m`: resume validation.
- Modify `Full Car Models/utilities/doeAnalyze.m`: cached metrics, GP validation, and preferred-model metadata.
- Create `Full Car Models/utilities/predictDOEModel.m`: common prediction interface for linear or GP models.
- Create `Full Car Models/utilities/doeSobol.m`: surrogate-based first-order and total-order indices.
- Modify `Full Car Models/utilities/doePlotCatalog.m`, `plotDOEMetrics.m`, and `DOE_Fitting.m`: model comparison and Sobol plots.
- Add focused scripts under `Full Car Models/tests/` for every task below.

---

### Task 1: Study configuration and range resolution

**Files:**
- Create: `Full Car Models/DOEStudyConfig.m`
- Create: `Full Car Models/utilities/doeResolveStudy.m`
- Create: `Full Car Models/tests/test_doe_study_config.m`

**Interfaces:**
- Consumes: one-row baseline table whose variable names match `parameters_loop` names.
- Produces: `study = DOEStudyConfig()` and `resolved = doeResolveStudy(study,baselineTable)`.
- `resolved.parameters` contains `name`, `rangeType`, `lowerInput`, `upperInput`, `baseline`, `lowerPhysical`, and `upperPhysical`.
- `resolved.toPhysical(U)` maps normalized rows to a table containing only varied physical parameters.
- `resolved.signature` is a deterministic JSON string over names, baseline values, range types, and resolved bounds.

- [ ] **Step 1: Write the failing range-resolution test**

```matlab
%% test_doe_study_config
setup_paths
baseline = table(162,0.34,-1,0,'VariableNames', ...
    {'mass','R_sf','gamma_f','static_r_toe'});
study = struct();
study.parameters = table( ...
    ["mass";"R_sf";"gamma_f";"static_r_toe"], ...
    [-5;-10;-2;-0.2],[5;10;0;0.2], ...
    ["percent";"percent";"absolute";"absolute"], ...
    'VariableNames',{'name','lower','upper','rangeType'});
R = doeResolveStudy(study,baseline);
assert(abs(R.parameters.lowerPhysical(1)-153.9) < 1e-12)
assert(abs(R.parameters.upperPhysical(2)-0.374) < 1e-12)
T = R.toPhysical([0 0 0 0; 1 1 1 1]);
assert(isequal(T.Properties.VariableNames,cellstr(study.parameters.name)'))
assert(abs(T.gamma_f(1)+2) < 1e-12 && abs(T.gamma_f(2)) < 1e-12)
assert(R.signature == doeResolveStudy(study,baseline).signature)

bad = study; bad.parameters.rangeType(3) = "percent";
assertError(@() doeResolveStudy(bad,baseline),'doeResolveStudy:badPercentBaseline')
bad = study; bad.parameters.name(1) = "not_a_car_parameter";
assertError(@() doeResolveStudy(bad,baseline),'doeResolveStudy:unknownParameter')

function assertError(f,id)
try, f(); error('test:missingError','Expected %s',id)
catch ME, assert(strcmp(ME.identifier,id),ME.message)
end
end
```

- [ ] **Step 2: Run the test and verify the missing API failure**

Run:

```powershell
matlab -batch "cd('Full Car Models'); run('tests/test_doe_study_config.m')"
```

Expected: failure because `doeResolveStudy` is undefined.

- [ ] **Step 3: Implement the user-editable study defaults**

Define all settings explicitly in `DOEStudyConfig.m`, including:

```matlab
function study = DOEStudyConfig()
study.name = "vehicle_design_2026";
study.mode = "sensitivity";
study.randomSeed = 1;
study.initialCases = 128;
study.batchSize = 64;
study.maxCases = 512;
study.numWorkers = 16;
study.allowSerialFallback = false;
study.resume = true;
study.parameters = table( ...
    ["mass";"wheelbase";"weight_dist";"track_width";"cg_height"; ...
     "R_sf";"accel_cda";"accel_cla";"final_drive"; ...
     "gamma_f";"gamma_r";"p_i"], ...
    [-5;-5;-5;-5;-10;-15;-15;-15;-8;-2;-2;-10], ...
    [ 5; 5; 5; 5; 10; 15; 15; 15; 8; 0; 0; 10], ...
    [repmat("percent",9,1);"absolute";"absolute";"percent"], ...
    'VariableNames',{'name','lower','upper','rangeType'});
study.adaptive = struct('responses', ...
    {{'modeled_dynamic_points','t_autox','t_accel','t_skid', ...
      'total_work_kJ','understeer_gradient_10_deg_per_g', ...
      'understeer_gradient_25_deg_per_g'}}, ...
    'hybridSensitivityFraction',0.70,'candidatePoolSize',20000, ...
    'minimumDistance',0.03,'kernel','ardmatern32');
study.events = ["skidpad","accel","autocross","endurance"];
study.ramps = struct('enabled',true,'speeds',[10 25],'nRamp',8, ...
    'nBisect',3,'mode','balanced','saveFullPoints',false);
study.output = struct('directory',fullfile(fileparts(which('DOEStudyConfig')), ...
    'DOE_output'),'checkpoint','DOE_checkpoint.mat', ...
    'results','DOE_results.mat','analysis','DOE_analysis.mat');
study.objective = defaultObjective();
end
```

Use this explicit penalty schema in local `defaultObjective`:

```matlab
function objective = defaultObjective()
objective.primary = "modeled_dynamic_points";
objective.penalties.invalidCase = struct( ...
    'enabled',true,'fixedPoints',575);
objective.penalties.solveFailure = struct( ...
    'enabled',false,'thresholdFraction',0,'pointsPerFraction',0);
objective.penalties.wheelLift = struct( ...
    'enabled',false,'threshold_N',0,'pointsPerN',0);
objective.penalties.energy = struct( ...
    'enabled',false,'threshold_kJ',Inf,'pointsPerKJ',0);
objective.penalties.understeer = struct( ...
    'enabled',false,'lower_deg_per_g',-Inf,'upper_deg_per_g',Inf, ...
    'pointsPerDegPerG',0,'metric',"understeer_gradient_25_deg_per_g");
objective.penalties.rebalance = struct( ...
    'enabled',false,'target_mps',NaN,'pointsPerMps',0);
end
```

For `optimization` and `hybrid` modes, validate that `study.events` contains `skidpad`, `accel`, `autocross`, and `endurance`. Sensitivity mode may omit events, but point fields are then NaN and cannot appear in `study.adaptive.responses`.

- [ ] **Step 4: Implement validation and normalized mapping**

In `doeResolveStudy`, validate scalar settings and parameter rows, calculate physical bounds, create function handles for normalized-to-physical mapping, and form the signature with:

```matlab
payload = struct('name',{cellstr(P.name)},'baseline',P.baseline, ...
    'rangeType',{cellstr(P.rangeType)}, ...
    'lower',P.lowerPhysical,'upper',P.upperPhysical);
resolved.signature = string(jsonencode(payload));
```

- [ ] **Step 5: Run the focused test**

Run the command from Step 2. Expected: script exits successfully with no assertion failures.

- [ ] **Step 6: Commit the configuration boundary**

```powershell
git add -- "Full Car Models/DOEStudyConfig.m" "Full Car Models/utilities/doeResolveStudy.m" "Full Car Models/tests/test_doe_study_config.m"
git commit -m "feat: define external DOE study configuration"
```

---

### Task 2: Scalar baseline and explicit car-case construction

**Files:**
- Modify: `Full Car Models/carConfig.m`
- Modify: `Full Car Models/utilities/parameters_loop.m`
- Create: `Full Car Models/tests/test_doe_explicit_cars.m`

**Interfaces:**
- Produces: `[carCell,eventParams,designTable,baselineTable] = carConfig()` with one baseline row.
- Produces: `[carCell,eventParams,designTable] = carConfig("Explicit",overrideTable)`.
- `parameters_loop(...,"Explicit",overrideTable)` overlays named columns onto one-row scalar baseline values and returns a complete design table.

- [ ] **Step 1: Capture baseline and explicit-override expectations in a failing test**

```matlab
%% test_doe_explicit_cars
setup_paths
[base,eventParams,X0,B] = carConfig();
assert(size(base,1)==1 && size(base,2)==2)
assert(height(X0)==1 && height(B)==1)
assert(abs(base{1,1}.M-226) < 1e-10)
assert(abs(base{1,1}.W_b-62*0.0254) < 1e-12)
assert(abs(base{1,1}.aero.C_lA-3.969) < 1e-12)
assert(isfield(eventParams,'winning_time'))

Q = table([155;170],[0.30;0.50],[3.5;4.2], ...
    'VariableNames',{'mass','R_sf','cla'});
[cars,~,X] = carConfig("Explicit",Q);
assert(size(cars,1)==2)
assert(abs(cars{1,1}.M-(155+64)) < 1e-10)
assert(abs(cars{2,1}.R_sf-0.50) < 1e-12)
assert(abs(cars{1,1}.aero.C_lA-3.5) < 1e-12)
assert(all(abs(X.p_i-12) < 1e-12))
```

- [ ] **Step 2: Run the test and verify it fails on the current multi-row baseline or missing explicit mode**

```powershell
matlab -batch "cd('Full Car Models'); run('tests/test_doe_explicit_cars.m')"
```

- [ ] **Step 3: Make baseline assignments scalar and remove the embedded DOE range block**

Keep nominal values, including the first nominal element of the current mass, wheelbase, weight distribution, and track-width vectors. Do not change calibration numbers, tire scaling, inertia, rolling resistance, event parameters, or loaded tire fits.

- [ ] **Step 4: Add explicit sampling to `parameters_loop`**

Extend `sampleValues` so `samplingType="Explicit"` starts from the scalar baseline vector and replaces only provided columns:

```matlab
P = repmat(cellfun(@(x) x(1),values),height(sampleInput),1);
for k = 1:width(sampleInput)
    j = find(strcmp(names,sampleInput.Properties.VariableNames{k}),1);
    if isempty(j)
        error('parameters_loop:unknownExplicitParameter', ...
            'Unknown explicit parameter %s.',sampleInput.Properties.VariableNames{k});
    end
    P(:,j) = sampleInput{:,k};
end
```

Validate that every baseline value is scalar before explicit construction.

- [ ] **Step 5: Add the `carConfig("Explicit",table)` dispatcher and baseline-table output**

Retain no-argument compatibility. Implement deprecated `LHS` and `Random` calls by loading `DOEStudyConfig`, resolving its ranges against `baselineTable`, creating one non-adaptive design, and warning with identifier `carConfig:legacyDOE`.

- [ ] **Step 6: Run baseline, explicit-car, and existing DOE metric tests**

```powershell
matlab -batch "cd('Full Car Models'); run('tests/test_doe_explicit_cars.m'); run('tests/test_doe_metrics.m')"
```

Expected: both scripts pass; the no-argument car has one row.

- [ ] **Step 7: Commit explicit construction**

```powershell
git add -- "Full Car Models/carConfig.m" "Full Car Models/utilities/parameters_loop.m" "Full Car Models/tests/test_doe_explicit_cars.m"
git commit -m "refactor: separate baseline car from DOE ranges"
```

---

### Task 3: Per-case metrics and auditable points penalties

**Files:**
- Create: `Full Car Models/utilities/doeCaseMetrics.m`
- Create: `Full Car Models/utilities/doeScoreCase.m`
- Modify: `Full Car Models/utilities/doeMetrics.m`
- Create: `Full Car Models/tests/test_doe_case_scoring.m`

**Interfaces:**
- Produces: `[row,details] = doeCaseMetrics(car,caseIndex,rampResult)` where `row` is a one-row table.
- Produces: `[score,breakdown] = doeScoreCase(points,row,objectiveConfig)`.
- `breakdown` contains raw event points, `modeled_dynamic_points`, one column per penalty, `total_penalty_points`, and `objective_score`.

- [ ] **Step 1: Write a failing synthetic scoring test**

```matlab
%% test_doe_case_scoring
setup_paths
points = struct('skidpad',60,'accel',80,'autocross',100, ...
    'endurance',220,'total',460);
M = table(true,-12,310,2.5,18, ...
    'VariableNames',{'valid','min_Fz_N','total_work_kJ', ...
    'understeer_gradient_25_deg_per_g','rebalance_speed_mps'});
study = DOEStudyConfig();
cfg = study.objective;
cfg.penalties.wheelLift.enabled = true;
cfg.penalties.wheelLift.threshold_N = 0;
cfg.penalties.wheelLift.pointsPerN = 0.25;
cfg.penalties.energy.enabled = true;
cfg.penalties.energy.threshold_kJ = 300;
cfg.penalties.energy.pointsPerKJ = 0.1;
[score,B] = doeScoreCase(points,M,cfg);
assert(abs(B.modeled_dynamic_points-460) < 1e-12)
assert(abs(B.penalty_wheel_lift-3) < 1e-12)
assert(abs(B.penalty_energy-1) < 1e-12)
assert(abs(score-456) < 1e-12)
assert(abs(B.total_penalty_points-4) < 1e-12)

M.valid = false;
[score,B] = doeScoreCase(points,M,cfg);
assert(B.penalty_invalid_case == cfg.penalties.invalidCase.fixedPoints)
assert(score == B.modeled_dynamic_points-B.total_penalty_points)
```

- [ ] **Step 2: Run the test and verify `doeScoreCase` is missing**

```powershell
matlab -batch "cd('Full Car Models'); run('tests/test_doe_case_scoring.m')"
```

- [ ] **Step 3: Extract one-car logic from `doeMetrics`**

Move the body of the current car loop into `doeCaseMetrics`. Preserve every existing metric name. When `rampResult` is empty, return NaN for ramp-specific metrics; when present, extract understeer gradients and rebalance speed from `rampResult.perSpeed`.

- [ ] **Step 4: Implement deterministic penalty accounting**

Implement fixed invalid-case penalty and linear threshold penalties using explicit formulas such as:

```matlab
penaltyWheelLift = enabled * pointsPerN * ...
    max(threshold_N-row.min_Fz_N,0);
penaltyEnergy = enabled * pointsPerKJ * ...
    max(row.total_work_kJ-threshold_kJ,0);
score = modeledPoints - sumPenalty;
```

Understeer bounds use distance outside `[lower,upper]`; rebalance uses absolute error from its target.

- [ ] **Step 5: Rewrite `doeMetrics` as the compatibility all-cars wrapper**

Keep its current public signature. It calls `doeCaseMetrics` for each car and retains optional serial ramp behavior for older callers. Do not add parallelism here; Task 6 moves ramps into the outer batch loop.

- [ ] **Step 6: Run scoring and existing metric tests**

```powershell
matlab -batch "cd('Full Car Models'); run('tests/test_doe_case_scoring.m'); run('tests/test_doe_metrics.m')"
```

- [ ] **Step 7: Commit metric and scoring boundaries**

```powershell
git add -- "Full Car Models/utilities/doeCaseMetrics.m" "Full Car Models/utilities/doeScoreCase.m" "Full Car Models/utilities/doeMetrics.m" "Full Car Models/tests/test_doe_case_scoring.m"
git commit -m "feat: add auditable DOE points objective"
```

---

### Task 4: Initial design, Gaussian-process fitting, and acquisition

**Files:**
- Create: `Full Car Models/utilities/doeInitialDesign.m`
- Create: `Full Car Models/utilities/doeFitSurrogates.m`
- Create: `Full Car Models/utilities/doeExpectedImprovement.m`
- Create: `Full Car Models/utilities/doeSelectAdaptiveBatch.m`
- Create: `Full Car Models/tests/test_doe_adaptive_sampling.m`

**Interfaces:**
- Produces: `[U,T] = doeInitialDesign(resolved,n,seed)`.
- Produces: `models = doeFitSurrogates(U,metricTable,responseNames)` with `models.responses`, `models.objective`, and `models.feasibility`.
- Produces: `[Unew,selection] = doeSelectAdaptiveBatch(state,study,candidateU)`; optional `candidateU` makes tests deterministic.
- `selection` records mode, acquisition value, source (`sensitivity` or `optimization`), and nearest existing-point distance.

- [ ] **Step 1: Write failing deterministic sampling and acquisition tests**

```matlab
%% test_doe_adaptive_sampling
setup_paths
baseline = table(100,2,'VariableNames',{'mass','cla'});
study = DOEStudyConfig();
study.parameters = table(["mass";"cla"],[-10;-20],[10;20], ...
    ["percent";"percent"], ...
    'VariableNames',{'name','lower','upper','rangeType'});
study.adaptive.responses = {'response_a','response_b'};
study.batchSize = 4;
R = doeResolveStudy(study,baseline);
[U1,T1] = doeInitialDesign(R,12,7);
[U2,T2] = doeInitialDesign(R,12,7);
assert(isequal(U1,U2) && isequal(T1,T2))
assert(any(vecnorm(U1-R.baselineNormalized,2,2) < 1e-12))

M = table(true(12,1),U1(:,1)+U1(:,2),U1(:,1)-U1(:,2), ...
    500-(U1(:,1)-0.8).^2*100, ...
    'VariableNames',{'valid','response_a','response_b','objective_score'});
state = struct('U',U1,'metricTable',M);
C = [0.05 0.95;0.95 0.05;0.80 0.50;0.50 0.80;0.2 0.2;0.7 0.7;0.9 0.9];

study.mode = "sensitivity";
[Us,Ss] = doeSelectAdaptiveBatch(state,study,C);
assert(size(Us,1)==4 && all(Ss.source=="sensitivity"))
study.mode = "optimization";
[Uo,So] = doeSelectAdaptiveBatch(state,study,C);
assert(size(Uo,1)==4 && all(So.source=="optimization"))
study.mode = "hybrid"; study.adaptive.hybridSensitivityFraction=0.5;
[Uh,Sh] = doeSelectAdaptiveBatch(state,study,C);
assert(sum(Sh.source=="sensitivity")==2)
assert(sum(Sh.source=="optimization")==2)
assert(size(unique(Uh,'rows'),1)==4)
```

- [ ] **Step 2: Run the test and verify missing-function failure**

```powershell
matlab -batch "cd('Full Car Models'); run('tests/test_doe_adaptive_sampling.m')"
```

- [ ] **Step 3: Implement repeatable LHS and physical mapping**

Use the repository's base-MATLAB stratified LHS pattern, add `resolved.baselineNormalized`, reject duplicates, and map through `resolved.toPhysical`.

When `candidateU` is omitted, generate `study.adaptive.candidatePoolSize` normalized candidates from the saved random stream. Exclude completed rows before acquisition scoring and store the updated stream state in the checkpoint so resume produces the same future batches.

- [ ] **Step 4: Implement exact ARD GPs and feasibility fitting**

For each finite response:

```matlab
mdl = fitrgp(U(ok,:),y(ok),'KernelFunction','ardmatern32', ...
    'Standardize',true,'FitMethod','exact','PredictMethod','exact');
```

Fit `fitcensemble(U,valid)` only when both validity classes exist; otherwise store the scalar feasibility probability.

Skip an adaptive response whose finite observed spread is zero and record it in `models.skippedResponses`; fail only if every configured adaptive response is unavailable.

- [ ] **Step 5: Implement maximization expected improvement**

Use `z=(mu-best)./max(sigma,eps)` and:

```matlab
ei = (mu-best).*normcdf(z) + sigma.*normpdf(z);
ei(sigma <= eps) = max(mu(sigma <= eps)-best,0);
```

- [ ] **Step 6: Implement the three acquisitions and greedy diversity exclusion**

Sensitivity score is RMS predictive standard deviation after division by observed response spread. Optimization score is expected improvement times feasibility. Hybrid rounds `batchSize*hybridSensitivityFraction` to the sensitivity count and fills the remainder from optimization. Greedily remove candidates within normalized RMS distance `minimumDistance`, halving that threshold only when required to fill the batch.

- [ ] **Step 7: Run the adaptive-sampling test twice**

Run the command from Step 2 twice. Expected: both runs pass and select identical normalized rows.

- [ ] **Step 8: Commit adaptive selection**

```powershell
git add -- "Full Car Models/utilities/doeInitialDesign.m" "Full Car Models/utilities/doeFitSurrogates.m" "Full Car Models/utilities/doeExpectedImprovement.m" "Full Car Models/utilities/doeSelectAdaptiveBatch.m" "Full Car Models/tests/test_doe_adaptive_sampling.m"
git commit -m "feat: add adaptive DOE acquisition modes"
```

---

### Task 5: Atomic checkpoint and resume validation

**Files:**
- Create: `Full Car Models/utilities/saveDOECheckpoint.m`
- Create: `Full Car Models/utilities/loadDOECheckpoint.m`
- Create: `Full Car Models/tests/test_doe_checkpoint.m`

**Interfaces:**
- Produces: `saveDOECheckpoint(path,state)` using MAT v7.3 temporary replacement.
- Produces: `[state,resumeInfo] = loadDOECheckpoint(path,resolvedStudy)`.
- `resumeInfo.changed` lists allowed settings changed since the saved run.

- [ ] **Step 1: Write the failing checkpoint round-trip test**

```matlab
%% test_doe_checkpoint
setup_paths
d = tempname; mkdir(d); cleaner=onCleanup(@() rmdir(d,'s'));
path = fullfile(d,'DOE_checkpoint.mat');
R = struct('signature',"fixed-space",'mode',"sensitivity", ...
    'maxCases',16,'numWorkers',2,'batchSize',4);
state = struct('resolvedStudy',R,'U',[0.2;0.8], ...
    'metricTable',table([1;2],'VariableNames',{'case_index'}), ...
    'batchNumber',1);
saveDOECheckpoint(path,state)
[L,info] = loadDOECheckpoint(path,R);
assert(isequal(L.U,state.U) && L.batchNumber==1)
assert(isempty(info.changed))

R2=R; R2.mode="optimization"; R2.maxCases=32; R2.numWorkers=8;
[~,info] = loadDOECheckpoint(path,R2);
assert(all(ismember(["mode","maxCases","numWorkers"],info.changed)))
Rbad=R; Rbad.signature="different-space";
assertError(@() loadDOECheckpoint(path,Rbad),'loadDOECheckpoint:designMismatch')

function assertError(f,id)
try, f(); error('test:missingError','Expected %s',id)
catch ME, assert(strcmp(ME.identifier,id),ME.message)
end
end
```

- [ ] **Step 2: Run the test and verify missing checkpoint functions**

```powershell
matlab -batch "cd('Full Car Models'); run('tests/test_doe_checkpoint.m')"
```

- [ ] **Step 3: Implement atomic save**

Save a variable named `checkpoint` to `path + ".tmp.mat"` with `-v7.3`, verify the temporary file can be loaded, then replace the target with `movefile(tmp,path,'f')`. Delete the temporary file through an `onCleanup` only if it still exists.

- [ ] **Step 4: Implement strict design-space validation and allowed setting changes**

Require exact signature equality. Report differences in mode, maximum cases, workers, batch size, objective, and ramp settings without rejecting them. Return the saved state with those allowed fields updated from the requested study.

- [ ] **Step 5: Run the checkpoint test**

Run the command from Step 2. Expected: pass with the temporary directory removed by cleanup.

- [ ] **Step 6: Commit persistence**

```powershell
git add -- "Full Car Models/utilities/saveDOECheckpoint.m" "Full Car Models/utilities/loadDOECheckpoint.m" "Full Car Models/tests/test_doe_checkpoint.m"
git commit -m "feat: checkpoint and resume adaptive DOE"
```

---

### Task 6: One-case execution and parallel ramp batches

**Files:**
- Create: `Full Car Models/utilities/doeRunCase.m`
- Create: `Full Car Models/utilities/doeRunBatch.m`
- Create: `Full Car Models/tests/test_doe_parallel_batch.m`

**Interfaces:**
- Consumes: `result = doeRunCase(car,accelCar,eventParams,study,caseIndex,runMode)` where `runMode` defaults to `"full"` and also accepts `"rampOnly"` for already-solved cars.
- Produces: result fields `car`, `accelCar`, `metricRow`, `rampSummary`, `rampPoints`, `points`, `scoreBreakdown`, `status`, `errorIdentifier`, `errorMessage`, and `elapsed`.
- Produces: `results = doeRunBatch(carCell,eventParams,study,caseIndices,runMode)` with one result cell per case; `runMode` defaults to `"full"`.

- [ ] **Step 1: Write a failing serial-versus-parallel ramp equivalence test**

```matlab
%% test_doe_parallel_batch
setup_paths
[cars,eventParams] = carConfig();
study = DOEStudyConfig();
study.events = ["skidpad","accel","autocross","endurance"];
study.ramps.enabled = true;
study.ramps.speeds = 10;
study.ramps.nRamp = 4;
study.ramps.nBisect = 1;
study.ramps.saveFullPoints = false;
study.numWorkers = 0;
serial = doeRunBatch(cars,eventParams,study,1);
assert(serial{1}.status=="complete")
assert(isfinite(serial{1}.scoreBreakdown.objective_score))
assert(~isempty(serial{1}.rampSummary))

study.numWorkers = min(2,feature('numcores'));
parallel = doeRunBatch([cars;cars],eventParams,study,[1;2]);
assert(all(cellfun(@(x) x.status=="complete",parallel)))
a = serial{1}.rampSummary;
b = parallel{1}.rampSummary;
assert(max(abs(a.K_linear-b.K_linear),[],'omitnan') < 1e-8)
assert(abs(serial{1}.scoreBreakdown.objective_score- ...
    parallel{1}.scoreBreakdown.objective_score) < 1e-8)
```

- [ ] **Step 2: Run the test and verify missing batch APIs**

```powershell
matlab -batch "cd('Full Car Models'); run('tests/test_doe_parallel_batch.m')"
```

- [ ] **Step 3: Implement `doeRunCase` using the existing model chain**

The success path is exactly:

```matlab
rawGG = gg2(car,0);
car = makeGG(rawGG,car);
comp = Events2(car,accelCar,eventParams);
if any(study.events=="skidpad"), comp.Skidpad(); end
if any(study.events=="accel"), comp.Accel(); end
if any(study.events=="autocross"), comp.Autocross(); end
if any(study.events=="endurance"), comp.Endurance(); end
if all(ismember(["skidpad","accel","autocross","endurance"],study.events))
    points = comp.computePoints();
else
    points = struct('skidpad',NaN,'accel',NaN,'autocross',NaN, ...
        'endurance',NaN,'total',NaN);
end
car.comp = comp;
```

Run `rampSweep` afterward on the same worker when enabled, reduce metrics with `doeCaseMetrics`, and score with `doeScoreCase`. Catch case-level exceptions and return diagnostics without throwing out of the batch.

For `runMode="rampOnly"`, require a car that already contains solved g-g data and `car.comp`; skip `gg2`, event execution, and point recomputation, run only `rampSweep` and ramp metric reduction, and preserve the previously cached point/score fields on append.

- [ ] **Step 4: Implement the only outer parallel loop**

Use:

```matlab
results = cell(n,1);
parfor (j=1:n,study.numWorkers)
    results{j} = doeRunCase(carCell{j,1},carCell{j,2}, ...
        eventParams,study,caseIndices(j),runMode);
end
```

Set `numWorkers=0` for deterministic serial execution in tests and debugging. Do not start a nested pool or write shared files from workers.

- [ ] **Step 5: Run the equivalence test and record elapsed serial/parallel ramp time**

Run the command from Step 2. Expected: both modes produce matching point totals and ramp summaries. Record timings in the test output without asserting a speedup on a two-case smoke test.

- [ ] **Step 6: Commit parallel case execution**

```powershell
git add -- "Full Car Models/utilities/doeRunCase.m" "Full Car Models/utilities/doeRunBatch.m" "Full Car Models/tests/test_doe_parallel_batch.m"
git commit -m "feat: run DOE ramps across case workers"
```

---

### Task 7: Adaptive orchestration, batch append, and entry point

**Files:**
- Create: `Full Car Models/runAdaptiveDOE.m`
- Modify: `Full Car Models/SteadyStateLapsim.m`
- Create: `Full Car Models/tests/test_adaptive_doe_resume.m`

**Interfaces:**
- Produces: `state = runAdaptiveDOE(study)`.
- `state` contains `resolvedStudy`, `U`, `designTable`, `carCell`, `metricTable`, `rampData`, `caseStatus`, `selectionHistory`, `batchNumber`, `elapsed`, and `randomState`.
- Final `DOE_results.mat` contains existing `carCell`, `designTable`, and `eventParams`, plus cached `metricTable`, `rampData`, `study`, and `selectionHistory`.

- [ ] **Step 1: Write a failing small resume test with one varying parameter**

```matlab
%% test_adaptive_doe_resume
setup_paths
d=tempname; mkdir(d); cleaner=onCleanup(@() rmdir(d,'s'));
study=DOEStudyConfig();
study.name="resume_smoke";
study.parameters=table("mass",-2,2,"percent", ...
    'VariableNames',{'name','lower','upper','rangeType'});
study.initialCases=4; study.batchSize=2; study.maxCases=4;
study.numWorkers=0; study.mode="sensitivity";
study.events=["skidpad","accel","autocross","endurance"];
study.ramps.enabled=false;
study.output.directory=d;
state1=runAdaptiveDOE(study);
assert(height(state1.designTable)==4)

study.mode="optimization"; study.maxCases=6; study.resume=true;
state2=runAdaptiveDOE(study);
assert(height(state2.designTable)==6)
assert(isequal(state2.U(1:4,:),state1.U))
assert(all(state2.selectionHistory.source(end-1:end)=="optimization"))
assert(isfile(fullfile(d,study.output.results)))
```

- [ ] **Step 2: Run the test and verify `runAdaptiveDOE` is missing**

```powershell
matlab -batch "cd('Full Car Models'); run('tests/test_adaptive_doe_resume.m')"
```

- [ ] **Step 3: Implement new-study initialization**

Check `fitrgp` before any simulation. If Parallel Computing Toolbox is unavailable and `numWorkers>0`, throw `runAdaptiveDOE:noParallelToolbox` unless `allowSerialFallback=true`, in which case set the effective workers to zero and record the fallback. Then call `carConfig()` for the baseline table and event parameters, resolve the study, create the initial design, initialize empty typed state fields, and save a checkpoint before executing the first batch.

- [ ] **Step 4: Implement the batch loop**

For each batch:

1. map normalized points to an explicit table;
2. build cars with `carConfig("Explicit",batchTable)`;
3. call `doeRunBatch`;
4. append cars, metrics, ramps, statuses, and acquisition metadata;
5. recompute penalties from raw cached metrics;
6. atomically save the checkpoint; and
7. fit/select the next batch only after that save succeeds.

If valid cases are not greater than the predictor count, select another space-filling batch instead of fitting GPs.

- [ ] **Step 5: Implement resume and ramp backfill**

Load and validate the checkpoint. Never rebuild completed cases. If ramps have become enabled and old cases have no ramp summaries, call `doeRunBatch(...,"rampOnly")` for those cars, append the results, and checkpoint before continuing adaptive selection.

- [ ] **Step 6: Export compatibility results and simulation log**

Write `DOE_results.mat` with `-v7.3`, call `simLog.start/finish`, and include mode, events, workers, batches, valid/invalid counts, output paths, and elapsed time.

- [ ] **Step 7: Replace `SteadyStateLapsim` settings with the thin wrapper**

```matlab
%% SteadyStateLapsim - adaptive DOE runner
clear classes
setup_paths
study = DOEStudyConfig();
state = runAdaptiveDOE(study); %#ok<NASGU>
```

- [ ] **Step 8: Run the resume smoke test**

Run the command from Step 2. Expected: six unique cases after resume, the first four unchanged, and the final two selected in optimization mode.

- [ ] **Step 9: Commit orchestration**

```powershell
git add -- "Full Car Models/runAdaptiveDOE.m" "Full Car Models/SteadyStateLapsim.m" "Full Car Models/tests/test_adaptive_doe_resume.m"
git commit -m "feat: orchestrate resumable adaptive DOE batches"
```

---

### Task 8: Cached analysis and GP-versus-quadratic validation

**Files:**
- Modify: `Full Car Models/utilities/doeAnalyze.m`
- Create: `Full Car Models/utilities/predictDOEModel.m`
- Modify: `Full Car Models/DOE_Fitting.m`
- Create: `Full Car Models/tests/test_doe_model_selection.m`

**Interfaces:**
- `doeAnalyze` preserves `analysis.models.(response)` as the existing stepwise linear/quadratic model.
- Adds `analysis.gpModels.(response)`, `analysis.validation`, `analysis.preferredModel`, and `analysis.predictorBounds` with physical lower and upper values for every predictor.
- Produces: `[y,sd] = predictDOEModel(analysis,response,predictorTable)` using the selected model.
- Loads `metricTable` from the result file when present; otherwise uses `doeMetrics` for legacy files.

- [ ] **Step 1: Write a failing synthetic model-selection test**

```matlab
%% test_doe_model_selection
setup_paths
rng(3); n=80;
X=table(rand(n,1),rand(n,1),'VariableNames',{'x1','x2'});
M=table((1:n)',true(n,1),repmat("",n,1), ...
    45+3*sin(4*X.x1)+2*X.x2.^2, ...
    'VariableNames',{'case_index','valid','error_message','t_autox'});
d=tempname; mkdir(d); cleaner=onCleanup(@() rmdir(d,'s'));
carCell=cell(n,2); designTable=X; metricTable=M; %#ok<NASGU>
path=fullfile(d,'DOE_results.mat');
save(path,'carCell','designTable','metricTable')
A=doeAnalyze(path,struct('responsesWanted',{{'t_autox'}}, ...
    'plots',strings(0,1),'visible','off'));
assert(isfield(A.gpModels,'t_autox'))
assert(any(A.preferredModel.t_autox==["quadratic","gp"]))
[yp,sd]=predictDOEModel(A,'t_autox',X(1:5,:));
assert(numel(yp)==5 && numel(sd)==5 && all(isfinite(yp)))
assert(A.settings.usedCachedMetrics)
```

- [ ] **Step 2: Run the test and verify missing GP analysis fields**

```powershell
matlab -batch "cd('Full Car Models'); run('tests/test_doe_model_selection.m')"
```

- [ ] **Step 3: Prefer cached metric rows**

Load `metricTable` if it exists and its height matches `designTable`. Set `settings.usedCachedMetrics=true`. Otherwise retain the current extraction behavior and set it false.

- [ ] **Step 4: Fit and cross-validate both model families**

Keep current BIC stepwise models. Fit one ARD Matérn GP per response. Use a repeatable five-fold partition, reduced to the number of finite rows when fewer than five. Store RMSE, normalized RMSE, MAE, and rank correlation for both model families.

Store each predictor's physical minimum and maximum from the clean design table in `analysis.predictorBounds`; Sobol sampling uses these values to convert normalized Monte Carlo rows back to model input units.

Choose GP only when its cross-validated normalized RMSE is at least 5% lower than the quadratic model; otherwise prefer the simpler quadratic model.

- [ ] **Step 5: Implement common prediction**

For a preferred GP, return `predict(gp,T{:,:})` mean and standard deviation. For a preferred stepwise model, return `predict(mdl,T)` and a NaN standard-deviation vector.

- [ ] **Step 6: Add fitting configuration to `DOE_Fitting.m`**

Expose `fitGaussianProcesses=true`, selected responses, and checkpoint/result paths without rerunning simulations. Include the three ramp responses when available.

- [ ] **Step 7: Run focused and existing analysis tests**

```powershell
matlab -batch "cd('Full Car Models'); run('tests/test_doe_model_selection.m'); run('tests/test_doe_metrics.m')"
```

- [ ] **Step 8: Commit model comparison**

```powershell
git add -- "Full Car Models/utilities/doeAnalyze.m" "Full Car Models/utilities/predictDOEModel.m" "Full Car Models/DOE_Fitting.m" "Full Car Models/tests/test_doe_model_selection.m"
git commit -m "feat: compare GP and quadratic DOE models"
```

---

### Task 9: Surrogate Sobol sensitivity, plots, and full verification

**Files:**
- Create: `Full Car Models/utilities/doeSobol.m`
- Modify: `Full Car Models/utilities/doeAnalyze.m`
- Modify: `Full Car Models/utilities/doePlotCatalog.m`
- Modify: `Full Car Models/utilities/plotDOEMetrics.m`
- Modify: `Full Car Models/DOE_Fitting.m`
- Create: `Full Car Models/tests/test_doe_sobol.m`

**Interfaces:**
- Produces: `S = doeSobol(analysis,response,nSamples,seed)`.
- `S` is a table with `parameter`, `firstOrder`, `totalOrder`, and bootstrap uncertainty columns.
- `analysis.sobol.(response)` stores the result for each selected response.

- [ ] **Step 1: Write a failing analytical Sobol test**

```matlab
%% test_doe_sobol
setup_paths
rng(4); n=250;
X=table(rand(n,1),rand(n,1),'VariableNames',{'x1','x2'});
M=table((1:n)',true(n,1),repmat("",n,1),X.x1+0.1*X.x2, ...
    'VariableNames',{'case_index','valid','error_message','response'});
d=tempname; mkdir(d); cleaner=onCleanup(@() rmdir(d,'s'));
carCell=cell(n,2); designTable=X; metricTable=M; %#ok<NASGU>
path=fullfile(d,'DOE_results.mat'); save(path,'carCell','designTable','metricTable')
A=doeAnalyze(path,struct('responsesWanted',{{'response'}}, ...
    'plots',strings(0,1),'visible','off'));
S=doeSobol(A,'response',20000,9);
assert(S.firstOrder(S.parameter=="x1") > 0.95)
assert(S.totalOrder(S.parameter=="x1") > S.totalOrder(S.parameter=="x2"))
assert(abs(sum(S.firstOrder)-1) < 0.08)
```

- [ ] **Step 2: Run the test and verify `doeSobol` is missing**

```powershell
matlab -batch "cd('Full Car Models'); run('tests/test_doe_sobol.m')"
```

- [ ] **Step 3: Implement Jansen/Saltelli surrogate sampling**

Generate repeatable independent normalized matrices `A` and `B`, map them to physical predictor tables using `analysis.predictorBounds`, construct one physical `A_Bi` table per predictor, predict through `predictDOEModel`, and calculate:

```matlab
V = var([fA;fB],1);
firstOrder(i) = 1 - mean((fB-fABi).^2)/(2*V);
totalOrder(i) = mean((fA-fABi).^2)/(2*V);
```

Clamp only tiny numerical excursions outside `[0,1]`; retain larger excursions and flag them in the output. Bootstrap rows with the configured seed to generate `firstOrderLow`, `firstOrderHigh`, `totalOrderLow`, and `totalOrderHigh` columns.

- [ ] **Step 4: Add Sobol analysis and figures**

Add `sobol` to `doePlotCatalog`. Produce grouped first-order/total-order bars for points, event times, energy, understeer gradients, and rebalance speed. Label the surrogate family and cross-validation error on each panel. Continue using `setPlotFont` and `saveFigures` so the existing font and output path rules apply.

- [ ] **Step 5: Run all focused DOE tests**

```powershell
matlab -batch "cd('Full Car Models'); setup_paths; run('tests/test_doe_study_config.m'); run('tests/test_doe_explicit_cars.m'); run('tests/test_doe_case_scoring.m'); run('tests/test_doe_adaptive_sampling.m'); run('tests/test_doe_checkpoint.m'); run('tests/test_doe_parallel_batch.m'); run('tests/test_adaptive_doe_resume.m'); run('tests/test_doe_model_selection.m'); run('tests/test_doe_sobol.m'); run('tests/test_doe_metrics.m')"
```

Expected: every script exits successfully with no incomplete or failed assertion.

- [ ] **Step 6: Run MATLAB static checks on every changed production file**

```powershell
matlab -batch "cd('Full Car Models'); setup_paths; F={'DOEStudyConfig.m','runAdaptiveDOE.m','SteadyStateLapsim.m','carConfig.m'}; F=[F,cellstr(string(fullfile('utilities',{'doeResolveStudy.m','doeInitialDesign.m','doeFitSurrogates.m','doeExpectedImprovement.m','doeSelectAdaptiveBatch.m','doeScoreCase.m','doeCaseMetrics.m','doeRunCase.m','doeRunBatch.m','saveDOECheckpoint.m','loadDOECheckpoint.m','doeAnalyze.m','predictDOEModel.m','doeSobol.m'})))]; bad=0; for i=1:numel(F), m=checkcode(F{i},'-id'); bad=bad+sum(ismember({m.id},{'PARSE_ERROR','NASGU','NODEF'})); end; assert(bad==0)"
```

- [ ] **Step 7: Benchmark a small production-shaped batch**

Run 16 LHS cases with two ramp speeds on up to 16 workers. Record total case time, ramp time per car, valid-case count, g-g solve failures, checkpoint size, and estimated 512-case runtime in the simulation log. Verify a plotting-only `DOE_Fitting` run reads cached metrics and does not execute `rampSweep`.

- [ ] **Step 8: Review the final diff for physics-model isolation and user edits**

```powershell
git diff --check
git status --short
git diff -- "Full Car Models/carComponents/Car.m" "Full Car Models/carComponents/Tire2.m"
```

Expected: no adaptive-DOE task changes in `Car.m` or `Tire2.m`; unrelated pre-existing changes remain present and unmodified.

- [ ] **Step 9: Commit Sobol reporting and verification updates**

```powershell
git add -- "Full Car Models/utilities/doeSobol.m" "Full Car Models/utilities/doeAnalyze.m" "Full Car Models/utilities/doePlotCatalog.m" "Full Car Models/utilities/plotDOEMetrics.m" "Full Car Models/DOE_Fitting.m" "Full Car Models/tests/test_doe_sobol.m"
git commit -m "feat: add surrogate Sobol sensitivity reporting"
```

- [ ] **Step 10: Request final code review before merging or running 512 cases**

Use `superpowers:requesting-code-review`, address every correctness finding, rerun Steps 5-8, and only then recommend the production study.
