%% test_doe_case_scoring
testDir = fileparts(mfilename('fullpath'));
modelRoot = fileparts(testDir);
cd(modelRoot); setup_paths
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

cfg.penalties.understeer.enabled = true;
cfg.penalties.understeer.lower_deg_per_g = 1;
cfg.penalties.understeer.upper_deg_per_g = 2;
cfg.penalties.understeer.pointsPerDegPerG = 4;
[~,B] = doeScoreCase(points,M,cfg);
assert(abs(B.penalty_understeer-2) < 1e-12)
cfg.penalties.rebalance.enabled = true;
cfg.penalties.rebalance.target_mps = 16;
cfg.penalties.rebalance.pointsPerMps = 3;
[~,B] = doeScoreCase(points,M,cfg);
assert(abs(B.penalty_rebalance-6) < 1e-12)

emptyMetrics = doeMetrics(cell(0,1));
invalidMetrics = doeMetrics({[]});
assert(height(emptyMetrics) == 0)
assert(isequal(emptyMetrics,invalidMetrics([],:)))
