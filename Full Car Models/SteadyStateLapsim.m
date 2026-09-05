%% SteadyStateLapsim - adaptive DOE runner
clear classes
setup_paths
study = DOEStudyConfig();

study.name = "design_sensitivity_small_v1";
study.mode = "sensitivity";
study.randomSeed = 13;

study.initialCases = 32;   % space-filling initial design
study.batchSize = 8;       % four adaptive batches after the initial set
study.maxCases = 64;       % small directional study
study.numWorkers = 8;      % set 0 for serial
study.resume = true;

study.events = ["skidpad","accel","autocross"];

% Keep the first study fast; no balance ramps yet.
study.ramps.enabled = false;

% Only select on responses available without ramps.
study.adaptive.responses = { ...
    'modeled_dynamic_points','t_autox','t_accel','t_skid','total_work_kJ'};

study.output.directory = fullfile( ...
    fileparts(which('DOEStudyConfig')), ...
    "DOE_output_design_sensitivity_small_v1");

state = runAdaptiveDOE(study);
