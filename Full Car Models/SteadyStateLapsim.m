%% SteadyStateLapsim - adaptive DOE runner
clear classes
setup_paths
study = DOEStudyConfig();
state = runAdaptiveDOE(study); %#ok<NASGU>
