function cfg = vdTireConfig(modelType,uncertaintyMode,loadParameters)
%VDTIRECONFIG Shared tire settings and explicit repository asset locations.
% The default preserves the existing legacy model and calibration.
% Pass false as the third argument when only settings or asset paths are
% needed. Historical models then do not depend on the full-car fit files.
if nargin < 1 || isempty(modelType), modelType = 'legacy'; end
if nargin < 2 || isempty(uncertaintyMode), uncertaintyMode = 'nominal'; end
if nargin < 3, loadParameters = true; end
fullCarRoot = fileparts(fileparts(mfilename('fullpath')));
repoRoot = fileparts(fullCarRoot);
components = fullfile(fullCarRoot,'carComponents');
cfg = struct();
cfg.gamma_f = -1;
cfg.gamma_r = -1;
cfg.p_i = 12; % psi (native fit pressure unit)
cfg.model_type = lower(char(modelType));
cfg.model_uncertainty = lower(char(uncertaintyMode));
cfg.model_artifact = fullfile(repoRoot,'Magic Formula','experimental', ...
    'lc0_nd_tire','results','lc0_nd_pacejka_comparison.mat');
cfg.Fx_file = fullfile(components,'Fx_combined_parameters_run38_30.mat');
cfg.Fy_file = fullfile(components,'Lapsim_Fy_combined_parameters_1965run15.mat');
cfg.camber_ratio_file = fullfile(components,'camberratiossmoothed.mat');
cfg.camber_model_file = fullfile(components,'camber_models_fast.mat');
if loadParameters
    cfg.Fx_parameters = vdLoadTireParameters(cfg.Fx_file,26);
    cfg.Fy_parameters = vdLoadTireParameters(cfg.Fy_file,42);
end
cfg.friction_scaling_factor = 1;
cfg.grip_scaling_front = 0.6125;
cfg.grip_scaling_rear = 0.62;
% Existing vertical stiffness choices for the suspension studies (N/m).
cfg.vertical_stiffness.modal_spring_multiplier = 2;
cfg.vertical_stiffness.transient_N_per_m = 724/0.0057101;
cfg.vertical_stiffness.vibration_spring_multiplier = 10;

% Existing alternate fit formats for historical standalone studies.
cfg.legacyFits.fit1 = fullfile(repoRoot,'Bicycle Model','utilities','Fy_pure_parameters_run24_new2.mat');
cfg.legacyFits.fit2 = fullfile(repoRoot,'Magic Formula','Parameters','Fy_pure_parameters_run24_17.mat');
cfg.legacyFits.fit3 = fullfile(repoRoot,'Magic Formula','Parameters','Fy_pure_parameters_run24_18.mat');
cfg.legacyFits.fit4 = fullfile(repoRoot,'Full Car Models','carComponents','Fy_combined_parameters_16x10x7.5_LC0.mat');
cfg.legacyFits.fit5 = fullfile(repoRoot,'Full Car Models','CamberCurveGen','Fy_pure_parameters_1965run24.mat');

% Standalone historical studies retain their original inputs and native units.
% They are separate named presets, not a claim that they describe one car.
% Bicycle Model/toe_compliance.m
cfg.studies.Bicycle_Model_toe_compliance.fit_1 = cfg.legacyFits.fit1;
% Bicycle Model/BicycleModel.m
cfg.studies.Bicycle_Model_BicycleModel.fit_1 = cfg.legacyFits.fit1;
% Bicycle Model/Simulink Models/BicycleModelSimv2.m
cfg.studies.Bicycle_Model_Simulink_Models_BicycleModelSimv2.C_f = 136*4.448*180/pi/2;
cfg.studies.Bicycle_Model_Simulink_Models_BicycleModelSimv2.C_r = 136*4.448*180/pi/2;
% Bicycle Model/Simulink Models/NonLin_BicycleModel.m
cfg.studies.Bicycle_Model_Simulink_Models_NonLin_BicycleModel.C_f = 136*4.448*180/pi/2;
cfg.studies.Bicycle_Model_Simulink_Models_NonLin_BicycleModel.C_r = 136*4.448*180/pi/2;
cfg.studies.Bicycle_Model_Simulink_Models_NonLin_BicycleModel.fit_1 = cfg.legacyFits.fit2;
% Bicycle Model/Simulink Models/FullCarModel.m
cfg.studies.Bicycle_Model_Simulink_Models_FullCarModel.C_f = 136*4.448*180/pi/2;
cfg.studies.Bicycle_Model_Simulink_Models_FullCarModel.C_r = 136*4.448*180/pi/2;
cfg.studies.Bicycle_Model_Simulink_Models_FullCarModel.fit_1 = cfg.legacyFits.fit3;
% Bicycle Model/Simulink Models/BicycleModelSim.m
cfg.studies.Bicycle_Model_Simulink_Models_BicycleModelSim.C_f = 55000;
cfg.studies.Bicycle_Model_Simulink_Models_BicycleModelSim.C_r = 60000;
% Bicycle Model/ODE Formulation/BicycleModelODE.m
cfg.studies.Bicycle_Model_ODE_Formulation_BicycleModelODE.C_f = 136*4.448*180/pi;
cfg.studies.Bicycle_Model_ODE_Formulation_BicycleModelODE.C_r = 136*4.448*180/pi;
cfg.studies.Bicycle_Model_ODE_Formulation_BicycleModelODE.fit_1 = cfg.legacyFits.fit2;
% understeer_bode.m
cfg.studies.understeer_bode.Csum = 1.6e5;
% understeer_budget.m
cfg.studies.understeer_budget.CAF = 136*4.448;
cfg.studies.understeer_budget.CCF = 0.313;
cfg.studies.understeer_budget.NAF = 50.082;
cfg.studies.understeer_budget.CGF = 47.657;
cfg.studies.understeer_budget.NGF = 4.766;
% Full Car Models/CamberCurveGen/generateIdealCamberCurves.m
cfg.studies.camber1_generateIdealCamberCurves.friction_scaling_factor = 1.05*0.55;
cfg.studies.camber1_generateIdealCamberCurves.p_i = cfg.p_i;
cfg.studies.camber1_generateIdealCamberCurves.fit_1 = cfg.Fx_file;
cfg.studies.camber1_generateIdealCamberCurves.fit_2 = cfg.legacyFits.fit4;
% Full Car Models/CamberCurveGen/test.m
cfg.studies.camber1_test.friction_scaling_factor = 1.05*0.55;
cfg.studies.camber1_test.p_i = cfg.p_i;
cfg.studies.camber1_test.fit_1 = cfg.Fx_file;
cfg.studies.camber1_test.fit_2 = cfg.legacyFits.fit4;
% Full Car Models/CamberCurveGen/Testcamber.m
cfg.studies.camber1_Testcamber.friction_scaling_factor = 1.05*0.55;
cfg.studies.camber1_Testcamber.p_i = cfg.p_i;
cfg.studies.camber1_Testcamber.fit_1 = cfg.Fx_file;
cfg.studies.camber1_Testcamber.fit_2 = cfg.legacyFits.fit5;
% Full Car Models/CamberCurveGen/generateIdealCamberCurves_simplified.m
cfg.studies.camber1_generateIdealCamberCurves_simplified.friction_scaling_factor = 1.05*0.55;
cfg.studies.camber1_generateIdealCamberCurves_simplified.p_i = cfg.p_i;
cfg.studies.camber1_generateIdealCamberCurves_simplified.fit_1 = cfg.Fx_file;
cfg.studies.camber1_generateIdealCamberCurves_simplified.fit_2 = cfg.legacyFits.fit4;
% Full Car Models/CamberCurveGen2/generateIdealCamberCurves.m
cfg.studies.camber2_generateIdealCamberCurves.friction_scaling_factor = 1.05*0.55;
cfg.studies.camber2_generateIdealCamberCurves.p_i = cfg.p_i;
cfg.studies.camber2_generateIdealCamberCurves.fit_1 = cfg.Fx_file;
cfg.studies.camber2_generateIdealCamberCurves.fit_2 = cfg.legacyFits.fit4;
% Full Car Models/CamberCurveGen2/test.m
cfg.studies.camber2_test.friction_scaling_factor = 1.05*0.55;
cfg.studies.camber2_test.p_i = cfg.p_i;
cfg.studies.camber2_test.fit_1 = cfg.Fx_file;
cfg.studies.camber2_test.fit_2 = cfg.legacyFits.fit4;
% Full Car Models/CamberCurveGen2/Testcamber.m
cfg.studies.camber2_Testcamber.friction_scaling_factor = 1.05*0.55;
cfg.studies.camber2_Testcamber.p_i = cfg.p_i;
cfg.studies.camber2_Testcamber.fit_1 = cfg.Fx_file;
cfg.studies.camber2_Testcamber.fit_2 = cfg.legacyFits.fit4;
% Full Car Models/CamberCurveGen2/generateIdealCamberCurves_simplified.m
cfg.studies.camber2_generateIdealCamberCurves_simplified.friction_scaling_factor = 1.05*0.55;
cfg.studies.camber2_generateIdealCamberCurves_simplified.p_i = cfg.p_i;
cfg.studies.camber2_generateIdealCamberCurves_simplified.fit_1 = cfg.Fx_file;
cfg.studies.camber2_generateIdealCamberCurves_simplified.fit_2 = cfg.legacyFits.fit4;

% Remaining paired and embedded inputs from the same historical studies.
cfg.studies.understeer_budget.CAR = 136*4.448;
cfg.studies.understeer_budget.CCR = 0.323;
cfg.studies.understeer_budget.NAR = 42.897;
cfg.studies.understeer_budget.CGR = 43.392;
cfg.studies.understeer_budget.NGR = 4.339;
cfg.studies.understeer_bode.K_deg_per_g = 0.18;
cfg.studies.Bicycle_Model_BicycleModel.cornering_force_scale = 0.55;
cfg.studies.Bicycle_Model_BicycleModel.cornering_load_mass_lb = 535;
cfg.studies.Bicycle_Model_BicycleModel.cornering_pressure_psi = 12;
cfg.studies.Bicycle_Model_toe_compliance.cornering_force_scale = 0.55;
cfg.studies.Bicycle_Model_toe_compliance.cornering_load_mass_lb = 535;
cfg.studies.Bicycle_Model_toe_compliance.cornering_pressure_psi = 12;

end
