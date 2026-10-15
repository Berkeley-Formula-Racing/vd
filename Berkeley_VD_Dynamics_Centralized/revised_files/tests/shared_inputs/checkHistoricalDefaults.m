function count = checkHistoricalDefaults(cfg,tireCfg)
% Compare migrated standalone literals with upstream 0b02ea46 defaults.
assert(abs(cfg.studies.Bicycle_Model_BicycleModel.weight_dist-(0.54)) < 1e-10); % Bicycle Model/BicycleModel.m: weight_dist
assert(abs(cfg.studies.Bicycle_Model_BicycleModel.m-(535*0.453592)) < 1e-10); % Bicycle Model/BicycleModel.m: m
assert(abs(cfg.studies.Bicycle_Model_BicycleModel.L-(1.6)) < 1e-10); % Bicycle Model/BicycleModel.m: L
assert(abs(cfg.studies.Bicycle_Model_BicycleModel.I_zz-(83.28)) < 1e-10); % Bicycle Model/BicycleModel.m: I_zz
assert(abs(cfg.studies.Bicycle_Model_BicycleModel.N-(4)) < 1e-10); % Bicycle Model/BicycleModel.m: N
assert(abs(cfg.studies.Bicycle_Model_LatAccelResponseTimePlotter.weight_dist-(0.51)) < 1e-10); % Bicycle Model/LatAccelResponseTimePlotter.m: weight_dist
assert(abs(cfg.studies.Bicycle_Model_LatAccelResponseTimePlotter.m-(535*0.453592)) < 1e-10); % Bicycle Model/LatAccelResponseTimePlotter.m: m
assert(abs(cfg.studies.Bicycle_Model_LatAccelResponseTimePlotter.L-(1.524)) < 1e-10); % Bicycle Model/LatAccelResponseTimePlotter.m: L
assert(abs(cfg.studies.Bicycle_Model_LatAccelResponseTimePlotter.I_zz-(83.28)) < 1e-10); % Bicycle Model/LatAccelResponseTimePlotter.m: I_zz
assert(abs(cfg.studies.Bicycle_Model_LatAccelResponseTimePlotter.N-(4)) < 1e-10); % Bicycle Model/LatAccelResponseTimePlotter.m: N
assert(abs(cfg.studies.Bicycle_Model_ODE_Formulation_BicycleModelODE.weight_dist-(0.3)) < 1e-10); % Bicycle Model/ODE Formulation/BicycleModelODE.m: weight_dist
assert(abs(cfg.studies.Bicycle_Model_ODE_Formulation_BicycleModelODE.m-(242.6718)) < 1e-10); % Bicycle Model/ODE Formulation/BicycleModelODE.m: m
assert(abs(cfg.studies.Bicycle_Model_ODE_Formulation_BicycleModelODE.L-(1.524)) < 1e-10); % Bicycle Model/ODE Formulation/BicycleModelODE.m: L
assert(abs(cfg.studies.Bicycle_Model_ODE_Formulation_BicycleModelODE.I_zz-(83.28)) < 1e-10); % Bicycle Model/ODE Formulation/BicycleModelODE.m: I_zz
assert(abs(tireCfg.studies.Bicycle_Model_ODE_Formulation_BicycleModelODE.C_f-(136*4.448*180/pi)) < 1e-10); % Bicycle Model/ODE Formulation/BicycleModelODE.m: C_f
assert(abs(tireCfg.studies.Bicycle_Model_ODE_Formulation_BicycleModelODE.C_r-(136*4.448*180/pi)) < 1e-10); % Bicycle Model/ODE Formulation/BicycleModelODE.m: C_r
assert(abs(cfg.studies.Bicycle_Model_ODE_Formulation_BicycleModelODE.N-(16)) < 1e-10); % Bicycle Model/ODE Formulation/BicycleModelODE.m: N
assert(abs(cfg.studies.Bicycle_Model_Simulink_Models_BicycleModelSim.m-(1500)) < 1e-10); % Bicycle Model/Simulink Models/BicycleModelSim.m: m
assert(abs(cfg.studies.Bicycle_Model_Simulink_Models_BicycleModelSim.L-(2500)) < 1e-10); % Bicycle Model/Simulink Models/BicycleModelSim.m: L
assert(abs(cfg.studies.Bicycle_Model_Simulink_Models_BicycleModelSim.I_zz-(2500)) < 1e-10); % Bicycle Model/Simulink Models/BicycleModelSim.m: I_zz
assert(abs(tireCfg.studies.Bicycle_Model_Simulink_Models_BicycleModelSim.C_f-(55000)) < 1e-10); % Bicycle Model/Simulink Models/BicycleModelSim.m: C_f
assert(abs(tireCfg.studies.Bicycle_Model_Simulink_Models_BicycleModelSim.C_r-(60000)) < 1e-10); % Bicycle Model/Simulink Models/BicycleModelSim.m: C_r
assert(abs(cfg.studies.Bicycle_Model_Simulink_Models_BicycleModelSimv2.weight_dist-(0.3)) < 1e-10); % Bicycle Model/Simulink Models/BicycleModelSimv2.m: weight_dist
assert(abs(cfg.studies.Bicycle_Model_Simulink_Models_BicycleModelSimv2.m-(242.6718)) < 1e-10); % Bicycle Model/Simulink Models/BicycleModelSimv2.m: m
assert(abs(cfg.studies.Bicycle_Model_Simulink_Models_BicycleModelSimv2.L-(1.524)) < 1e-10); % Bicycle Model/Simulink Models/BicycleModelSimv2.m: L
assert(abs(cfg.studies.Bicycle_Model_Simulink_Models_BicycleModelSimv2.I_zz-(83.28)) < 1e-10); % Bicycle Model/Simulink Models/BicycleModelSimv2.m: I_zz
assert(abs(tireCfg.studies.Bicycle_Model_Simulink_Models_BicycleModelSimv2.C_f-(136*4.448*180/pi/2)) < 1e-10); % Bicycle Model/Simulink Models/BicycleModelSimv2.m: C_f
assert(abs(tireCfg.studies.Bicycle_Model_Simulink_Models_BicycleModelSimv2.C_r-(136*4.448*180/pi/2)) < 1e-10); % Bicycle Model/Simulink Models/BicycleModelSimv2.m: C_r
assert(abs(cfg.studies.Bicycle_Model_Simulink_Models_FullCarModel.weight_dist-(0.6)) < 1e-10); % Bicycle Model/Simulink Models/FullCarModel.m: weight_dist
assert(abs(cfg.studies.Bicycle_Model_Simulink_Models_FullCarModel.m-(242.6718)) < 1e-10); % Bicycle Model/Simulink Models/FullCarModel.m: m
assert(abs(cfg.studies.Bicycle_Model_Simulink_Models_FullCarModel.L-(1.524)) < 1e-10); % Bicycle Model/Simulink Models/FullCarModel.m: L
assert(abs(cfg.studies.Bicycle_Model_Simulink_Models_FullCarModel.I_zz-(83.28)) < 1e-10); % Bicycle Model/Simulink Models/FullCarModel.m: I_zz
assert(abs(tireCfg.studies.Bicycle_Model_Simulink_Models_FullCarModel.C_f-(136*4.448*180/pi/2)) < 1e-10); % Bicycle Model/Simulink Models/FullCarModel.m: C_f
assert(abs(tireCfg.studies.Bicycle_Model_Simulink_Models_FullCarModel.C_r-(136*4.448*180/pi/2)) < 1e-10); % Bicycle Model/Simulink Models/FullCarModel.m: C_r
assert(abs(cfg.studies.Bicycle_Model_Simulink_Models_FullCarModel.track_width_f-(1.1938)) < 1e-10); % Bicycle Model/Simulink Models/FullCarModel.m: track_width_f
assert(abs(cfg.studies.Bicycle_Model_Simulink_Models_FullCarModel.track_width_r-(1.1938)) < 1e-10); % Bicycle Model/Simulink Models/FullCarModel.m: track_width_r
assert(abs(cfg.studies.Bicycle_Model_Simulink_Models_FullCarModel.cla-(1.771)) < 1e-10); % Bicycle Model/Simulink Models/FullCarModel.m: cla
assert(abs(cfg.studies.Bicycle_Model_Simulink_Models_FullCarModel.cp-(0.5)) < 1e-10); % Bicycle Model/Simulink Models/FullCarModel.m: cp
assert(abs(cfg.studies.Bicycle_Model_Simulink_Models_FullCarModel.LLTD-(0.6)) < 1e-10); % Bicycle Model/Simulink Models/FullCarModel.m: LLTD
assert(abs(cfg.studies.Bicycle_Model_Simulink_Models_FullCarModel.cg_height-(0.3048)) < 1e-10); % Bicycle Model/Simulink Models/FullCarModel.m: cg_height
assert(abs(cfg.studies.Bicycle_Model_Simulink_Models_NonLin_BicycleModel.weight_dist-(0.55)) < 1e-10); % Bicycle Model/Simulink Models/NonLin_BicycleModel.m: weight_dist
assert(abs(cfg.studies.Bicycle_Model_Simulink_Models_NonLin_BicycleModel.m-(242.6718)) < 1e-10); % Bicycle Model/Simulink Models/NonLin_BicycleModel.m: m
assert(abs(cfg.studies.Bicycle_Model_Simulink_Models_NonLin_BicycleModel.L-(1.524)) < 1e-10); % Bicycle Model/Simulink Models/NonLin_BicycleModel.m: L
assert(abs(cfg.studies.Bicycle_Model_Simulink_Models_NonLin_BicycleModel.I_zz-(83.28)) < 1e-10); % Bicycle Model/Simulink Models/NonLin_BicycleModel.m: I_zz
assert(abs(tireCfg.studies.Bicycle_Model_Simulink_Models_NonLin_BicycleModel.C_f-(136*4.448*180/pi/2)) < 1e-10); % Bicycle Model/Simulink Models/NonLin_BicycleModel.m: C_f
assert(abs(tireCfg.studies.Bicycle_Model_Simulink_Models_NonLin_BicycleModel.C_r-(136*4.448*180/pi/2)) < 1e-10); % Bicycle Model/Simulink Models/NonLin_BicycleModel.m: C_r
assert(abs(cfg.studies.Bicycle_Model_YawResponseTimePlotter.weight_dist-(0.51)) < 1e-10); % Bicycle Model/YawResponseTimePlotter.m: weight_dist
assert(abs(cfg.studies.Bicycle_Model_YawResponseTimePlotter.m-(535*0.453592)) < 1e-10); % Bicycle Model/YawResponseTimePlotter.m: m
assert(abs(cfg.studies.Bicycle_Model_YawResponseTimePlotter.L-(1.524)) < 1e-10); % Bicycle Model/YawResponseTimePlotter.m: L
assert(abs(cfg.studies.Bicycle_Model_YawResponseTimePlotter.I_zz-(83.28)) < 1e-10); % Bicycle Model/YawResponseTimePlotter.m: I_zz
assert(abs(cfg.studies.Bicycle_Model_YawResponseTimePlotter.N-(4)) < 1e-10); % Bicycle Model/YawResponseTimePlotter.m: N
assert(abs(cfg.studies.Bicycle_Model_toe_compliance.weight_dist-(0.51)) < 1e-10); % Bicycle Model/toe_compliance.m: weight_dist
assert(abs(cfg.studies.Bicycle_Model_toe_compliance.m-(520*0.453592)) < 1e-10); % Bicycle Model/toe_compliance.m: m
assert(abs(cfg.studies.Bicycle_Model_toe_compliance.L-(62 * 0.0254)) < 1e-10); % Bicycle Model/toe_compliance.m: L
assert(abs(cfg.studies.Bicycle_Model_toe_compliance.I_zz-(83.28)) < 1e-10); % Bicycle Model/toe_compliance.m: I_zz
assert(abs(cfg.studies.Bicycle_Model_toe_compliance.N-(4)) < 1e-10); % Bicycle Model/toe_compliance.m: N
assert(abs(tireCfg.studies.camber1_Testcamber.friction_scaling_factor-(1.05*0.55)) < 1e-10); % Full Car Models/CamberCurveGen/Testcamber.m: friction_scaling_factor
assert(abs(tireCfg.studies.camber1_Testcamber.p_i-(12)) < 1e-10); % Full Car Models/CamberCurveGen/Testcamber.m: p_i
assert(abs(cfg.studies.camber1_Testcamber.C_front_width-(47)) < 1e-10); % Full Car Models/CamberCurveGen/Testcamber.m: C.front_width
assert(abs(cfg.studies.camber1_Testcamber.C_rear_width-(47)) < 1e-10); % Full Car Models/CamberCurveGen/Testcamber.m: C.rear_width
assert(abs(cfg.studies.camber1_Testcamber.C_wheelbase-(62)) < 1e-10); % Full Car Models/CamberCurveGen/Testcamber.m: C.wheelbase
assert(abs(cfg.studies.camber1_Testcamber.C_front_spring_roll_stiffness-(2860)) < 1e-10); % Full Car Models/CamberCurveGen/Testcamber.m: C.front_spring_roll_stiffness
assert(abs(cfg.studies.camber1_Testcamber.C_rear_spring_roll_stiffness-(2911)) < 1e-10); % Full Car Models/CamberCurveGen/Testcamber.m: C.rear_spring_roll_stiffness
assert(abs(cfg.studies.camber1_Testcamber.C_front_bump_stiffness-(2*191.66)) < 1e-10); % Full Car Models/CamberCurveGen/Testcamber.m: C.front_bump_stiffness
assert(abs(cfg.studies.camber1_Testcamber.C_rear_bump_stiffness-(2*196.02)) < 1e-10); % Full Car Models/CamberCurveGen/Testcamber.m: C.rear_bump_stiffness
assert(abs(cfg.studies.camber1_Testcamber.C_front_ARB_roll_stiffness-(0)) < 1e-10); % Full Car Models/CamberCurveGen/Testcamber.m: C.front_ARB_roll_stiffness
assert(abs(cfg.studies.camber1_Testcamber.C_rear_ARB_roll_stiffness-(986)) < 1e-10); % Full Car Models/CamberCurveGen/Testcamber.m: C.rear_ARB_roll_stiffness
assert(abs(cfg.studies.camber1_Testcamber.C_weight_dist-(0.54)) < 1e-10); % Full Car Models/CamberCurveGen/Testcamber.m: C.weight_dist
assert(abs(cfg.studies.camber1_Testcamber.C_mass-((395 + 150)*4.44)) < 1e-10); % Full Car Models/CamberCurveGen/Testcamber.m: C.mass
assert(abs(cfg.studies.camber1_Testcamber.front_anti_dive-(.545)) < 1e-10); % Full Car Models/CamberCurveGen/Testcamber.m: front_anti_dive
assert(abs(cfg.studies.camber1_Testcamber.rear_anti_lift-(1.02)) < 1e-10); % Full Car Models/CamberCurveGen/Testcamber.m: rear_anti_lift
assert(abs(cfg.studies.camber1_Testcamber.rear_anti_squat-(.457)) < 1e-10); % Full Car Models/CamberCurveGen/Testcamber.m: rear_anti_squat
assert(abs(tireCfg.studies.camber1_generateIdealCamberCurves.friction_scaling_factor-(1.05*0.55)) < 1e-10); % Full Car Models/CamberCurveGen/generateIdealCamberCurves.m: friction_scaling_factor
assert(abs(tireCfg.studies.camber1_generateIdealCamberCurves.p_i-(12)) < 1e-10); % Full Car Models/CamberCurveGen/generateIdealCamberCurves.m: p_i
assert(abs(cfg.studies.camber1_generateIdealCamberCurves.C_front_width-(47)) < 1e-10); % Full Car Models/CamberCurveGen/generateIdealCamberCurves.m: C.front_width
assert(abs(cfg.studies.camber1_generateIdealCamberCurves.C_rear_width-(47)) < 1e-10); % Full Car Models/CamberCurveGen/generateIdealCamberCurves.m: C.rear_width
assert(abs(cfg.studies.camber1_generateIdealCamberCurves.C_front_spring_roll_stiffness-(3000)) < 1e-10); % Full Car Models/CamberCurveGen/generateIdealCamberCurves.m: C.front_spring_roll_stiffness
assert(abs(cfg.studies.camber1_generateIdealCamberCurves.C_rear_spring_roll_stiffness-(3000)) < 1e-10); % Full Car Models/CamberCurveGen/generateIdealCamberCurves.m: C.rear_spring_roll_stiffness
assert(abs(cfg.studies.camber1_generateIdealCamberCurves.C_front_ARB_roll_stiffness-(0)) < 1e-10); % Full Car Models/CamberCurveGen/generateIdealCamberCurves.m: C.front_ARB_roll_stiffness
assert(abs(cfg.studies.camber1_generateIdealCamberCurves.C_rear_ARB_roll_stiffness-(0)) < 1e-10); % Full Car Models/CamberCurveGen/generateIdealCamberCurves.m: C.rear_ARB_roll_stiffness
assert(abs(cfg.studies.camber1_generateIdealCamberCurves.C_weight_dist-(0.54)) < 1e-10); % Full Car Models/CamberCurveGen/generateIdealCamberCurves.m: C.weight_dist
assert(abs(cfg.studies.camber1_generateIdealCamberCurves.C_mass-((395 + 150))) < 1e-10); % Full Car Models/CamberCurveGen/generateIdealCamberCurves.m: C.mass
assert(abs(tireCfg.studies.camber1_generateIdealCamberCurves_simplified.friction_scaling_factor-(1.05*0.55)) < 1e-10); % Full Car Models/CamberCurveGen/generateIdealCamberCurves_simplified.m: friction_scaling_factor
assert(abs(tireCfg.studies.camber1_generateIdealCamberCurves_simplified.p_i-(12)) < 1e-10); % Full Car Models/CamberCurveGen/generateIdealCamberCurves_simplified.m: p_i
assert(abs(cfg.studies.camber1_generateIdealCamberCurves_simplified.C_front_width-(47)) < 1e-10); % Full Car Models/CamberCurveGen/generateIdealCamberCurves_simplified.m: C.front_width
assert(abs(cfg.studies.camber1_generateIdealCamberCurves_simplified.C_rear_width-(47)) < 1e-10); % Full Car Models/CamberCurveGen/generateIdealCamberCurves_simplified.m: C.rear_width
assert(abs(cfg.studies.camber1_generateIdealCamberCurves_simplified.C_front_spring_roll_stiffness-(3000)) < 1e-10); % Full Car Models/CamberCurveGen/generateIdealCamberCurves_simplified.m: C.front_spring_roll_stiffness
assert(abs(cfg.studies.camber1_generateIdealCamberCurves_simplified.C_rear_spring_roll_stiffness-(0)) < 1e-10); % Full Car Models/CamberCurveGen/generateIdealCamberCurves_simplified.m: C.rear_spring_roll_stiffness
assert(abs(cfg.studies.camber1_generateIdealCamberCurves_simplified.C_front_ARB_roll_stiffness-(0)) < 1e-10); % Full Car Models/CamberCurveGen/generateIdealCamberCurves_simplified.m: C.front_ARB_roll_stiffness
assert(abs(cfg.studies.camber1_generateIdealCamberCurves_simplified.C_rear_ARB_roll_stiffness-(0)) < 1e-10); % Full Car Models/CamberCurveGen/generateIdealCamberCurves_simplified.m: C.rear_ARB_roll_stiffness
assert(abs(cfg.studies.camber1_generateIdealCamberCurves_simplified.C_weight_dist-(0.54)) < 1e-10); % Full Car Models/CamberCurveGen/generateIdealCamberCurves_simplified.m: C.weight_dist
assert(abs(cfg.studies.camber1_generateIdealCamberCurves_simplified.C_mass-((395 + 150))) < 1e-10); % Full Car Models/CamberCurveGen/generateIdealCamberCurves_simplified.m: C.mass
assert(abs(tireCfg.studies.camber1_test.friction_scaling_factor-(1.05*0.55)) < 1e-10); % Full Car Models/CamberCurveGen/test.m: friction_scaling_factor
assert(abs(tireCfg.studies.camber1_test.p_i-(12)) < 1e-10); % Full Car Models/CamberCurveGen/test.m: p_i
assert(abs(tireCfg.studies.camber2_Testcamber.friction_scaling_factor-(1.05*0.55)) < 1e-10); % Full Car Models/CamberCurveGen2/Testcamber.m: friction_scaling_factor
assert(abs(tireCfg.studies.camber2_Testcamber.p_i-(12)) < 1e-10); % Full Car Models/CamberCurveGen2/Testcamber.m: p_i
assert(abs(cfg.studies.camber2_Testcamber.C_front_width-(47)) < 1e-10); % Full Car Models/CamberCurveGen2/Testcamber.m: C.front_width
assert(abs(cfg.studies.camber2_Testcamber.C_rear_width-(47)) < 1e-10); % Full Car Models/CamberCurveGen2/Testcamber.m: C.rear_width
assert(abs(cfg.studies.camber2_Testcamber.C_front_spring_roll_stiffness-(2860)) < 1e-10); % Full Car Models/CamberCurveGen2/Testcamber.m: C.front_spring_roll_stiffness
assert(abs(cfg.studies.camber2_Testcamber.C_rear_spring_roll_stiffness-(2911)) < 1e-10); % Full Car Models/CamberCurveGen2/Testcamber.m: C.rear_spring_roll_stiffness
assert(abs(cfg.studies.camber2_Testcamber.C_front_ARB_roll_stiffness-(0)) < 1e-10); % Full Car Models/CamberCurveGen2/Testcamber.m: C.front_ARB_roll_stiffness
assert(abs(cfg.studies.camber2_Testcamber.C_rear_ARB_roll_stiffness-(986)) < 1e-10); % Full Car Models/CamberCurveGen2/Testcamber.m: C.rear_ARB_roll_stiffness
assert(abs(cfg.studies.camber2_Testcamber.C_weight_dist-(0.54)) < 1e-10); % Full Car Models/CamberCurveGen2/Testcamber.m: C.weight_dist
assert(abs(cfg.studies.camber2_Testcamber.C_mass-(395 + 150)) < 1e-10); % Full Car Models/CamberCurveGen2/Testcamber.m: C.mass
assert(abs(tireCfg.studies.camber2_generateIdealCamberCurves.friction_scaling_factor-(1.05*0.55)) < 1e-10); % Full Car Models/CamberCurveGen2/generateIdealCamberCurves.m: friction_scaling_factor
assert(abs(tireCfg.studies.camber2_generateIdealCamberCurves.p_i-(12)) < 1e-10); % Full Car Models/CamberCurveGen2/generateIdealCamberCurves.m: p_i
assert(abs(cfg.studies.camber2_generateIdealCamberCurves.C_front_width-(47)) < 1e-10); % Full Car Models/CamberCurveGen2/generateIdealCamberCurves.m: C.front_width
assert(abs(cfg.studies.camber2_generateIdealCamberCurves.C_rear_width-(47)) < 1e-10); % Full Car Models/CamberCurveGen2/generateIdealCamberCurves.m: C.rear_width
assert(abs(cfg.studies.camber2_generateIdealCamberCurves.C_front_spring_roll_stiffness-(3000)) < 1e-10); % Full Car Models/CamberCurveGen2/generateIdealCamberCurves.m: C.front_spring_roll_stiffness
assert(abs(cfg.studies.camber2_generateIdealCamberCurves.C_rear_spring_roll_stiffness-(3000)) < 1e-10); % Full Car Models/CamberCurveGen2/generateIdealCamberCurves.m: C.rear_spring_roll_stiffness
assert(abs(cfg.studies.camber2_generateIdealCamberCurves.C_front_ARB_roll_stiffness-(0)) < 1e-10); % Full Car Models/CamberCurveGen2/generateIdealCamberCurves.m: C.front_ARB_roll_stiffness
assert(abs(cfg.studies.camber2_generateIdealCamberCurves.C_rear_ARB_roll_stiffness-(0)) < 1e-10); % Full Car Models/CamberCurveGen2/generateIdealCamberCurves.m: C.rear_ARB_roll_stiffness
assert(abs(cfg.studies.camber2_generateIdealCamberCurves.C_weight_dist-(0.54)) < 1e-10); % Full Car Models/CamberCurveGen2/generateIdealCamberCurves.m: C.weight_dist
assert(abs(cfg.studies.camber2_generateIdealCamberCurves.C_mass-((395 + 150))) < 1e-10); % Full Car Models/CamberCurveGen2/generateIdealCamberCurves.m: C.mass
assert(abs(tireCfg.studies.camber2_generateIdealCamberCurves_simplified.friction_scaling_factor-(1.05*0.55)) < 1e-10); % Full Car Models/CamberCurveGen2/generateIdealCamberCurves_simplified.m: friction_scaling_factor
assert(abs(tireCfg.studies.camber2_generateIdealCamberCurves_simplified.p_i-(12)) < 1e-10); % Full Car Models/CamberCurveGen2/generateIdealCamberCurves_simplified.m: p_i
assert(abs(cfg.studies.camber2_generateIdealCamberCurves_simplified.C_front_width-(47)) < 1e-10); % Full Car Models/CamberCurveGen2/generateIdealCamberCurves_simplified.m: C.front_width
assert(abs(cfg.studies.camber2_generateIdealCamberCurves_simplified.C_rear_width-(47)) < 1e-10); % Full Car Models/CamberCurveGen2/generateIdealCamberCurves_simplified.m: C.rear_width
assert(abs(cfg.studies.camber2_generateIdealCamberCurves_simplified.C_front_spring_roll_stiffness-(3000)) < 1e-10); % Full Car Models/CamberCurveGen2/generateIdealCamberCurves_simplified.m: C.front_spring_roll_stiffness
assert(abs(cfg.studies.camber2_generateIdealCamberCurves_simplified.C_rear_spring_roll_stiffness-(0)) < 1e-10); % Full Car Models/CamberCurveGen2/generateIdealCamberCurves_simplified.m: C.rear_spring_roll_stiffness
assert(abs(cfg.studies.camber2_generateIdealCamberCurves_simplified.C_front_ARB_roll_stiffness-(0)) < 1e-10); % Full Car Models/CamberCurveGen2/generateIdealCamberCurves_simplified.m: C.front_ARB_roll_stiffness
assert(abs(cfg.studies.camber2_generateIdealCamberCurves_simplified.C_rear_ARB_roll_stiffness-(0)) < 1e-10); % Full Car Models/CamberCurveGen2/generateIdealCamberCurves_simplified.m: C.rear_ARB_roll_stiffness
assert(abs(cfg.studies.camber2_generateIdealCamberCurves_simplified.C_weight_dist-(0.54)) < 1e-10); % Full Car Models/CamberCurveGen2/generateIdealCamberCurves_simplified.m: C.weight_dist
assert(abs(cfg.studies.camber2_generateIdealCamberCurves_simplified.C_mass-((395 + 150))) < 1e-10); % Full Car Models/CamberCurveGen2/generateIdealCamberCurves_simplified.m: C.mass
assert(abs(tireCfg.studies.camber2_test.friction_scaling_factor-(1.05*0.55)) < 1e-10); % Full Car Models/CamberCurveGen2/test.m: friction_scaling_factor
assert(abs(tireCfg.studies.camber2_test.p_i-(12)) < 1e-10); % Full Car Models/CamberCurveGen2/test.m: p_i
assert(abs(cfg.studies.calcWheelForcesAndDisplacements_pitch.front_anti_lift-(0)) < 1e-10); % calcWheelForcesAndDisplacements_pitch.m: front_anti_lift
assert(abs(cfg.studies.calcWheelForcesAndDisplacements_pitch.front_anti_dive-(.545)) < 1e-10); % calcWheelForcesAndDisplacements_pitch.m: front_anti_dive
assert(abs(cfg.studies.calcWheelForcesAndDisplacements_pitch.rear_anti_lift-(1.02)) < 1e-10); % calcWheelForcesAndDisplacements_pitch.m: rear_anti_lift
assert(abs(cfg.studies.calcWheelForcesAndDisplacements_pitch.rear_anti_squat-(.457)) < 1e-10); % calcWheelForcesAndDisplacements_pitch.m: rear_anti_squat
assert(abs(cfg.studies.roll.mass-(370 * 0.45359)) < 1e-10); % roll.m: mass
assert(abs(cfg.studies.roll.r_g-(0.85)) < 1e-10); % roll.m: r_g
assert(abs(cfg.studies.roll.k_ride_f-(23177.80762)) < 1e-10); % roll.m: k_ride_f
assert(abs(cfg.studies.roll.k_ride_r-(26147.44256)) < 1e-10); % roll.m: k_ride_r
assert(abs(cfg.studies.roll.cLa-(3.969)) < 1e-10); % roll.m: cLa
assert(abs(cfg.studies.roll.CoP-(0.418)) < 1e-10); % roll.m: CoP
assert(abs(cfg.studies.roll.rho-(1.2)) < 1e-10); % roll.m: rho
assert(abs(cfg.studies.understeer_bode.L-(62 * 0.0254)) < 1e-10); % understeer_bode.m: L
assert(abs(cfg.studies.understeer_bode.m-(535 * 0.453592)) < 1e-10); % understeer_bode.m: m
assert(abs(cfg.studies.understeer_bode.Iz-(83.28)) < 1e-10); % understeer_bode.m: Iz
assert(abs(tireCfg.studies.understeer_bode.Csum-(1.6e5)) < 1e-10); % understeer_bode.m: Csum
assert(abs(tireCfg.studies.understeer_budget.CAF-(136*4.448)) < 1e-10); % understeer_budget.m: CAF
assert(abs(tireCfg.studies.understeer_budget.CCF-(0.313)) < 1e-10); % understeer_budget.m: CCF
assert(abs(tireCfg.studies.understeer_budget.NAF-(50.082)) < 1e-10); % understeer_budget.m: NAF
assert(abs(tireCfg.studies.understeer_budget.CGF-(47.657)) < 1e-10); % understeer_budget.m: CGF
assert(abs(tireCfg.studies.understeer_budget.NGF-(4.766)) < 1e-10); % understeer_budget.m: NGF
assert(abs(cfg.studies.understeer_budget.m-(168)) < 1e-10); % understeer_budget.m: m
assert(abs(cfg.studies.understeer_budget.L-(62 * 25.4)) < 1e-10); % understeer_budget.m: L
assert(abs(cfg.studies.understeer_budget.weight_dist-(0.51)) < 1e-10); % understeer_budget.m: weight_dist
assert(abs(cfg.studies.understeer_budget.WUF-(17.07)) < 1e-10); % understeer_budget.m: WUF
assert(abs(cfg.studies.understeer_budget.EF-(0.001)) < 1e-10); % understeer_budget.m: EF
assert(abs(cfg.studies.understeer_budget.GF-(0.789)) < 1e-10); % understeer_budget.m: GF
assert(abs(cfg.studies.understeer_budget.KF-(2288)) < 1e-10); % understeer_budget.m: KF
assert(abs(cfg.studies.understeer_budget.EYF-(0.191)) < 1e-10); % understeer_budget.m: EYF
assert(abs(cfg.studies.understeer_budget.GYF-(0.291)) < 1e-10); % understeer_budget.m: GYF
assert(abs(cfg.studies.understeer_budget.ENF-(0.774)) < 1e-10); % understeer_budget.m: ENF
assert(abs(cfg.studies.understeer_budget.GNF-(0.022)) < 1e-10); % understeer_budget.m: GNF
assert(abs(cfg.studies.understeer_budget.ROLLpg-(4.47)) < 1e-10); % understeer_budget.m: ROLLpg
assert(abs(tireCfg.studies.understeer_budget.CAR-(136*4.448)) < 1e-10); % understeer_budget.m: CAR
assert(abs(tireCfg.studies.understeer_budget.CCR-(0.323)) < 1e-10); % understeer_budget.m: CCR
assert(abs(tireCfg.studies.understeer_budget.NAR-(42.897)) < 1e-10); % understeer_budget.m: NAR
assert(abs(tireCfg.studies.understeer_budget.CGR-(43.392)) < 1e-10); % understeer_budget.m: CGR
assert(abs(tireCfg.studies.understeer_budget.NGR-(4.339)) < 1e-10); % understeer_budget.m: NGR
assert(abs(cfg.studies.understeer_budget.FYRperg-(8430.9)) < 1e-10); % understeer_budget.m: FYRperg
assert(abs(cfg.studies.understeer_budget.NFperg-(-403.16)) < 1e-10); % understeer_budget.m: NFperg
assert(abs(cfg.studies.understeer_budget.NRperg-(-317.93)) < 1e-10); % understeer_budget.m: NRperg
assert(abs(cfg.studies.understeer_budget.WR-(885.000)) < 1e-10); % understeer_budget.m: WR
assert(abs(cfg.studies.understeer_budget.WUR-(99.000)) < 1e-10); % understeer_budget.m: WUR
assert(abs(cfg.studies.understeer_budget.ER-(0.055)) < 1e-10); % understeer_budget.m: ER
assert(abs(cfg.studies.understeer_budget.GR-(-0.720)) < 1e-10); % understeer_budget.m: GR
assert(abs(cfg.studies.understeer_budget.KR-(367)) < 1e-10); % understeer_budget.m: KR
assert(abs(cfg.studies.understeer_budget.EYR-(0.045)) < 1e-10); % understeer_budget.m: EYR
assert(abs(cfg.studies.understeer_budget.GYR-(-0.328)) < 1e-10); % understeer_budget.m: GYR
assert(abs(cfg.studies.understeer_budget.ENR-(-0.113)) < 1e-10); % understeer_budget.m: ENR
assert(abs(cfg.studies.understeer_budget.GNR-(0.023)) < 1e-10); % understeer_budget.m: GNR
assert(abs(cfg.studies.Bicycle_Model_Simulink_Models_BicycleModelSim.a-(1.1)) < 1e-10); % Bicycle Model/Simulink Models/BicycleModelSim.m: a
assert(abs(cfg.studies.Bicycle_Model_Simulink_Models_BicycleModelSim.b-(1.6)) < 1e-10); % Bicycle Model/Simulink Models/BicycleModelSim.m: b
assert(abs(tireCfg.studies.understeer_bode.K_deg_per_g-(0.18)) < 1e-10); % understeer_bode.m: K_deg_per_g
assert(abs(cfg.studies.understeer_bode.weight_dist-0.512) < 1e-10);
assert(abs(tireCfg.studies.Bicycle_Model_BicycleModel.cornering_force_scale-(0.55)) < 1e-10);
assert(abs(tireCfg.studies.Bicycle_Model_BicycleModel.cornering_load_mass_lb-(535)) < 1e-10);
assert(abs(tireCfg.studies.Bicycle_Model_BicycleModel.cornering_pressure_psi-(12)) < 1e-10);
assert(abs(tireCfg.studies.Bicycle_Model_toe_compliance.cornering_force_scale-(0.55)) < 1e-10);
assert(abs(tireCfg.studies.Bicycle_Model_toe_compliance.cornering_load_mass_lb-(535)) < 1e-10);
assert(abs(tireCfg.studies.Bicycle_Model_toe_compliance.cornering_pressure_psi-(12)) < 1e-10);
count = 185;
end
