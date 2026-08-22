function initialize_tire()
    % tire parameters
    tireParams = struct();
    gamma = 0; % camber angle
    p_i = 12; % pressure
    % these parameters are non-iterable
    load('Fx_combined_parameters_run38_30.mat'); % F_x combined magic formula parameters
    Fx_parameters = cell2mat(Xbestcell);
    load('Lapsim_Fy_combined_parameters_1965run15.mat'); % F_y combined magic formula parameters
    Fy_parameters = cell2mat(Xbestcell);
    friction_scaling_factor = 1.05*0.55; % scales tire forces to account for test/road surface difference

    tire = BDITire2(gamma,p_i,Fx_parameters,Fy_parameters,friction_scaling_factor);
end
