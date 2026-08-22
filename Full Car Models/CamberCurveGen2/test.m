load('Fx_combined_parameters_run38_30.mat'); % F_x combined magic formula parameters
Fx_parameters = cell2mat(Xbestcell);
load('Fy_combined_parameters_16x10x7.5_LC0.mat'); % F_y combined magic formula parameters
Fy_parameters = cell2mat(Xbestcell);
friction_scaling_factor = 1.05*0.55;
p_i = 12; % tire pressure

%initialize tire object, used for camber evaluation
tire = BDITire2(3,p_i,Fx_parameters,Fy_parameters,friction_scaling_factor);

F_z_vec = 0:1:250;
F_y_max_vec = 0*F_z_vec;
for i = 1:numel(F_z_vec)
    alpha_vec = -20:20;
    F_y_vec = 0*alpha_vec;
    for j = 1:numel(alpha_vec)
        F_y_vec(j) = F_y(tire, alpha_vec(j), 0, F_z_vec(i));
    end
    %plot(alpha_vec, F_y_vec, 'displayName', string(F_z_vec(i)));
    hold on;
    F_y_max_vec(i) = max(F_y_vec);
end

%figure;

plot(F_z_vec,F_y_max_vec./F_z_vec);
hold on;




