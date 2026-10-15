addpath(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'Full Car Models','config'));
sharedStudyTire = vdTireConfig([],[],false);
sharedFitData = load(sharedStudyTire.studies.camber1_test.fit_1,'Xbestcell');
Xbestcell = sharedFitData.Xbestcell; % F_x combined magic formula parameters
Fx_parameters = cell2mat(Xbestcell);
sharedFitData = load(sharedStudyTire.studies.camber1_test.fit_2,'Xbestcell');
Xbestcell = sharedFitData.Xbestcell; % F_y combined magic formula parameters
Fy_parameters = cell2mat(Xbestcell);
friction_scaling_factor = sharedStudyTire.studies.camber1_test.friction_scaling_factor;
p_i = sharedStudyTire.studies.camber1_test.p_i; % tire pressure

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




