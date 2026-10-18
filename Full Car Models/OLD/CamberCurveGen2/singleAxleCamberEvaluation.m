function [F_y_tot, F_y_L, F_y_R, M_x_L, M_x_R, alpha_val] = singleAxleCamberEvaluation(F_z_L, F_z_R, gamma_L, gamma_R, tire)
    %for negative camber, gamma_L -> + & gamma_R -> -
    left_tire = tire;
    right_tire = tire;
    left_tire.gamma = gamma_L;
    right_tire.gamma = gamma_R;
   
    
    % sweep through slip angles to find maximum lateral force for the
    % entire axle
    F_y_tot = 0;
    F_y_L = 0;
    F_y_R = 0;
    alpha_val = 0;
    for alpha = -25:1:25
        kappa = 0;
        % turning right - for sign convention | F_y -> + | alpha -> - 
        F_y_L_temp = F_y(left_tire, alpha, kappa, F_z_L);
        F_y_R_temp = F_y(right_tire, alpha, kappa, F_z_R);
        if F_y_L_temp + F_y_R_temp > F_y_tot
            F_y_L = F_y_L_temp;
            F_y_R = F_y_R_temp;
            F_y_tot = F_y_L + F_y_R;
            alpha_val = alpha;
        end
    end
    
    %still have to add overturning moment
    M_x_L = 0;
    M_x_R = 0;
end