function [Fx_tot, Fx_front, Fx_rear, M_x_L, M_x_R, kappa_vec] = AccelerationCamberEvaluation(F_z_F, F_z_R, gamma_L, gamma_R, tire,state)
    %negative pitch angle = nose down
    left_tire = tire;
    right_tire = tire;
    left_tire.gamma = gamma_L;
    right_tire.gamma = gamma_R;
    
    
    if state == 'brake'

    end
    % sweep through slip angles to find maximum lateral force for the
    % entire axle
    Fx_tot = 0;
    Fx_front = 0;
    Fx_rear = 0;
    kappa_val_f = 0;
    kappa_val_r = 0;
    for kappa = 0:.05:2
        alpha = 0;
        % acceleration forward = + | Fx -> +
        Fx_front_temp = F_x(left_tire, alpha, kappa, F_z_F);
        Fx_rear_temp = F_x(right_tire, alpha, kappa, F_z_R);
        if Fx_front_temp > Fx_front
            Fx_front = Fx_front_temp;
            Fx_tot = Fx_front + Fx_rear;
            kappa_val_f = kappa;
        elseif Fx_rear_temp > Fx_rear
            Fx_rear = Fx_rear_temp;
            Fx_tot = Fx_front + Fx_rear;
            kappa_val_r = kappa;
        end
    end
    kappa_vec = [kappa_val_f,kappa_val_r];
    
    %still have to add overturning moment
    M_x_L = 0;
    M_x_R = 0;
end