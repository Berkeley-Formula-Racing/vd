function [Fx_tot, Fx_front, Fx_rear, M_x_L, M_x_R, kappa_vec] = BicycleCamberEvaluation(F_z_F, F_z_R, gamma_L, gamma_R, tire,state)
    %negative pitch angle = nose down
    front_tire = tire;
    rear_tire = tire;
    front_tire.gamma = gamma_L;
    rear_tire.gamma = gamma_R;
    %F_z_F = 4.44822*F_z_F;
    %F_z_R = 4.44822*F_z_R;
    
    
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
        Fx_front_temp = F_x(front_tire, alpha, kappa, F_z_F);
        Fx_rear_temp = F_x(rear_tire, alpha, kappa, F_z_R);

        %check if kappa vals give better acceleration
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
    if state == "brake"
        Fx_tot= -Fx_tot;

    else
        Fx_tot=Fx_tot - Fx_front;
    end
    Fx_tot = 2*Fx_tot;

    
    %still have to add overturning moment
    M_x_L = 0;
    M_x_R = 0;
end