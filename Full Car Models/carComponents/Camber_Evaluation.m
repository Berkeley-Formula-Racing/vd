function camber = Camber_Evaluation(long_vel, yaw_rate, steer_angle_1, steer_angle_2, ...
                                    static_camber_f, static_camber_r, ccVal_f, ccVal_r)

    roll_grad_deg_per_g = 0.68;
    g0 = 9.81;
    rear_out_deg_per_deg = 0.58;
    rear_in_deg_per_deg  = -0.592;

    dir     = sign(yaw_rate);
    ay_g    = (long_vel * yaw_rate) / g0;
    rolldeg = ay_g * roll_grad_deg_per_g;

    persistent betaL betaR
    if isempty(betaL)
        S     = load('camber_models_fast.mat','betaL','betaR');
        betaL = S.betaL;  betaR = S.betaR;
    end

    steerL = -dir * steer_angle_1;
    steerR =  dir * steer_angle_2;
    r = rolldeg;

    cam_FL = betaL(1) + betaL(2)*r + betaL(3)*steerL + betaL(4)*r.^2 + betaL(5)*steerL.^2 + betaL(6)*(r.*steerL);
    cam_FR = betaR(1) + betaR(2)*r + betaR(3)*steerR + betaR(4)*r.^2 + betaR(5)*steerR.^2 + betaR(6)*(r.*steerR);

    sign_comp  = [-dir;  dir; -dir;  dir];
    c_fromcomp = [ccVal_f; ccVal_f; ccVal_r; ccVal_r] .* ay_g .* sign_comp;

    if dir >= 0
        cam_RL_roll = rear_in_deg_per_deg  * r;
        cam_RR_roll = rear_out_deg_per_deg * r;
    else
        cam_RL_roll = rear_out_deg_per_deg * r;
        cam_RR_roll = rear_in_deg_per_deg  * r;
    end

    camber        = zeros(4,1);
    camber(1)     = static_camber_f + cam_FL      + c_fromcomp(1);
    camber(2)     = static_camber_f + cam_FR      + c_fromcomp(2);
    camber(3)     = static_camber_r + cam_RL_roll + c_fromcomp(3);
    camber(4)     = static_camber_r + cam_RR_roll + c_fromcomp(4);
end
