function camber = Camber_Evaluation(long_vel, yaw_rate, steer_angle_1, steer_angle_2, ...
                                    static_camber_f, static_camber_r, ccVal_f, ccVal_r)

    % --- constants & tunables ---
    roll_grad_deg_per_g = 0.68;
    g0 = 9.81;
    rear_out_deg_per_deg = 0.58;
    rear_in_deg_per_deg  = -0.592;

    % --- direction & roll ---
    dir     = sign(yaw_rate);
    ay_g    = (long_vel * yaw_rate) / g0;
    rolldeg = ay_g * roll_grad_deg_per_g;

    % --- load compact params once ---
    persistent betaL betaR
    if isempty(betaL)
        S     = load('camber_models_fast.mat','betaL','betaR');
        betaL = S.betaL;  betaR = S.betaR;
    end

    % --- fast polynomial eval: b0 + b1*r + b2*s + b3*r^2 + b4*s^2 + b5*r*s ---
    % note the steer sign mapping you had:
    steerL = -dir * steer_angle_1;
    steerR =  dir * steer_angle_2;
    r = rolldeg;

    cam_FL = betaL(1) + betaL(2)*r + betaL(3)*steerL + betaL(4)*r.^2 + betaL(5)*steerL.^2 + betaL(6)*(r.*steerL);
    cam_FR = betaR(1) + betaR(2)*r + betaR(3)*steerR + betaR(4)*r.^2 + betaR(5)*steerR.^2 + betaR(6)*(r.*steerR);

    % --- compliance (inner +, outer -) ---
    % abs(ay_g), NOT ay_g. sign_comp already carries the turn direction
    % through dir, and ay_g carries it a second time -- multiplying both
    % cancels the direction out entirely, so the compliance camber came out
    % identical in left and right turns. That applies it by SIDE (left wheels
    % always negative) rather than by inner/outer, which is not what
    % compliance does: it comes from lateral force at the contact patch, and
    % that reverses when the car turns the other way.
    %
    % Measured before the fix, at cc = 1.5 deg/g: right turn gave
    % [-0.917 +0.917 -0.917 +0.917] and the mirrored left turn gave exactly
    % the same, instead of the negated set.
    %
    % Same class of error as the static camber sign fixed below -- a signed
    % quantity multiplied by a sign that already encodes the same thing.
    sign_comp  = [-dir;  dir; -dir;  dir];
    c_fromcomp = [ccVal_f; ccVal_f; ccVal_r; ccVal_r] .* ay_g .* sign_comp;

    % --- rear roll->camber ---
    % SAE convention (see yaw moment equation in Car.equations): tires 1/3 are
    % left, positive yaw rate is a right turn, so dir >= 0 puts the LEFT wheels
    % on the outside. This matches the sign_comp assignment above.
    if dir >= 0
        cam_RL_roll = rear_out_deg_per_deg * r;
        cam_RR_roll = rear_in_deg_per_deg  * r;
    else
        cam_RL_roll = rear_in_deg_per_deg  * r;
        cam_RR_roll = rear_out_deg_per_deg * r;
    end

    % --- assemble ---
    % Static camber is MIRRORED left/right. The tire model takes a signed
    % inclination in a single-tire frame, so a symmetric setup (e.g. -1 deg
    % on both sides, both leaning inboard) must enter as +g on the left and
    % -g on the right -- otherwise the two camber thrusts add instead of
    % cancelling and the car generates a side force and yaw moment while
    % driving dead straight. That made the straight-line constraint set in
    % constraint1 infeasible, which is why max_long_accel never converged.
    %
    % Sign check: in a right turn (dir>0) the left wheel is outer; roll makes
    % it lean outboard, and the compliance term gives it a NEGATIVE increment
    % (sign_comp(1) = -dir). So "leaning outboard" is negative on the left,
    % hence static negative camber -- which leans the left wheel INBOARD --
    % must be positive there. Every dynamic term above is already antisymmetric.
    staticSign = [-1; 1; -1; 1];
    camber        = zeros(4,1);
    camber(1)     = staticSign(1)*static_camber_f + cam_FL      + c_fromcomp(1);
    camber(2)     = staticSign(2)*static_camber_f + cam_FR      + c_fromcomp(2);
    camber(3)     = staticSign(3)*static_camber_r + cam_RL_roll + c_fromcomp(3);
    camber(4)     = staticSign(4)*static_camber_r + cam_RR_roll + c_fromcomp(4);
end
