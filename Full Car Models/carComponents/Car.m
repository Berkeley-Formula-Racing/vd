classdef Car
    % 3 DOF Model
    % equations adapted from Casanova Appendix

    properties
        M %mass
        W_b %wheelbase
        l_f %dist from cg to front axle
        l_r %dist from cg to rear axle
        t_f %trackwidth, front
        t_r %trackwidth, rear
        h_rr %roll center at rear
        h_rf %roll center at front
        h_rc %roll center at cg (approx)
        R %wheel radius
        h_g %cg height
        R_sf %roll stiffness in front
        I_zz %polar moment of inertia, z axis
        static_gamma_f
        static_gamma_r %static camber
        static_r_toe
        camber_compliance_f
        camber_compliance_r
        % Values are supplied by carConfig. The defaults retain the previous
        % model for old saved cars that predate this configuration struct.
        camberKinematics = struct('roll_gradient_deg_per_g',0.68, ...
            'rear_roll_camber_outer_deg_per_deg',0.58, ...
            'rear_roll_camber_inner_deg_per_deg',-0.592, ...
            'ride_camber_front_deg_per_in',0, ...
            'ride_camber_rear_deg_per_in',0)
        aero
        powertrain
        tire
        ackermann
        g = 9.81;

        % Per-axle grip scaling, applied to the tire forces in tireForce.
        % The model carries ONE Tire2 for all four corners, so this is how
        % front and rear grip are varied independently -- see gripSweep.
        % 1 = the tire model as fitted; 0.95 = 5% less grip on that axle.
        %
        % Defaults matter: a Car saved before these existed loads with these
        % values, so old .mat sweep caches keep working unchanged.
        gripScaleF = 1;
        gripScaleR = 1;

        % Rotating inertia, kg-m^2. The longitudinal equation of motion used
        % to divide the net tyre force by M alone, so the wheels, brake discs,
        % crank and driveline all counted as zero -- the car accelerated as if
        % nothing had to be spun up. Against a measured 75 m run the model was
        % 5.8% quick and, tellingly, completely insensitive to tyre grip
        % (4.5222 to 4.5241 s across a whole grip sweep), which is the
        % signature of a missing inertia rather than a missing force.
        %
        % Default 0 reproduces the old behaviour exactly, so a Car sitting in
        % an old sweep cache loads and solves as it always did. carConfig sets
        % real values, so anything built from carConfig gets them.
        I_wheel = 0;      % PER WHEEL: tyre, rim, hub, brake disc
        I_driveline = 0;  % crank + clutch + primary, referred to the CRANK

        % Rolling resistance coefficient. A rearward force Crr*(M*g + downforce)
        % opposing travel, added to the longitudinal balance in equations()
        % alongside aero drag. Absent before this -- the model had no rolling
        % resistance at all, which at 8 m/s is 92% of the aero drag it did have.
        % Default 0 reproduces the old behaviour, so old sweep caches are
        % unaffected; carConfig sets the real value.
        Crr = 0;

        % Ride-height-dependent aero configuration. carConfig stores spring
        % rates at the damper and converts them to wheel rates using MR^2.
        % Disabled by default so saved cars retain static-aero behaviour.
        rideHeightAero = struct('enabled',false);

        Iyy
        Ixx   % needs to be about roll center
        k     % spring rate (assumed same over all tires) N/m
        k_c   % Chassis TS (Nm/deg)
        rs_total % Total Roll Stiffness (Nm/deg)
        k_tf % tire stiffness (N/m)
        k_tr % tire stiffness (N/m)
        k_rf  % front arb roll stiffness (Nm/rad)
        k_rr  % rear arb roll stiffness (Nm/rad)
        c_compression % damper curves ([in/s, lbf])
        c_rebound % damper curves ([in/s, lbf])
        MR_F % front motion ratio curve ([in,MR])
        MR_R % rear motion ratio curve ([in,MR])
        TSmpc % mpc timestep
        TSdyn % dynamics timestep

        % Decoupled Suspension Parameters
        k_f_r % front roll spring stiffness (N/m)
        k_f_b % front bounce spring stiffness (N/m)
        k_r_r % rear roll spring stiffness (N/m)
        k_r_b % rear bounce spring stiffness (N/m)

        Jm %engine polar moi
        Jw %wheel polar moi

        ggPoints %g-g diagram points for car instance
        ss_info % information (x vector including normal loads, etc) for pure cornering
        accel_info % information (x vector including normal loads, etc) when accelerating
        decel_info % information (x vector including normal loads, etc) when decelerating
        longAccelLookup %maxLongAccel = f(latAccel,velocity)
        longDecelLookup %maxLongDecel = f(latAccel,velocity)
        comp
    end

    methods
        function obj = Car(mass,wheelbase,weight_dist,track_width,wheel_radius,cg_height,...
                roll_center_height_front,roll_center_height_rear,R_sf,I_zz,static_gamma_f,static_gamma_r,camber_compliance_f,camber_compliance_r,aero,powertrain,tire,ackermann,static_r_toe,...
                gripScaleF,gripScaleR,I_wheel,I_driveline,Crr,rideHeightAero,camberKinematics)
            % gripScaleF/R and the two inertias are optional and default to
            % the values declared above, so a call written before they existed
            % still builds a car that behaves as it did then
            obj.M = mass;
            obj.W_b = wheelbase;
            obj.l_f = wheelbase*weight_dist; % distance from cg to front
            obj.l_r = wheelbase*(1-weight_dist); % distance from cg to rear
            obj.t_f = track_width;
            obj.t_r = track_width;
            obj.h_rr = roll_center_height_rear;
            obj.h_rf = roll_center_height_front;
            obj.h_rc = (obj.h_rf+obj.h_rr)/2; % approximation of roll center height at cg
            obj.R = wheel_radius;
            obj.h_g = cg_height;
            obj.R_sf = R_sf;
            obj.I_zz = I_zz;
            obj.static_gamma_f = static_gamma_f;
            obj.static_gamma_r = static_gamma_r;
            obj.aero = aero;
            obj.powertrain = powertrain;
            obj.tire = tire;
            obj.ackermann = ackermann;
            obj.camber_compliance_f = camber_compliance_f;
            obj.camber_compliance_r = camber_compliance_r;
            obj.static_r_toe = static_r_toe;
            if nargin >= 20 && ~isempty(gripScaleF),   obj.gripScaleF   = gripScaleF;   end
            if nargin >= 21 && ~isempty(gripScaleR),   obj.gripScaleR   = gripScaleR;   end
            if nargin >= 22 && ~isempty(I_wheel),      obj.I_wheel      = I_wheel;      end
            if nargin >= 23 && ~isempty(I_driveline),  obj.I_driveline  = I_driveline;  end
            if nargin >= 24 && ~isempty(Crr),          obj.Crr          = Crr;          end
            if nargin >= 25 && ~isempty(rideHeightAero)
                obj.rideHeightAero = rideHeightAero;
            end
            if nargin >= 26 && ~isempty(camberKinematics)
                obj.camberKinematics = camberKinematics;
            end
        end

        function m = rotatingMass(obj,current_gear)
            % Equivalent translating mass, in kg, of everything that has to be
            % angularly accelerated along with the car.
            %
            % Each wheel contributes I/R^2. The driveline sits UPSTREAM of the
            % gearbox, so its inertia is referred through the square of the
            % total reduction and is therefore strongly gear-dependent -- on
            % this car the same crank is worth several times more in first
            % than in fifth. That is exactly why the effect is largest where
            % the acceleration event spends its time.
            %
            % Only the inertia is added here. The WEIGHT of these parts is
            % already in M, so load transfer and the static axle loads must
            % keep using M on its own.
            n = obj.powertrain.drivetrain_reduction(current_gear);
            m = (4*obj.I_wheel + obj.I_driveline*n^2)/obj.R^2;
        end

        function [engine_rpm,beta,lat_accel,long_accel,yaw_accel,wheel_accel,omega,current_gear,...
                Fzvirtual,Fz,alpha,T,Fy, gamma, Fx, ssInfo] = equations(obj,P)

            % inputs: vehicle parameters
            % outputs: vehicle accelerations and other properties

            % state and control matrix
            steer_angle = P(1);
            throttle = P(2); % -1 for full braking, 1 for full throttle
            long_vel = P(3); % m/s
            lat_vel = P(4); % m/s
            yaw_rate = P(5); % equal to long_vel/radius (v/r)
            kappa = P(6:9);
            % note: 1 = front left tire, 2 = front right tire
            %       3 = rear left tire, 4 = rear right tire

            % Powertrain
            omega = zeros(1,4);
            omega(1) = (kappa(1)+1)/obj.R*(long_vel+yaw_rate*obj.t_f/2);
            omega(2) = (kappa(2)+1)/obj.R*(long_vel-yaw_rate*obj.t_f/2);
            omega(3) = (kappa(3)+1)/obj.R*(long_vel+yaw_rate*obj.t_r/2);
            omega(4) = (kappa(4)+1)/obj.R*(long_vel-yaw_rate*obj.t_r/2);

            [engine_rpm,current_gear] = obj.powertrain.engine_rpm(omega(3),omega(4),long_vel);
            [T_1,T_2,T_3,T_4] = obj.powertrain.wheel_torques(engine_rpm, omega(3), omega(4), throttle, current_gear, long_vel);
            T = [T_1,T_2,T_3,T_4];


            % Tire Slips
            beta = atan(lat_vel/long_vel)*180/pi; % vehicle slip angle in deg

            a_poly = [1.000170221169974e-05,-6.101610470854148e-05,0.009416280591935,0.999321611021116];

            if steer_angle >= 0
                steer_angle_1 = steer_angle;
                steer_angle_2 = abs(obj.ackermann)*polyval(a_poly, steer_angle)*steer_angle;%(0.0094*steer_angle + 0.9993)*steer_angle;
            else
                steer_angle_2 = steer_angle;
                steer_angle_1 = -(abs(obj.ackermann)*polyval(a_poly, -steer_angle)*(-steer_angle));%-(0.0094*(-steer_angle) + 0.9993)*(-steer_angle);
            end

            if obj.ackermann < 0
                [steer_angle_1, steer_angle_2] = deal(steer_angle_2, steer_angle_1);
            end


            %steer_angle_1 = steer_angle; % could be modified for ackermann steering
            %steer_angle_2 = steer_angle;

            %disp([steer_angle_1 steer_angle_2]);

            % only assemble the diagnostics struct when a caller asks for it:
            % this runs inside every fmincon function evaluation
            needRideCamber = obj.hasRideCamber();
            if nargout >= 16
                [Fz,Fzvirtual,downforce,drag,wheelRideHeights,ssInfo] = ssForces(obj,long_vel,yaw_rate,T,(1/2)*(steer_angle_1+steer_angle_2)*pi/180);
            elseif needRideCamber
                [Fz,Fzvirtual,downforce,drag,wheelRideHeights] = ssForces(obj,long_vel,yaw_rate,T,(1/2)*(steer_angle_1+steer_angle_2)*pi/180);
            else
                [Fz,Fzvirtual,downforce,drag] = ssForces(obj,long_vel,yaw_rate,T,(1/2)*(steer_angle_1+steer_angle_2)*pi/180);
                wheelRideHeights = [];
            end

            % slip angles (small angle assumption)
            alpha(1) = -steer_angle_1+(lat_vel+obj.l_f*yaw_rate)/(long_vel+yaw_rate*obj.t_f/2)*180/pi; %deg
            alpha(2) = -steer_angle_2+(lat_vel+obj.l_f*yaw_rate)/(long_vel-yaw_rate*obj.t_f/2)*180/pi; %deg
            % Rear static toe is MIRRORED left to right. Both wheels used to
            % get -static_r_toe, which is toe-in on one side and toe-out on
            % the other: the car then makes a side force and yaw moment while
            % driving straight, exactly as the static camber term did before
            % it was given staticSign. Dormant while static_r_toe = 0, but at
            % 0.5 deg it put 2.5 m/s^2 of lateral asymmetry into the model.
            % Sign order matches the camber convention: 3 = RL, 4 = RR.
            alpha(3) = -obj.static_r_toe + (lat_vel-obj.l_r*yaw_rate)/(long_vel+yaw_rate*obj.t_r/2)*180/pi;
            alpha(4) = +obj.static_r_toe + (lat_vel-obj.l_r*yaw_rate)/(long_vel-yaw_rate*obj.t_r/2)*180/pi;


            %gamma = [obj.static_gamma obj.static_gamma obj.static_gamma obj.static_gamma];
            gamma = Camber_Evaluation(long_vel, yaw_rate, steer_angle_1, steer_angle_2, ...
                obj.static_gamma_f, obj.static_gamma_r, obj.camber_compliance_f, ...
                obj.camber_compliance_r,obj.camberKinematics, ...
                obj.rideCompression(wheelRideHeights)).';

            %disp(gamma);
            %disp("------------");

            % Tire Forces
            %steer_angle = steer_angle_1*pi/180;
            %disp(gamma);
            [Fx,Fy,Fxw] = obj.tireForce(steer_angle_1,steer_angle_2,alpha,kappa,Fz, gamma);

            % Equations of Motion
            lat_accel = sum(Fy)*(1/obj.M)-yaw_rate*long_vel;
            % Rolling resistance: a rearward force scaling with the total
            % vertical load, downforce included, opposing travel. long_vel > 0
            % everywhere in the g-g and every event, so the sign is fixed --
            % it subtracts from tractive force and ADDS to braking, both
            % correct. This is the single place longitudinal force is summed,
            % so putting it here reaches the g-g solvers, the accel-event
            % lookup and the lap solver alike. Zero Crr is a no-op.
            F_rr = obj.Crr * (obj.M*obj.g + downforce);

            % M + rotatingMass, not M: the tyre force has to accelerate the
            % spinning parts as well as the car. Zero inertias reduce this to
            % the old expression exactly.
            long_accel = (sum(Fx)-drag-F_rr)*(1/(obj.M+obj.rotatingMass(current_gear)))+yaw_rate*lat_vel;
            yaw_accel = ((Fx(1)-Fx(2))*obj.t_f/2+(Fx(3)-Fx(4))*obj.t_r/2+(Fy(1)+Fy(2))*obj.l_f-(Fy(3)+Fy(4))*obj.l_r)*(1/obj.I_zz);
            %yaw_accel = ((Fy(1)+Fy(2))*obj.l_f-(Fy(3)+Fy(4))*obj.l_r)*(1/obj.I_zz);

            % neglects wheel rotational dynamics: for justification see Koutrik p.16
            wheel_accel(1) = (T(1)-Fxw(1)*obj.R);
            wheel_accel(2) = (T(2)-Fxw(2)*obj.R);
            wheel_accel(3) = (T(3)-Fx(3)*obj.R);
            wheel_accel(4) = (T(4)-Fx(4)*obj.R);
        end

        function m = metrics(obj,P)
            % Full named record of one solved operating point.
            % P is the 9-element state/control vector (same as equations).
            % Returns a scalar struct; struct arrays of these concatenate
            % straight into a table via struct2table.

            [engine_rpm,beta,lat_accel,long_accel,yaw_accel,wheel_accel,omega,current_gear,...
                Fzvirtual,Fz,alpha,T,Fy,gamma,Fx,ss] = obj.equations(P);

            % --- inputs / state ---
            m.steer_angle = P(1);
            m.throttle    = P(2);
            m.vCar        = P(3);
            m.lat_vel     = P(4);
            m.yaw_rate    = P(5);
            m.beta        = beta;

            % --- accelerations (g), gLat/gLong are the steady-state values ---
            % lat_accel/long_accel out of equations are the residuals that the
            % optimiser drives to zero, so the achieved accelerations are the
            % kinematic ones: a_y = v*r, a_x = the residual plus v*r*beta terms
            m.gLat   = (P(3)*P(5))/obj.g;
            m.gLong  = long_accel/obj.g;
            m.gMag   = hypot(m.gLat,m.gLong);
            % g-g polar angle: 0 deg = pure acceleration, +90 = pure cornering,
            % 180 = pure braking
            m.theta  = atan2d(m.gLat,m.gLong);
            m.lat_accel_residual = lat_accel;
            m.yaw_accel_residual = yaw_accel;

            % --- powertrain ---
            m.engine_rpm   = engine_rpm;
            m.current_gear = current_gear;

            % --- aero ---
            m.downforce = ss.downforce;
            m.drag      = ss.drag;
            m.ClA       = ss.ClA;
            m.CdA       = ss.CdA;
            m.CoP       = ss.CoP;             % front downforce fraction
            m.LoD       = ss.downforce/max(ss.drag,eps);
            m.aero_downforce_front_N = ss.aero_downforce_front_N;
            m.aero_downforce_rear_N  = ss.aero_downforce_rear_N;
            m.front_ride_height_in = ss.front_ride_height_in;
            m.rear_ride_height_in  = ss.rear_ride_height_in;
            m.front_ride_height_offset_in = ss.front_ride_height_offset_in;
            m.rear_ride_height_offset_in  = ss.rear_ride_height_offset_in;
            m.ride_height_FL_in = ss.wheel_ride_height_in(1);
            m.ride_height_FR_in = ss.wheel_ride_height_in(2);
            m.ride_height_RL_in = ss.wheel_ride_height_in(3);
            m.ride_height_RR_in = ss.wheel_ride_height_in(4);
            m.aero_iterations = ss.aero_iterations;
            m.aero_residual_in = ss.aero_residual_in;
            m.aero_outside_map = ss.aero_outside_map;

            % --- load distribution ---
            m.Fz_front_axle = ss.Fz_front_axle;
            m.Fz_rear_axle  = ss.Fz_rear_axle;
            m.front_Fz_frac = ss.Fz_front_axle/(ss.Fz_front_axle+ss.Fz_rear_axle);
            m.LLT_front     = ss.lat_load_transfer_front;
            m.LLT_rear      = ss.lat_load_transfer_rear;
            m.LLTD          = ss.LLTD;
            m.long_load_transfer = ss.long_load_transfer;
            m.min_Fz        = min(Fzvirtual);   % negative => wheel lift

            % --- per corner (1 FL, 2 FR, 3 RL, 4 RR) ---
            for i = 1:4
                s = num2str(i);
                m.(['Fz_' s])    = Fz(i);
                m.(['Fy_' s])    = Fy(i);
                m.(['Fx_' s])    = Fx(i);
                m.(['alpha_' s]) = alpha(i);
                m.(['gamma_' s]) = gamma(i);
                m.(['kappa_' s]) = P(5+i);
                m.(['T_' s])     = T(i);
                m.(['omega_' s]) = omega(i);
                m.(['wheel_accel_residual_' s]) = wheel_accel(i);
            end

            % --- axle summaries ---
            m.Fy_front = Fy(1)+Fy(2);
            m.Fy_rear  = Fy(3)+Fy(4);
            % front share of total lateral force: >l_r/W_b implies the front
            % is working harder than its static load share (understeer margin)
            m.Fy_front_frac = m.Fy_front/max(abs(m.Fy_front+m.Fy_rear),eps);
        end

        function [Fx,Fy,F_xw] = tireForce(obj,steer_angle_1,steer_angle_2,alpha,kappa,Fz,gamma)
            %radians
            % Per-axle grip scaling. Applied to the tire FORCES, not to Fz, so
            % it is a friction multiplier rather than a load change -- same
            % thing tire.friction_scaling_factor does globally, split by axle.
            % Scales Fx as well as Fy: a friction change is not lateral-only,
            % and F_xw feeds the wheel torque balance, so scaling one without
            % the other would leave wheel_accel inconsistent.
            % isempty guard covers Car objects saved before these properties
            % existed, which load with [] rather than the class default.
            sF = obj.gripScaleF; if isempty(sF), sF = 1; end
            sR = obj.gripScaleR; if isempty(sR), sR = 1; end

            % forces in tire frame of reference
            F_xw1 = sF*obj.tire.F_x(alpha(1),kappa(1),Fz(1),gamma(1));
            F_yw1 = sF*obj.tire.F_y(alpha(1),kappa(1),Fz(1),gamma(1));
            F_xw2 = sF*obj.tire.F_x(alpha(2),kappa(2),Fz(2),gamma(2));
            F_yw2 = sF*obj.tire.F_y(alpha(2),kappa(2),Fz(2),gamma(2));
            F_xw = [F_xw1; F_xw2];

            % forces in vehicle frame of reference
            F_x1 = F_xw1*cosd(steer_angle_1)-F_yw1*sind(steer_angle_1);
            F_y1 = F_xw1*sind(steer_angle_1)+F_yw1*cosd(steer_angle_1);
            F_x2 = F_xw2*cosd(steer_angle_2)-F_yw2*sind(steer_angle_2);
            F_y2 = F_xw2*sind(steer_angle_2)+F_yw2*cosd(steer_angle_2);

            F_x3 = sR*obj.tire.F_x(alpha(3),kappa(3),Fz(3),gamma(3));
            F_y3 = sR*obj.tire.F_y(alpha(3),kappa(3),Fz(3),gamma(3));
            F_x4 = sR*obj.tire.F_x(alpha(4),kappa(4),Fz(4),gamma(4));
            F_y4 = sR*obj.tire.F_y(alpha(4),kappa(4),Fz(4),gamma(4));

            Fx = [F_x1; F_x2; F_x3; F_x4];
            Fy = [F_y1; F_y2; F_y3; F_y4];
        end

        function [forces, Gr] = calcForces(obj,x,u,forces)
            % takes spring-damper forces, adds powertrain, aero, tireXY
            throttle = u(2); %[-1,1] max braking to max throttle
            longVel = x(3); %m/s

            % powertrain
            omega = [x(8); x(10); x(12); x(14)];
            [engineRPM,currentGear] = obj.powertrain.engine_rpm(omega(3),omega(4),longVel);
            [T1,T2,T3,T4] = obj.powertrain.wheel_torques(engineRPM, omega(3), omega(4), throttle, currentGear);
            T = [T1,T2,T3,T4];
            Gr = obj.powertrain.drivetrain_reduction(currentGear);
            forces.T = T;

            % aero
            F_lift = [0 0 obj.aero.lift(longVel)];
            R_lift = [obj.aero.D_f*obj.W_b-obj.l_r 0 0];
            F_drag = [-obj.aero.drag(longVel) 0 0];
            R_drag = [0 0 0];
            lift = [F_lift R_lift];
            drag = [F_drag R_drag];
            forces.F = [forces.F; lift; drag];
        end

        function forces = calcTireForces(obj,x,u,forces)
            steerAngle = u(1); %steering angle, radians
            yawRate = x(2); %rad/s
            longVel = x(3);
            latVel = x(4); %m/s
            Fz = forces.Ftires(:,3);

            % slip angles
            alphaR = [steerAngle-atan((latVel+obj.l_f*yawRate)/abs(longVel-yawRate*obj.t_f/2));
                steerAngle-atan((latVel+obj.l_f*yawRate)/abs(longVel+yawRate*obj.t_f/2));
                -atan((latVel-obj.l_r*yawRate)/abs(longVel-yawRate*obj.t_r/2));
                -atan((latVel-obj.l_r*yawRate)/abs(longVel+yawRate*obj.t_r/2))];
            alphaR = -alphaR;
            alphaD = rad2deg(alphaR);

            % slip ratios
            k1 = (obj.R*x(8)/(x(3)-x(2)*obj.t_f/2))-1;
            k2 = (obj.R*x(10)/(x(3)+x(2)*obj.t_f/2))-1;
            k3 = (obj.R*x(12)/(x(3)-x(2)*obj.t_f/2))-1;
            k4 = (obj.R*x(14)/(x(3)+x(2)*obj.t_f/2))-1;
            kappa = [k1; k2; k3; k4];

            % calculate tire forces
            [Fx,Fy,Fxw] = tireForce(obj,steerAngle,alphaD,kappa,Fz,gamma);
            Rtire = [obj.l_f -obj.t_f/2 0;   %tire 1
                     obj.l_f obj.t_f/2 0;   %tire 2
                     -obj.l_r -obj.t_f/2 0;   %tire 3
                     -obj.l_r obj.t_f/2 0]; %tire 4
            forces.alpha = alphaR;

            Ftires = [Fx Fy Fz Rtire];
            forces.Ftires = Ftires;
            forces.Fxw = Fxw;
            forces.Fx = Fx;
        end

        function [xdot, forces] = dynamics(obj,x,forces,Gr)
            % initial vehicle states (vector of 14 values)
            % 1: yaw angle 2: yaw rate 3: long velocity 4: lat velocity
            % 5: x position of cg 6: y position of cg
            % 7:  FL angular position 8:  FL angular velocity
            % 9:  FR angular position 10: FR angular velocity
            % 11: RL angular position 12: RL angular velocity
            % 13: RR angular position 14: RR angular velocity
            % u(1): steering input u(2): throttle

            longVel = x(3); %m/s
            latVel = x(4); %m/s
            psi = x(1);
            psid = x(2);
            beta = rad2deg(atan(latVel/longVel)); % vehicle slip angle in deg

            Fx = forces.Fx;
            T = forces.T;
            Fxw = forces.Fxw;
            Fapplied = forces.F(:,1:3); % applied Fxyz: car frame
            xF = forces.F(:,4:6); % position vectors Xxyz: car frame
            psiMoments = 0;

            %add up all applied moments, using given position vectors
            for i = 1:size(Fapplied,1)
                psiMoments = psiMoments + det([xF(i,1:2);Fapplied(i,1:2)]);
            end
            Ftires = forces.Ftires(:,1:3);
            rTires = forces.Ftires(:,4:6);
            for i = 1:size(Ftires,1)
                psiMoments = psiMoments + det([rTires(i,1:2);Ftires(i,1:2)]);
            end

            %total matrix of forces in vehicle axes (e1, e2)
            allForces = [Fapplied(:,1:2); Ftires(:,1:2)];

            % total acceleration vector
            sumA = sum(allForces,1)/obj.M;

            xdot = zeros(14,1);

            %yaw velocity
            xdot(1) = psid;
            xdot(2) = (1/obj.I_zz)*psiMoments;

            %long accel, lat accel. Vehicle coordinates
            xdot(3) = sumA(1)+psid*latVel;
            if longVel <=0
                xdot(3) = max(0,sumA(1)+psid*latVel);
            end

            %xdot(3) = 0; % for pure cornering studies ONLY

            xdot(4) = sumA(2)-psid*longVel;
            %X velocity, Y velocity. Global coordinates
            xdot(5) = longVel*cos(-psi)-latVel*sin(-psi);
            xdot(6) = longVel*sin(-psi)+latVel*cos(-psi);

            %tires: angular velocity, acceleration, 1-4
            xdot(7) = x(8);
            xdot(8) = (T(1) - Fxw(1)*obj.R)/obj.Jw;
            xdot(9) = x(10);
            xdot(10) = (T(2) - Fxw(2)*obj.R)/obj.Jw;
            denom = (obj.Jw^2+2*obj.Jw*obj.Jm*(Gr/2)^2);
            xdot(11) = x(12);
            xdot(12) = ((T(3)-Fx(3)*obj.R)*(obj.Jw+obj.Jm*(Gr/2)^2) - (T(4)-Fx(4)*obj.R)*obj.Jm*(Gr/2)^2)*(1/denom);
            xdot(13) = x(14);
            xdot(14) = ((T(4)-Fx(4)*obj.R)*(obj.Jw+obj.Jm*(Gr/2)^2) - (T(3)-Fx(3)*obj.R)*obj.Jm*(Gr/2)^2)*(1/denom);
        end

        function [Fz_f,Fz_r,aeroInfo] = FzForces(obj,longVel,T)
            % Coupled F/R ride-height aero. The static branch intentionally
            % keeps legacy behaviour for saved cars and acceleration aero.
            if obj.hasRideHeightAero()
                aeroInfo = obj.solveRideHeightAero(longVel,T);
            else
                coeff = obj.aero.coefficients(0,0);
                aeroInfo = obj.aeroLoadsAtHeights(longVel,T,coeff,NaN,NaN,0,0);
            end
            Fz_f = aeroInfo.Fz_front_axle;
            Fz_r = aeroInfo.Fz_rear_axle;
        end

        function [Fz,Fzvirtual,downforce,drag,wheelRideHeights,ssInfo] = ssForces(obj,longVel,yawRate,T,steer_angle)
            % Fifth output is the numeric corner ride-height vector needed by
            % ride-camber. Sixth exposes diagnostic load-transfer terms.

            [Fz_front,Fz_rear,aeroInfo] = FzForces(obj,longVel,T);
            downforce = aeroInfo.downforce;
            drag = aeroInfo.drag;


            lat_load_transfer_front = (yawRate*longVel*obj.M)/obj.t_f*((obj.l_r*obj.h_rf)/obj.W_b+...
                obj.R_sf*(obj.h_g-obj.h_rc));
            lat_load_transfer_rear = (yawRate*longVel*obj.M)/obj.t_r*((obj.l_f*obj.h_rr)/obj.W_b+...
                (1-obj.R_sf)*(obj.h_g-obj.h_rc));

%             LLTD_caster = obj.R_sf-0.06*steer_angle/(25*pi/180);
%             lat_load_transfer_front = (yawRate*longVel*obj.M)/obj.t_f*((obj.l_r*obj.h_rf)/obj.W_b+...
%                 LLTD_caster*(obj.h_g-obj.h_rc));
%             lat_load_transfer_rear = (yawRate*longVel*obj.M)/obj.t_r*((obj.l_r*obj.h_rr)/obj.W_b+...
%                 (1-LLTD_caster)*(obj.h_g-obj.h_rc));
%
            % wheel load constraint method from Kelly
            Fzvirtual = zeros(1,4);
            Fzvirtual(1) = 0.5*Fz_front+lat_load_transfer_front;
            Fzvirtual(2) = 0.5*Fz_front-lat_load_transfer_front;
            Fzvirtual(3) = 0.5*Fz_rear+lat_load_transfer_rear;
            Fzvirtual(4) = 0.5*Fz_rear-lat_load_transfer_rear;

            % smooth approximation of max function
            epsilon = 10;
            Fz = (Fzvirtual + sqrt(Fzvirtual.^2 + epsilon))./2;

            wheelRideHeights = [];
            if nargout > 4
                wheelRideHeights = obj.wheelRideHeights(Fzvirtual,aeroInfo);
            end
            % Keep struct construction out of the fmincon hot path, which
            % requests only Fz/Fzvirtual/downforce/drag unless ride camber is on.
            if nargout > 5
                ssInfo.Fz_front_axle = Fz_front;
                ssInfo.Fz_rear_axle = Fz_rear;
                ssInfo.lat_load_transfer_front = lat_load_transfer_front;
                ssInfo.lat_load_transfer_rear = lat_load_transfer_rear;
                % LLTD: front share of total lateral load transfer
                totalLLT = lat_load_transfer_front + lat_load_transfer_rear;
                if abs(totalLLT) < eps
                    ssInfo.LLTD = NaN;
                else
                    ssInfo.LLTD = lat_load_transfer_front/totalLLT;
                end
                ssInfo.downforce = downforce;
                ssInfo.drag = drag;
                ssInfo.long_load_transfer = aeroInfo.long_load_transfer;
                ssInfo.ClA = aeroInfo.cla;
                ssInfo.CdA = aeroInfo.cda;
                ssInfo.CoP = aeroInfo.D_f;
                ssInfo.aero_downforce_front_N = aeroInfo.downforce_front;
                ssInfo.aero_downforce_rear_N = aeroInfo.downforce_rear;
                ssInfo.front_ride_height_in = aeroInfo.frontRideHeightIn;
                ssInfo.rear_ride_height_in = aeroInfo.rearRideHeightIn;
                ssInfo.front_ride_height_offset_in = aeroInfo.frontOffsetIn;
                ssInfo.rear_ride_height_offset_in = aeroInfo.rearOffsetIn;
                ssInfo.wheel_ride_height_in = wheelRideHeights;
                ssInfo.aero_iterations = aeroInfo.iterations;
                ssInfo.aero_residual_in = aeroInfo.residualIn;
                ssInfo.aero_outside_map = aeroInfo.outsideMap;
            end
        end

        function tf = hasRideHeightAero(obj)
            cfg = obj.rideHeightAero;
            tf = isstruct(cfg) && isfield(cfg,'enabled') && cfg.enabled && ...
                isfield(cfg,'wheel_rate_front_Npm') && ...
                isfield(cfg,'wheel_rate_rear_Npm') && ...
                obj.aero.hasMap();
        end

        function tf = hasRideCamber(obj)
            cfg = obj.camberKinematics;
            tf = obj.hasRideHeightAero() && isstruct(cfg) && ...
                isfield(cfg,'ride_camber_front_deg_per_in') && ...
                isfield(cfg,'ride_camber_rear_deg_per_in') && ...
                any([cfg.ride_camber_front_deg_per_in,cfg.ride_camber_rear_deg_per_in] ~= 0);
        end

        function compressionIn = rideCompression(obj,wheelRideHeights)
            if isempty(wheelRideHeights) || any(~isfinite(wheelRideHeights)) || ...
                    ~obj.hasRideCamber()
                compressionIn = [0;0];
                return
            end
            cfg = obj.rideHeightAero;
            staticHeights = [cfg.static_front_ride_height_in; ...
                cfg.static_front_ride_height_in;cfg.static_rear_ride_height_in; ...
                cfg.static_rear_ride_height_in];
            compressionIn = staticHeights-wheelRideHeights(:);
        end

        function info = solveRideHeightAero(obj,longVel,T)
            % Solve r = height - height(load(height)) = 0 with a damped
            % two-variable Newton method. This avoids the unconverged one-step
            % pitch update the old model used.
            cfg = obj.rideHeightAero;
            heights = [cfg.static_front_ride_height_in; ...
                       cfg.static_rear_ride_height_in];
            [res,info] = obj.aeroRideResidual(heights,longVel,T);
            toleranceIn = 1e-7;
            maxIterations = 12;
            iterations = 0;

            while max(abs(res)) > toleranceIn && iterations < maxIterations
                iterations = iterations + 1;
                h = 1e-4;
                J = zeros(2,2);
                for j = 1:2
                    plus = heights;  plus(j) = plus(j) + h;
                    minus = heights; minus(j) = minus(j) - h;
                    rPlus = obj.aeroRideResidual(plus,longVel,T);
                    rMinus = obj.aeroRideResidual(minus,longVel,T);
                    J(:,j) = (rPlus-rMinus)/(2*h);
                end
                if all(isfinite(J),'all') && rcond(J) > 1e-10
                    step = -J\res;
                else
                    step = -0.5*res;
                end
                step = max(min(step,0.25),-0.25); % inches per iteration

                oldNorm = norm(res,inf);
                accepted = false;
                for lineSearch = 1:6
                    candidate = heights + step;
                    [candidateRes,candidateInfo] = obj.aeroRideResidual(candidate,longVel,T);
                    if norm(candidateRes,inf) < oldNorm
                        heights = candidate;
                        res = candidateRes;
                        info = candidateInfo;
                        accepted = true;
                        break
                    end
                    step = 0.5*step;
                end
                if ~accepted
                    heights = heights - 0.5*res;
                    [res,info] = obj.aeroRideResidual(heights,longVel,T);
                end
            end
            info.iterations = iterations;
            info.residualIn = max(abs(res));
            info.converged = info.residualIn <= toleranceIn;
        end

        function [res,info] = aeroRideResidual(obj,heights,longVel,T)
            cfg = obj.rideHeightAero;
            coeff = obj.aero.coefficients( ...
                heights(1)-cfg.map_reference_front_ride_height_in, ...
                heights(2)-cfg.map_reference_rear_ride_height_in);
            info = obj.aeroLoadsAtHeights(longVel,T,coeff,heights(1),heights(2),0,0);
            targetFront = cfg.static_front_ride_height_in - ...
                ((info.Fz_front_axle-info.static_front_load_N)/2) / ...
                cfg.wheel_rate_front_Npm / 0.0254;
            targetRear = cfg.static_rear_ride_height_in - ...
                ((info.Fz_rear_axle-info.static_rear_load_N)/2) / ...
                cfg.wheel_rate_rear_Npm / 0.0254;
            res = heights - [targetFront;targetRear];
        end

        function info = aeroLoadsAtHeights(obj,longVel,T,coeff,frontHeightIn,rearHeightIn,iterations,residualIn)
            staticFront = obj.M*obj.g*obj.l_r/obj.W_b;
            staticRear  = obj.M*obj.g*obj.l_f/obj.W_b;
            dynamicPressure = obj.aero.rho/2*longVel^2;
            downforce = dynamicPressure*coeff.cla;
            drag = dynamicPressure*coeff.cda;
            longLoadTransfer = (sum(T)/obj.R-drag)*(obj.h_g/obj.W_b);

            info.static_front_load_N = staticFront;
            info.static_rear_load_N = staticRear;
            info.Fz_front_axle = staticFront + downforce*coeff.D_f - longLoadTransfer;
            info.Fz_rear_axle  = staticRear + downforce*coeff.D_r + longLoadTransfer;
            info.downforce = downforce;
            info.drag = drag;
            info.downforce_front = downforce*coeff.D_f;
            info.downforce_rear = downforce*coeff.D_r;
            info.long_load_transfer = longLoadTransfer;
            info.cla = coeff.cla;
            info.cda = coeff.cda;
            info.D_f = coeff.D_f;
            info.D_r = coeff.D_r;
            info.frontRideHeightIn = frontHeightIn;
            info.rearRideHeightIn = rearHeightIn;
            info.frontOffsetIn = coeff.frontOffsetIn;
            info.rearOffsetIn = coeff.rearOffsetIn;
            info.outsideMap = coeff.outsideMap;
            info.iterations = iterations;
            info.residualIn = residualIn;
            info.converged = false;
        end

        function heights = wheelRideHeights(obj,Fzvirtual,aeroInfo)
            if ~obj.hasRideHeightAero()
                heights = nan(1,4);
                return
            end
            cfg = obj.rideHeightAero;
            staticLoads = [aeroInfo.static_front_load_N/2, ...
                aeroInfo.static_front_load_N/2,aeroInfo.static_rear_load_N/2, ...
                aeroInfo.static_rear_load_N/2];
            wheelRates = [cfg.wheel_rate_front_Npm,cfg.wheel_rate_front_Npm, ...
                cfg.wheel_rate_rear_Npm,cfg.wheel_rate_rear_Npm];
            staticHeights = [cfg.static_front_ride_height_in, ...
                cfg.static_front_ride_height_in,cfg.static_rear_ride_height_in, ...
                cfg.static_rear_ride_height_in];
            % Fzvirtual already contains the TOTAL corner load change from
            % static, including the axle-average aero/load-transfer term.
            % Starting from axleHeights would subtract that average deflection
            % a second time and make the wheel-height mean disagree with the
            % coupled F/R height used by the aeromap.
            heights = staticHeights - (Fzvirtual-staticLoads)./wheelRates/0.0254;
        end

        function plotGG(car)
            figure(123);clf;
            scatter3(car.ggPoints(:,1),car.ggPoints(:,2),car.ggPoints(:,3),'+')
            xlabel('Long Accel (m/s^2)');
            ylabel('Lat Accel (m/s^2)');
            zlabel('Velocity (m/s)');
        end

        % These functions are used to set constraints for fmincon
        % input P: state and control vector containing:
        %   steer angle,throttle position,longitudinal velocity,
        %   lateral velocity,yaw rate,wheel rotational speeds

        % output c: limits vehicle slip angle to less than 20 degrees (stability purposes)
        %   also limits engine rpm to below 13000
        %   also limits wheel loads to positive values (no wheel lift)
        % output ceq: constrains certain accelerations to 0 to satisfy
        %   steady-state conditions

        function [c,ceq] = constraint1(obj,P)
            % no lateral acceleration constraint
            % used for optimizing longitudinal acceleration/braking
            % note: callers that want symmetric rear slip (kappa_4 = kappa_3)
            % must impose it through Aeq, not here -- the objective and the
            % constraints have to evaluate the same state vector

            [engine_rpm,beta,lat_accel,long_accel,yaw_accel,wheel_accel,omega,current_gear,...
                Fzvirtual,Fz,alpha,T] = obj.equations(P);
            c = [engine_rpm-13000,abs(beta)-20,-Fzvirtual(1:4)];
            ceq = [lat_accel,yaw_accel,wheel_accel(1:4)];
        end

        function [c,ceq] = constraint2(obj,P,long_accel_value)
            % longitudinal acceleration constrained to equal long_accel_value
            % used for optimizing lateral force for given longitudinal acceleration

            [engine_rpm,beta,lat_accel,long_accel,yaw_accel,wheel_accel,omega,current_gear,...
                Fzvirtual,Fz,alpha,T]...
                = obj.equations(P);
            c = [engine_rpm-13000,abs(beta)-20,-Fzvirtual(1:4)];
            ceq = [lat_accel,long_accel-long_accel_value,yaw_accel,wheel_accel(1:4)];
        end

        function [c,ceq] = constraint3(obj,P,radius)
            % longitudinal acceleration constrained to equal zero
            % velocity divided by yaw rate constrained to equal inputted radius
            % used for solving skidpad (optimizing velocity for zero longitudinal acceleration

            [engine_rpm,beta,lat_accel,long_accel,yaw_accel,wheel_accel,omega,current_gear,...
                Fzvirtual,Fz,alpha,T] = obj.equations(P);
            c = [engine_rpm-13000,abs(beta)-20,-Fzvirtual(1:4)];
            ceq = [P(3)/(P(5))-radius,lat_accel,long_accel,yaw_accel,wheel_accel(1:4)];
        end

        function [c,ceq] = constraint4(obj,P,lat_accel_value)
            % lateral acceleration constrained to equal lat_accel_value
            % used for optimizing longitudinal acceleration for given lateral acceleration

            [engine_rpm,beta,lat_accel,long_accel,yaw_accel,wheel_accel,omega,current_gear,...
                Fzvirtual,Fz,alpha,T]...
                = obj.equations(P);
            c = [engine_rpm-13000,abs(beta)-20,-Fzvirtual(1:4)];
            ceq = [P(3)*P(5)-lat_accel_value,lat_accel,yaw_accel,wheel_accel(1:4)];
        end

        function [c,ceq] = constraint5(obj,P,radius)
            % velocity divided by yaw rate constrained to equal inputted radius
            % used for calculating max velocity the car can corner at for given radius

            [engine_rpm,beta,lat_accel,long_accel,yaw_accel,wheel_accel,omega,current_gear,...
                Fzvirtual,Fz,alpha,T]...
                = obj.equations(P);
            c = [engine_rpm-13000,abs(beta)-20,-Fzvirtual(1:4)];
            ceq = [P(3)/(P(5))-radius,lat_accel,yaw_accel,wheel_accel(1:4)];
        end

        function [c,ceq] = constraint6(obj,P)
            % no longitudinal acceleration constraint
            % used for optimizing lateral acceleration

           [engine_rpm,beta,lat_accel,long_accel,yaw_accel,wheel_accel,omega,current_gear,...
                Fzvirtual,Fz,alpha,T]...
                = obj.equations(P);
            c = [engine_rpm-13000,abs(beta)-20,-Fzvirtual(1:4)];
            ceq = [lat_accel,yaw_accel,wheel_accel(1:4)];
        end

        function [c,ceq] = constraint7(obj,P)
            % longitudinal acceleration constrained to 0
            % no yaw accel constraint
            % used to determine terminal under/oversteer

            [engine_rpm,beta,lat_accel,long_accel,yaw_accel,wheel_accel,omega,current_gear,...
                Fzvirtual,Fz,alpha,T]...
                = obj.equations(P);
            c = [engine_rpm-13000,abs(beta)-20,-Fzvirtual(1:4)];
            ceq = [lat_accel,long_accel,wheel_accel(1:4)];
        end

        function [c,ceq] = constraint8(obj,P,radius,long_vel_value)
            % velocity divided by yaw rate constrained to equal inputted radius
            % velocity constrained to equal long_vel_value
            % longitudinal acceleration constrained to 0
            % used for constant radius test

           [engine_rpm,beta,lat_accel,long_accel,yaw_accel,wheel_accel,omega,current_gear,...
                Fzvirtual,Fz,alpha,T]...
                = obj.equations(P);
            c = [engine_rpm-13000,abs(beta)-20,-Fzvirtual(1:4)];
            ceq = [P(3)/(P(5))-radius,P(3)-long_vel_value,lat_accel,...
                long_accel, yaw_accel,wheel_accel(1:4)];
        end

        function [c,ceq] = constraint9(obj,P,lat_accel_value, radius)
            % lateral acceleration constrained to equal lat_accel_value and
            % v^2/r
            % used for optimizing longitudinal acceleration for given
            % lateral acceleration AND radius

            [engine_rpm,beta,lat_accel,long_accel,yaw_accel,wheel_accel,omega,current_gear,...
                Fzvirtual,Fz,alpha,T]...
                = obj.equations(P);
            c = [engine_rpm-13000,abs(beta)-20,-Fzvirtual(1:4)];
            ceq = [P(3)*P(5)-lat_accel_value,P(3)/P(5)-radius, yaw_accel, wheel_accel(1:4)];
        end

        % objective function
        function out = long_accel(obj,P)
            % used for optimizing longitudinal acceleration
            % evaluates P exactly as given: fmincon's objective must see the
            % same state as its nonlinear constraints

            [engine_rpm,beta,lat_accel,long_accel,yaw_accel,wheel_accel,omega,current_gear,...
                Fzvirtual,Fz,alpha,T]...
                = obj.equations(P);
            out = long_accel;
        end

        function dxdt = PhasePlaneODE(obj,x,steer_angle)

%             % state and control matrix
%             steer_angle = P(1);
%             throttle = P(2); % -1 for full braking, 1 for full throttle
%             long_vel = P(3); % m/s
%             lat_vel = P(4); % m/s
%             yaw_rate = P(5); % equal to long_vel/radius (v/r)
%             kappa = P(6:9);

            P(1) = steer_angle;
            P(2) = 0;
            P(3) = x(1); % m/s
            P(4) = x(2); % m/s
            P(5) = x(3); % equal to long_vel/radius (v/r)
            P(6) = 0;
            P(7) = 0;
            P(8) = 0;
            P(9) = 0;

            [~,~,lat_accel,long_accel,yaw_accel,~,~,~,...
                ~,~,~,~,~] = equations(obj,P);

            lat_accel = lat_accel+P(3)*P(5);

            long_accel = 0;
            dxdt = [long_accel; lat_accel; yaw_accel];
        end

        % maximum possible car velocity
        function out = max_vel(obj)
            out = obj.powertrain.redline*pi/30*obj.R/...
                obj.powertrain.drivetrain_reduction(numel(obj.powertrain.gears))-0.001;
        end

        function printState(obj,x,xdot)
            fprintf("1. yaw angle: %0.2f\n",x(1));
            fprintf("2. yaw rate : %0.2f\n",x(2));
            fprintf("3. long velo: %0.2f\n",x(3));
            fprintf("  long accel: %0.2f\n",xdot(3));
            fprintf("4. lat velo : %0.2f\n",x(4));
            fprintf("   lat accel: %0.2f\n",xdot(4));
            fprintf("5.  Xcg     : %0.2f\n",x(5));
            fprintf("6.  Ycg     : %0.2f\n",x(6));
            fprintf("7.  FL theta: %0.2f\n",rad2deg(x(7)));
            fprintf("8.  FL w    : %0.2f\n",rad2deg(x(8)));
            fprintf("    FL a    : %0.2f\n",rad2deg(xdot(8)));
            fprintf("9.  FR theta: %0.2f\n",rad2deg(x(9)));
            fprintf("10. FR w    : %0.2f\n",rad2deg(x(10)));
            fprintf("    FR a    : %0.2f\n",rad2deg(xdot(10)));
            fprintf("11. RL theta: %0.2f\n",rad2deg(x(11)));
            fprintf("12. RL w    : %0.2f\n",rad2deg(x(12)));
            fprintf("    RL a    : %0.2f\n",rad2deg(xdot(12)));
            fprintf("13. RR theta: %0.2f\n",rad2deg(x(13)));
            fprintf("14. RR w    : %0.2f\n",rad2deg(x(14)));
            fprintf("    RR a    : %0.2f\n",rad2deg(xdot(14)));
            fprintf("\n");
        end
    end

end
