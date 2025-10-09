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
        static_gamma %static camber
        camber_compliance
        aero
        powertrain
        tire
        ackermann
        g = 9.81;
        
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
                roll_center_height_front,roll_center_height_rear,R_sf,I_zz,static_gamma,camber_compliance,aero,powertrain,tire,ackermann);
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
            obj.static_gamma = static_gamma;
            obj.aero = aero;
            obj.powertrain = powertrain;
            obj.tire = tire;
            obj.ackermann = ackermann;
            obj.camber_compliance = camber_compliance;
        end
        
        function [engine_rpm,beta,lat_accel,long_accel,yaw_accel,wheel_accel,omega,current_gear,...
                Fzvirtual,Fz,alpha,T,Fy, gamma] = equations(obj,P)           
            
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
            omega(3) = (kappa(3)+1)/obj.R*(long_vel+yaw_rate*obj.t_f/2);
            omega(4) = (kappa(4)+1)/obj.R*(long_vel-yaw_rate*obj.t_f/2);
                        
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

            [Fz, Fzvirtual] = ssForces(obj,long_vel,yaw_rate,T,(1/2)*(steer_angle_1+steer_angle_2)*pi/180);
            
            % slip angles (small angle assumption)
            alpha(1) = -steer_angle_1+(lat_vel+obj.l_f*yaw_rate)/(long_vel+yaw_rate*obj.t_f/2)*180/pi; %deg
            alpha(2) = -steer_angle_2+(lat_vel+obj.l_f*yaw_rate)/(long_vel-yaw_rate*obj.t_f/2)*180/pi; %deg
            alpha(3) = (lat_vel-obj.l_r*yaw_rate)/(long_vel+yaw_rate*obj.t_r/2)*180/pi;
            alpha(4) = (lat_vel-obj.l_r*yaw_rate)/(long_vel-yaw_rate*obj.t_r/2)*180/pi;
         
            %tire camber using tire.Fy with static camber 
            % (kind of stupid, very inneficient, but it will work for now
            % maybe?)
            %maybe theres a more efficient way to approximate Fy with slip
            %angle? We have lat vel, long vel, yaw rate, slip angles
            fyApprox = zeros(1,4);
            
            fyApprox = (Fz./(sum(Fz))) * (yaw_rate) * (long_vel) * obj.M;
            %disp(fyApprox);
            %{
            fyApprox(1) = obj.tire.F_y(alpha(1),kappa(1),Fz(1),-obj.static_gamma); 

            fyApprox(2) = obj.tire.F_y(alpha(2),kappa(2),Fz(2),obj.static_gamma);
            fyApprox(3) = obj.tire.F_y(alpha(3),kappa(3),Fz(3),-obj.static_gamma);
            fyApprox(4) = obj.tire.F_y(alpha(4),kappa(4),Fz(4),obj.static_gamma);
            %}

            %gamma = [obj.static_gamma obj.static_gamma obj.static_gamma obj.static_gamma];
            gamma = Camber_Evaluation(long_vel, yaw_rate, steer_angle_1, steer_angle_2, obj.static_gamma, fyApprox, obj.camber_compliance).';
            %disp(gamma);
            %disp("------------");
            
            % Tire Forces
            %steer_angle = steer_angle_1*pi/180;
            %disp(gamma);
            [Fx,Fy,Fxw] = obj.tireForce(steer_angle_1,steer_angle_2,alpha,kappa,Fz, gamma);
                        
            % Equations of Motion
            lat_accel = sum(Fy)*(1/obj.M)-yaw_rate*long_vel;
            long_accel = (sum(Fx)-obj.aero.drag(long_vel))*(1/(obj.M-10))+yaw_rate*lat_vel;
            yaw_accel = ((Fx(1)-Fx(2))*obj.t_f/2+(Fx(3)-Fx(4))*obj.t_r/2+(Fy(1)+Fy(2))*obj.l_f-(Fy(3)+Fy(4))*obj.l_r)*(1/obj.I_zz);
            %yaw_accel = ((Fy(1)+Fy(2))*obj.l_f-(Fy(3)+Fy(4))*obj.l_r)*(1/obj.I_zz);
    
            % neglects wheel rotational dynamics: for justification see Koutrik p.16
            wheel_accel(1) = (T(1)-Fxw(1)*obj.R);
            wheel_accel(2) = (T(2)-Fxw(2)*obj.R);
            wheel_accel(3) = (T(3)-Fx(3)*obj.R);
            wheel_accel(4) = (T(4)-Fx(4)*obj.R); 
        end
        
        function [Fx,Fy,F_xw] = tireForce(obj,steer_angle_1,steer_angle_2,alpha,kappa,Fz,gamma)
            %radians
            % forces in tire frame of reference
            F_xw1 = obj.tire.F_x(alpha(1),kappa(1),Fz(1),-gamma(1)); 
            F_yw1 = obj.tire.F_y(alpha(1),kappa(1),Fz(1),-gamma(1));
            F_xw2 = obj.tire.F_x(alpha(2),kappa(2),Fz(2),gamma(2));
            F_yw2 = obj.tire.F_y(alpha(2),kappa(2),Fz(2),gamma(2));
            F_xw = [F_xw1; F_xw2];

            % forces in vehicle frame of reference
            F_x1 = F_xw1*cosd(steer_angle_1)-F_yw1*sind(steer_angle_1);
            F_y1 = F_xw1*sind(steer_angle_1)+F_yw1*cosd(steer_angle_1);
            F_x2 = F_xw2*cosd(steer_angle_2)-F_yw2*sind(steer_angle_2);
            F_y2 = F_xw2*sind(steer_angle_2)+F_yw2*cosd(steer_angle_2);
            
            F_x3 = obj.tire.F_x(alpha(3),kappa(3),Fz(3),-gamma(3));
            F_y3 = obj.tire.F_y(alpha(3),kappa(3),Fz(3),-gamma(3));
            F_x4 = obj.tire.F_x(alpha(4),kappa(4),Fz(4),gamma(4));
            F_y4 = obj.tire.F_y(alpha(4),kappa(4),Fz(4),gamma(4));

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
        
        function [Fz, Fzvirtual] = ssForces(obj,longVel,yawRate,T,steer_angle)
            Fz_front_static = (obj.M*9.81*obj.l_r+obj.aero.lift(longVel)*obj.aero.D_f)/obj.W_b;
            Fz_rear_static = (obj.M*9.81*obj.l_f+obj.aero.lift(longVel)*obj.aero.D_r)/obj.W_b;
            
            Fy_front_approx = longVel*yawRate*obj.M*(obj.l_r/obj.W_b);
            Fx_fromFy = -Fy_front_approx*sin(steer_angle);
            %(F_x1+F_x2+F_x3+F_x4)*h_g/W_b approximated (neglecting wheel dynamics) since longitudinal forces are unknown
            % Fx_fromFy term added to avoid overestimation of long load
            %   transfer. Approximates Fx component from front tire lateral
            %   forces
            %long_load_transfer = (sum(T)/obj.R+Fx_fromFy)*(obj.h_g/obj.W_b);            
            long_load_transfer = (sum(T)/obj.R)*(obj.h_g/obj.W_b);
            
            
            lat_load_transfer_front = (yawRate*longVel*obj.M)/obj.t_f*((obj.l_r*obj.h_rf)/obj.W_b+...
                obj.R_sf*(obj.h_g-obj.h_rc));
            lat_load_transfer_rear = (yawRate*longVel*obj.M)/obj.t_r*((obj.l_r*obj.h_rr)/obj.W_b+...
                (1-obj.R_sf)*(obj.h_g-obj.h_rc));
            
%             LLTD_caster = obj.R_sf-0.06*steer_angle/(25*pi/180);
%             lat_load_transfer_front = (yawRate*longVel*obj.M)/obj.t_f*((obj.l_r*obj.h_rf)/obj.W_b+...
%                 LLTD_caster*(obj.h_g-obj.h_rc));
%             lat_load_transfer_rear = (yawRate*longVel*obj.M)/obj.t_r*((obj.l_r*obj.h_rr)/obj.W_b+...
%                 (1-LLTD_caster)*(obj.h_g-obj.h_rc));
%             
            % wheel load constraint method from Kelly
            Fzvirtual = zeros(1,4);
            Fzvirtual(1) = 0.5*Fz_front_static-0.5*long_load_transfer+lat_load_transfer_front;
            Fzvirtual(2) = 0.5*Fz_front_static-0.5*long_load_transfer-lat_load_transfer_front;
            Fzvirtual(3) = 0.5*Fz_rear_static+0.5*long_load_transfer+lat_load_transfer_rear;
            Fzvirtual(4) = 0.5*Fz_rear_static+0.5*long_load_transfer-lat_load_transfer_rear;

            % smooth approximation of max function
            epsilon = 10; 
            Fz = (Fzvirtual + sqrt(Fzvirtual.^2 + epsilon))./2;
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
            P(9) = P(8);
            
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
            c = [engine_rpm-13000,abs(beta)-20,-Fzvirtual(1:2)];
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
            c = [engine_rpm-13000,abs(beta)-20,-Fzvirtual(1:2)];
            ceq = [P(3)*P(5)-lat_accel_value,P(3)/P(5)-radius, yaw_accel, wheel_accel(1:4)];
        end
        
        % objective function
        function out = long_accel(obj,P)
            % used for optimizing longitudinal acceleration            
            P(9) = P(8);
            
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

