classdef car_pm
    properties
        % basic
        m double
        cda double
        cla double
        cop double    
        weight_r double   
        wheel_radius double
        rho double = 1.2
        g   double = 9.81
        c_roll double = 0.025

        tireMu tire_pm
        ptrain powertrain_pm

        % braking
        max_brake_torque double = 850
        tire_limited_brake logical = true
    end

    methods
        function obj = car_pm(params, tireMu, ptrain)
            arguments
                params struct
                tireMu tire_pm
                ptrain powertrain_pm
            end
            obj.m            = params.m;
            obj.cda          = params.cda;
            obj.cla          = params.cla;
            obj.cop          = params.cop;
            obj.weight_r     = params.weight_r;
            obj.wheel_radius = params.wheel_radius;
            if isfield(params,'rho'),   obj.rho   = params.rho;   end
            if isfield(params,'g'),     obj.g     = params.g;     end
            if isfield(params,'c_roll'),obj.c_roll= params.c_roll;end
            if isfield(params,'max_brake_torque'), obj.max_brake_torque = params.max_brake_torque; end
            if isfield(params,'tire_limited_brake'), obj.tire_limited_brake = params.tire_limited_brake; end
            obj.tireMu = tireMu;
            obj.ptrain = ptrain;
        end

        % loads & aero/drag
        function DF = downforce(obj, v)
            DF = 0.5*obj.rho*obj.cla.*(v.^2);
        end
        function Drag = drag(obj, v)
            Drag = 0.5*obj.cda*obj.rho.*(v.^2);
        end
        function Fz_ax = axleLoads(obj, v)
            % returns [Fz_front_ax, Fz_rear_ax]
            wf_static = 1 - obj.weight_r;
            bf = obj.cop; br = 1 - bf;
            DF = obj.downforce(v);
            Fz_front_ax = wf_static*obj.m*obj.g + bf*DF;
            Fz_rear_ax  = obj.weight_r*obj.m*obj.g + br*DF;
            Fz_ax = [Fz_front_ax, Fz_rear_ax];
        end

        % tire capacity with μ(Fz_tire)
        function [Fcap_front_ax, Fcap_rear_ax, Fcap_total] = longCapacities(obj, v)
            Fz_ax = obj.axleLoads(v); % Nx2
            Fz_f_t = 0.5*Fz_ax(:,1);  % per tire
            Fz_r_t = 0.5*Fz_ax(:,2);
            lbf = 1/4.4482216153;
            mu_f = obj.tireMu.mu_of_Fz_lbf(Fz_f_t*lbf);
            mu_r = obj.tireMu.mu_of_Fz_lbf(Fz_r_t*lbf);
            Fcap_front_ax = 2.*mu_f.*Fz_f_t;
            Fcap_rear_ax  = 2.*mu_r.*Fz_r_t;
            Fcap_total    = Fcap_front_ax + Fcap_rear_ax;
        end

        % longitudinal force balances
        function Fx_roll = rollingDrag(obj, v)
            Fx_roll = obj.c_roll .* (obj.m*obj.g + obj.downforce(v));
        end

        function a_accel = accelMax(obj, v)
            [~, Fcap_rear_ax, ~] = obj.longCapacities(v);
            F_engine = min(obj.ptrain.bestWheelForce(v), Fcap_rear_ax);
            Fx = F_engine - (obj.drag(v) + obj.rollingDrag(v));
            a_accel = max(0, Fx ./ obj.m);
        end

        function a_brake = brakeMax(obj, v)
            [~, ~, Fcap_total] = obj.longCapacities(v);
            if obj.tire_limited_brake
                brake_cap_force = Inf;
            else
                brake_cap_force = obj.max_brake_torque / obj.wheel_radius;
            end
            Ft = min(brake_cap_force, Fcap_total);   % available tire Fx
            Fx = Ft + (obj.drag(v) + obj.rollingDrag(v));  % drag adds
            a_brake = max(0, Fx ./ obj.m);
        end

        function [accel_F, decel_F, v_grid, acc_vec, dec_vec] = accelInterpolants(obj, vmax, npts)
            if nargin < 2, vmax = 33; end
            if nargin < 3, npts = 901; end
            v_grid = linspace(0, vmax, npts).';
            acc_vec = obj.accelMax(v_grid);
            dec_vec = obj.brakeMax(v_grid);
            accel_F = griddedInterpolant(v_grid, acc_vec, 'linear', 'nearest');
            decel_F = griddedInterpolant(v_grid, dec_vec, 'linear', 'nearest');
        end

        %lateral limit solver
        function v_lat = lateralLimit(obj, track, vmax)
            if nargin < 3, vmax = 33; end
        
            s = track.s(:); 
            k = track.k(:);
            r_i = 1 ./ max(1e-9, abs(k));      % radii for each point
            Np = numel(r_i);
        
            % 1) build a speed grid
            Nv = 400; %number of samples
            v_grid = linspace(0, vmax, Nv).'; 
        
            % 2) precompute stuff that only depends on v
            % total normal
            Fz_tot = obj.m*obj.g + obj.downforce(v_grid); % [N]
            % per-tire load in lbf
            Fz_tire_lbf = (Fz_tot./4) / 4.4482216153;
            % mu at that load
            mu_v = obj.tireMu.mu_of_Fz_lbf(Fz_tire_lbf);  % unitless
            % total available combined force
            Ft_v = mu_v .* Fz_tot; % [N]
            % aero + rolling drag at that v
            Fx_v = obj.drag(v_grid) + obj.rollingDrag(v_grid); % [N]
        
            % 3) for each track point, test the inequality on this grid and pick the highest v that passes
            v_lat = zeros(Np,1);
            for i = 1:Np
                r = r_i(i);
                if r > 1e4
                    % basically a straight
                    v_lat(i) = vmax;
                    continue;
                end
                % lateral demand at all speeds for this radius
                Fy_v = obj.m .* (v_grid.^2) ./ r;        % [N]
                % check feasibility: Fy^2 + Fx^2 <= Ft^2
                feasible = (Fy_v.^2 + Fx_v.^2) <= (Ft_v.^2);
                % find last feasible speed
                idx = find(feasible, 1, 'last');
                if isempty(idx)
                    v_lat(i) = 0;
                else
                    v_lat(i) = v_grid(idx);
                end
            end
        end

    end
end
