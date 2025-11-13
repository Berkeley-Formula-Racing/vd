classdef lapsim_pm
    properties
        epsv double = 0.2
    end
    methods
        function [lap_time, v_profile, v_latlim] = run(obj, car, track, vmax)
            if nargin < 4, vmax = 30; end

            [accel_F, decel_F, ~, ~, ~] = car.accelInterpolants(vmax, 901);
            v_latlim = car.lateralLimit(track, vmax);

            s  = track.s(:);
            k  = track.k(:);
            N  = numel(s);
            ds = diff(s); if size(ds,1) < size(ds,2), ds = ds'; end

            v_profile = zeros(N,1);
            v_profile(1) = min(1.0, v_latlim(1));


            for i = 1:N-1
                a_here = accel_F( max(0, v_profile(i)) );
                v_next = sqrt( max(0, v_profile(i)^2 + 2*a_here*ds(i)) );
                v_profile(i+1) = min(v_next, v_latlim(i+1));
            end

            for i = N-1:-1:1
                ab_here = decel_F( max(0, v_profile(i+1)) );
                v_cap   = sqrt( max(0, v_profile(i+1)^2 + 2*ab_here*ds(i)) );
                v_profile(i) = min( v_profile(i), v_cap );
                v_profile(i) = min( v_profile(i), v_latlim(i) );
            end

            dt_segments = [ds ./ max(v_profile(1:end-1), obj.epsv); 0];
            lap_time = sum(dt_segments);
        end

        function skidpad = skidpad_pm(obj, car, R)
            if nargin < 3
                R = 7.625; 
            end  
            fun = @(v) mCentripetal(v) - FyTotal(v);
        
            v_guess = 12; 
            v_ss = fzero(fun, v_guess);
        
            g_lat = v_ss^2 / (R * car.g);
        
            skidpad.v_ss  = v_ss;
            skidpad.g     = g_lat;
            skidpad.time  = (2*pi*R) / v_ss;
        
            function F = mCentripetal(v)
                F = car.m * v^2 / R;
            end
        
            function Fy = FyTotal(v)
                Fz_total = car.m * car.g + 0.5 * car.rho * car.cla * v.^2;
                Fz_tire_N = Fz_total / 4;
                Fz_tire_lbf = Fz_tire_N / 4.4482216153;
                mu = car.tireMu.mu_of_Fz_lbf(Fz_tire_lbf);
                Fy_tire = mu .* Fz_tire_N;
                Fy = 4 * Fy_tire;
            end
        end

        function result = accel(obj, car)
            L_total = 75.3;   % 75 m timed + 0.3 m before start
            N = 1000;
        
            s = linspace(0, L_total, N)';
            k = zeros(N,1);
        
            track.s = s;
            track.k = k;
        
            lapsim = lapsim_pm();
            vmax = 80;   % high enough for accel pass
            [lap_time_full, v_profile, ~] = lapsim.run(car, track, vmax);
        
            % Time corresponding to first 0.3 m
            idx_start = find(s >= 0.3, 1, 'first');
        
            % Time up to that point
            ds = diff(s);
            v_seg = v_profile(1:end-1);
            dt_seg = ds ./ max(v_seg, lapsim.epsv);
            time_before = sum(dt_seg(1:idx_start-1));
        
            % official time
            accel_time = lap_time_full - time_before;
        
            % Speed at 75 m
            idx_75 = find(s >= 75.3, 1, 'first');
            v_final = v_profile(idx_75);
        
            result.time = accel_time;
            result.v_profile = v_profile;
            result.v_final   = v_final;
        end
    end
end
