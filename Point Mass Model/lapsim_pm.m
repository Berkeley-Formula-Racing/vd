classdef lapsim_pm
    properties
        epsv double = 0.2
    end
    methods
        function [lap_time, v_profile, v_latlim] = run(obj, car, track, vmax)
            if nargin < 4, vmax = 33; end

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
    end
end
