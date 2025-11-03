classdef powertrain_pm
    properties
        gears double     
        final_drive double
        eta double = 0.97
        wheel_radius double
        rpm_redline double
        rpm_idle double
        tqF 
    end
    methods
        function obj = powertrain_pm(gears, final_drive, eta, wheel_radius, ...
                                  rpm_redline, rpm_idle, rpm_vec, tq_lbft_vec)
            obj.gears = gears(:).';
            obj.final_drive = final_drive;
            obj.eta = eta;
            obj.wheel_radius = wheel_radius;
            obj.rpm_redline = rpm_redline;
            obj.rpm_idle = rpm_idle;

            % build torque interpolant 
            tqNm = tq_lbft_vec(:) * 1.3558179483314004; %lbf to N
            rpm  = rpm_vec(:);
            [rpm, iu] = unique(rpm,'stable'); tqNm = tqNm(iu);
            [rpm, is] = sort(rpm);            tqNm = tqNm(is);
            F = griddedInterpolant(rpm, tqNm, 'pchip', 'linear');
            obj.tqF = @(x) max(0, F(x));
        end

        function Fw = wheelForceForGear(obj, v, gr)
            % v in m/s
            rpm = (v(:)./obj.wheel_radius) .* gr .* obj.final_drive * 60/(2*pi);
            rpm = min(obj.rpm_redline, max(obj.rpm_idle, rpm));
            Tq  = obj.tqF(rpm);                            % N·m
            Fw  = (Tq .* gr .* obj.final_drive .* obj.eta) ./ obj.wheel_radius; % N
        end

        function Fw = bestWheelForce(obj, v)
            % pick max across all gears
            WG = arrayfun(@(gr) obj.wheelForceForGear(v, gr), obj.gears, 'UniformOutput', false);
            Fw = max(cell2mat(WG), [], 2);
        end
    end
end
