classdef BDI_Tire_Model
    properties
        tire
    end

    methods
        function obj = BDI_Tire_Model()
            gamma = 0; % camber angle
            p_i = 12; % pressure
            % these parameters are non-iterable
            load('Fx_combined_parameters_run38_30.mat'); % F_x combined magic formula parameters
            Fx_parameters = cell2mat(Xbestcell);
            load('Lapsim_Fy_combined_parameters_1965run15.mat'); % F_y combined magic formula parameters
            Fy_parameters = cell2mat(Xbestcell);
            friction_scaling_factor = 1.05*0.55; % scales tire forces to account for test/road surface difference
            obj.tire = BDITire2(gamma,p_i,Fx_parameters,Fy_parameters,friction_scaling_factor);
        end
        function [force, slip_ratio] = peak_longitudinal_force(obj, normal_force)
            %iterate through slip ratios to find peak longitudinal force
            kappa_range = 0.1:0.01:0.2;
            force_range = zeros(1,numel(kappa_range));
            for i = 1:numel(kappa_range)
                force_range(i) = F_x(obj.tire,0,kappa_range(i),normal_force);
            end
            %plot(kappa_range, forces);
            [force, index] = max(force_range);
            slip_ratio = kappa_range(index);
        end
    end
end
