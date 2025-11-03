classdef tire_pm
    properties
        coeff double = 2.727
        expn  double = -0.096
        scale double = 1.0   % track friction scale
    end
    methods
        function obj = tire_pm(coeff, expn, scale)
            if nargin >= 1, obj.coeff = coeff; end
            if nargin >= 2, obj.expn  = expn;  end
            if nargin >= 3, obj.scale = scale; end
        end
        function mu = mu_of_Fz_lbf(obj, Fz_lbf)
            mu = obj.scale .* obj.coeff .* (Fz_lbf./100).^obj.expn;
        end
    end
end
