classdef Aero
    % aero forces
    % update to return aero forces on front axle, rear axle
    
    properties
        cda
        cla
        D_f
        D_r
        cla_p_deg_p % cla per degree of pitch
        D_p_deg_p % Distibution per degree of pitch
        map % optional AeroMap, indexed by F/R ride-height offsets in inches
    end
    properties (Constant)
        rho = 1.2
    end
    
    methods
        function obj = Aero(cda,cla,distribution,cla_p_deg_p,D_p_deg_p,map)
            obj.cda = cda;
            obj.cla = cla;
            obj.D_f = distribution;
            obj.D_r = 1-distribution;
            obj.cla_p_deg_p = cla_p_deg_p; % 0.5
            obj.D_p_deg_p = D_p_deg_p; % 0.0677
            if nargin >= 6 && ~isempty(map), obj.map = map; end
        end
        
        function out = lift(obj,long_vel)            
            out = obj.rho/2*(long_vel^2)*obj.cla;
        end
        function out = pd_lift(obj,long_vel,pitch)            
            out = obj.rho/2*(long_vel^2)*(obj.cla+pitch*obj.cla_p_deg_p);
        end

        function out = drag(obj,long_vel)
            out = obj.rho/2*(long_vel^2)*obj.cda;
        end
        function out = pd_drag(obj,long_vel,~)
            % no cda pitch-sensitivity parameter exists yet, so pitch is ignored
            out = obj.rho/2*(long_vel^2)*obj.cda;
        end

        function tf = hasMap(obj)
            tf = ~isempty(obj.map);
        end

        function out = coefficients(obj,frontOffsetIn,rearOffsetIn)
            %COEFFICIENTS Aero values at a ride-height state, or static fallback.
            if obj.hasMap()
                out = obj.map.evaluate(frontOffsetIn,rearOffsetIn);
            else
                out = struct('cla',obj.cla,'cda',obj.cda, ...
                    'D_f',obj.D_f,'D_r',obj.D_r, ...
                    'frontOffsetIn',frontOffsetIn,'rearOffsetIn',rearOffsetIn, ...
                    'outsideMap',false);
            end
        end

        function [cla,cda,D_f,D_r,outsideMap] = coefficientsNumeric(obj,frontOffsetIn,rearOffsetIn)
            %COEFFICIENTSNUMERIC Numeric hot-path equivalent of coefficients().
            if obj.hasMap()
                [cla,cda,D_f,D_r,outsideMap] = obj.map.evaluateNumeric( ...
                    frontOffsetIn,rearOffsetIn);
            else
                if isscalar(frontOffsetIn) && ~isscalar(rearOffsetIn)
                    frontOffsetIn = repmat(frontOffsetIn,size(rearOffsetIn));
                elseif isscalar(rearOffsetIn) && ~isscalar(frontOffsetIn)
                    rearOffsetIn = repmat(rearOffsetIn,size(frontOffsetIn));
                elseif ~isequal(size(frontOffsetIn),size(rearOffsetIn))
                    error('Aero:querySizeMismatch', ...
                        'Front and rear offsets must have equal sizes or be scalar.');
                end
                cla = repmat(obj.cla,size(frontOffsetIn));
                cda = repmat(obj.cda,size(frontOffsetIn));
                D_f = repmat(obj.D_f,size(frontOffsetIn));
                D_r = repmat(obj.D_r,size(frontOffsetIn));
                outsideMap = false(size(frontOffsetIn));
            end
        end
    end
    
end