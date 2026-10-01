classdef ParamSet
    properties
        car
        longVel
        maxLatFlag
        maxLatx
        maxLatLatAccel
        maxLatLongAccel
        maxLatX0
        maxLatValid = true
        maxLatResidual = NaN
        maxLatFunctionEvaluations = 0
        maxLatEquationCalls = 0
        
        maxLongFlag
        maxLongAccelx
        maxLongLongAccel
        maxLongLatAccel
        maxLongAccelx0
        maxLongValid = true
        maxLongResidual = NaN
        maxLongFunctionEvaluations = 0
        maxLongEquationCalls = 0
        
        maxBrakeFlag
        maxBrakeDecelx
        maxBrakeLongDecel
        maxBrakeLatAccel
        maxBrakingDecelx0
        maxBrakeValid = true
        maxBrakeResidual = NaN
        maxBrakeFunctionEvaluations = 0
        maxBrakeEquationCalls = 0
        
    end
    methods
        function obj = ParamSet(inputCar,longVel)
            if nargin > 0
                obj.car = inputCar;
                obj.longVel = longVel;
            end
        end
        function obj = setMaxLatParams(obj,x_ss,latAccel,longAccel,x0,diagnostics)
            %steady state state vector
            obj.maxLatx = x_ss;
            obj.maxLatFlag = x_ss(1);
            %accelerations
            obj.maxLatLatAccel = latAccel;
            obj.maxLatLongAccel = longAccel;
            obj.maxLatX0 = x0;
            if nargin >= 6 && isstruct(diagnostics)
                obj.maxLatValid = diagnostics.valid;
                obj.maxLatResidual = diagnostics.residual;
                obj.maxLatFunctionEvaluations = diagnosticCount(diagnostics, ...
                    'functionEvaluations');
                obj.maxLatEquationCalls = diagnosticCount(diagnostics,'equationCalls');
            end
        end
        function obj = setMaxAccelParams(obj,xAccel,longAccel,latAccel,longAccelx0,diagnostics)
            obj.maxLongFlag = xAccel(1);
            obj.maxLongAccelx = xAccel;
            obj.maxLongLongAccel = longAccel;
            obj.maxLongLatAccel = latAccel;
            obj.maxLongAccelx0 = longAccelx0;
            if nargin >= 6 && isstruct(diagnostics)
                obj.maxLongValid = diagnostics.valid;
                obj.maxLongResidual = diagnostics.residual;
                obj.maxLongFunctionEvaluations = diagnosticCount(diagnostics, ...
                    'functionEvaluations');
                obj.maxLongEquationCalls = diagnosticCount(diagnostics,'equationCalls');
            end
        end
        function obj = setMaxDecelParams(obj,xBraking,longDecel,latAccel,longDecelx0,diagnostics)
            obj.maxBrakeFlag = xBraking(1);
            obj.maxBrakeDecelx = xBraking;
            obj.maxBrakeLongDecel = longDecel;
            obj.maxBrakeLatAccel = latAccel;
            obj.maxBrakingDecelx0 = longDecelx0;
            if nargin >= 6 && isstruct(diagnostics)
                obj.maxBrakeValid = diagnostics.valid;
                obj.maxBrakeResidual = diagnostics.residual;
                obj.maxBrakeFunctionEvaluations = diagnosticCount(diagnostics, ...
                    'functionEvaluations');
                obj.maxBrakeEquationCalls = diagnosticCount(diagnostics,'equationCalls');
            end
        end
    end
end

function value = diagnosticCount(diagnostics,name)
value = 0;
if isfield(diagnostics,name) && isscalar(diagnostics.(name)) && ...
        isfinite(diagnostics.(name))
    value = diagnostics.(name);
end
end

