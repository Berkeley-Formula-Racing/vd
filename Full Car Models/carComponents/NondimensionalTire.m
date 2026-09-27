classdef NondimensionalTire
    %NONDIMENSIONALTIRE Vehicle-facing adapter for the TTC target/donor model.
    %   Target lateral behavior comes from the exact target-size TTC data;
    %   longitudinal and combined-slip shape come from the configured larger
    %   same-compound donor. Unsupported slip queries are clamped for the
    %   vehicle solver and reported through EVALUATE's diagnostic struct.

    properties
        p_i
        friction_scaling_factor
        model
        uncertainty_mode
    end

    methods
        function obj = NondimensionalTire(model,p_i,friction_scaling_factor,uncertainty_mode)
            if nargin < 4 || isempty(uncertainty_mode)
                uncertainty_mode = "nominal";
            end
            if ~isstruct(model) || ~isscalar(model)
                error('NondimensionalTire:badModel', ...
                    'MODEL must be a scalar nondimensional tire-model struct.');
            end
            required = {'target_lateral','donor_longitudinal', ...
                'target_lateral_capacity','longitudinal_mu_scale', ...
                'coupling_exponent'};
            if ~all(isfield(model,required))
                error('NondimensionalTire:badModel', ...
                    'MODEL is missing one or more required model fields.');
            end

            packageRoot = fullfile(fileparts(fileparts(fileparts( ...
                mfilename('fullpath')))), 'Magic Formula','experimental', ...
                'lc0_nd_tire');
            if isfolder(packageRoot), addpath(packageRoot); end
            if ~isfield(model,'uncertainty_scenarios')
                base = struct('rhoMu',getModelField(model,'rho_mu',1), ...
                    'rhoStiff',getModelField(model,'rho_stiff',1), ...
                    'couplingExponent',model.coupling_exponent);
                model.uncertainty_scenarios = lc0NDUncertaintyScenarios(base);
            end
            if nargin < 2 || isempty(p_i), p_i = 12; end
            if nargin < 3 || isempty(friction_scaling_factor)
                friction_scaling_factor = 1;
            end
            if ~isscalar(p_i) || ~isfinite(p_i) || p_i <= 0
                error('NondimensionalTire:badPressure', ...
                    'p_i must be a finite positive pressure in psi.');
            end
            if ~isscalar(friction_scaling_factor) || ...
                    ~isfinite(friction_scaling_factor) || friction_scaling_factor <= 0
                error('NondimensionalTire:badScale', ...
                    'friction_scaling_factor must be finite and positive.');
            end
            obj.model = model;
            obj.p_i = p_i;
            obj.friction_scaling_factor = friction_scaling_factor;
            obj.uncertainty_mode = string(uncertainty_mode);
            validateScenario(obj.model,obj.uncertainty_mode);
        end

        function out = F_x(obj,alpha,kappa,F_z,gamma)
            if nargin < 5, gamma = zeros(size(alpha)); end
            [Fx,~,~] = obj.evaluate(alpha,kappa,F_z,gamma);
            out = Fx;
        end

        function out = F_y(obj,alpha,kappa,F_z,gamma)
            if nargin < 5, gamma = zeros(size(alpha)); end
            [~,Fy,~] = obj.evaluate(alpha,kappa,F_z,gamma);
            out = Fy;
        end

        function [Fx,Fy,info] = evaluate(obj,alpha,kappa,F_z,gamma)
            if nargin < 5, gamma = zeros(size(alpha)); end
            kappa = expandScalar(kappa,size(alpha));
            F_z = expandScalar(F_z,size(alpha));
            gamma = expandScalar(gamma,size(alpha));
            if ~(isequal(size(alpha),size(kappa),size(F_z)) && ...
                    isequal(size(alpha),size(gamma)))
                error('NondimensionalTire:sizeMismatch', ...
                    'alpha, kappa, F_z, and gamma must have matching sizes.');
            end
            if any(~isfinite(F_z(:))) || any(F_z(:) < 0)
                error('NondimensionalTire:negativeLoad', ...
                    'F_z must be finite and nonnegative.');
            end

            % The experimental evaluator is deliberately kept in its own
            % package. The constructor adds it lazily so a vehicle
            % simulation can select this tire without a global path-order
            % assumption.
            options = struct('outOfRange',"clamp",'scenario',obj.uncertainty_mode, ...
                'pressurePsi',obj.p_i,'camberDeg',gamma);
            [Fx,Fy,info] = lc0NDEvaluate(obj.model,alpha,kappa,F_z,options);
            Fx = Fx.*obj.friction_scaling_factor;
            Fy = Fy.*obj.friction_scaling_factor;

            info.pressure_psi = obj.p_i;
            info.camber_deg = gamma;
            info.is_camber_unmodeled = false(size(gamma));
            info.is_pressure_unmodeled = false(size(gamma));
            if isfield(obj.model.target_lateral,'options') && ...
                    isfield(obj.model.target_lateral.options,'pressure_psi')
                referencePressure = obj.model.target_lateral.options.pressure_psi;
                info.is_pressure_unmodeled = abs(obj.p_i-referencePressure) > 1e-9;
            end
            if any(abs(gamma(:)) > 1e-9)
                info.is_camber_unmodeled = true(size(gamma));
                info.is_extrapolated = info.is_extrapolated | ...
                    abs(gamma) > 1e-9;
            end
            info.is_extrapolated = info.is_extrapolated | ...
                info.is_pressure_unmodeled;
            info.clamp_count = sum(info.is_extrapolated(:));
            info.query_count = numel(F_z);
        end
    end
end

function value = expandScalar(value,targetSize)
if isscalar(value) && prod(targetSize) > 1
    value = repmat(value,targetSize);
end
end

function validateScenario(model,name)
if isfield(model,'uncertainty_scenarios') && ~isempty(model.uncertainty_scenarios)
    names = string({model.uncertainty_scenarios.name});
    if ~any(names == name)
        error('NondimensionalTire:unknownScenario', ...
            'Unknown uncertainty scenario %s.',name);
    end
elseif name ~= "nominal"
    error('NondimensionalTire:unknownScenario', ...
        'This model only contains the nominal scenario.');
end
end

function value = getModelField(model,name,defaultValue)
if isfield(model,name) && ~isempty(model.(name))
    value = model.(name);
else
    value = defaultValue;
end
end
