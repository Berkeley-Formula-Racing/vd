function [car,derived] = buildCar(setup,role,assets)
%BUILDCAR Build one fresh QSS lap or acceleration Car from editable sources.
if nargin < 2 || isempty(role), role = "lap"; end
if nargin < 3, assets = []; end

if isRampProfile(setup)
    [car,derived] = buildRampCar(setup,role,assets);
    return
end

[q,identity] = validateSetup(setup);
role = validateRole(role);
if isempty(assets) || (isstruct(assets) && isscalar(assets) && ...
        isempty(fieldnames(assets)))
    assets = vehicleConfig.loadAssets(setup);
end
assets = validateAssets(setup,identity,assets);

if role == "lap"
    map = [];
    if setup.aeroMap.enabled
        map = assets.map;
    end
    aero = Aero(q.cda,q.cla,q.distribution,q.cla_p_deg_p,q.D_p_deg_p,map);
    driverWeight = q.driver_weight;
else
    aero = Aero(q.accel_cda,q.accel_cla,q.distribution, ...
        q.acc_cla_p_deg_p,q.acc_D_p_deg_p);
    driverWeight = q.accel_driver_weight;
end

fixed = setup.fixed;
powertrain = Powertrain(q.redline,q.shift_point,fixed.gears, ...
    fixed.primary_reduction,fixed.torque_fn,q.shift_time,q.final_drive, ...
    q.wheel_radius,q.drivetrain_efficiency,q.G_d1,q.G_d2_overrun, ...
    q.G_d2_driving,q.brake_distribution,q.max_braking_torque);
tire = buildTire(fixed,q,assets);
rideHeightAero = makeRideHeightAero(q);
camberKinematics = makeCamberKinematics(q,assets.camberModelData);
car = buildCarFromParameterSet(q,driverWeight,aero,powertrain,tire, ...
    rideHeightAero,camberKinematics,q.R_sf);

derived = struct( ...
    'role',role, ...
    'driverWeight_kg',driverWeight, ...
    'totalMass_kg',car.M, ...
    'wheelbase_m',car.W_b, ...
    'frontAxleDistance_m',car.l_f, ...
    'rearAxleDistance_m',car.l_r, ...
    'wheelRadius_m',car.R, ...
    'wheelRateFront_Npm',rideHeightAero.wheel_rate_front_Npm, ...
    'wheelRateRear_Npm',rideHeightAero.wheel_rate_rear_Npm, ...
    'R_sf',car.R_sf, ...
    'aeroMapEnabled',car.aero.hasMap(), ...
    'assetIdentity',assets.identity);
if isfield(assets,'fingerprints')
    derived.assetFingerprints = assets.fingerprints;
else
    derived.assetFingerprints = struct();
end
end

function [q,identity] = validateSetup(setup)
if ~isstruct(setup) || ~isscalar(setup)
    error('vehicleConfig:badSchema','setup must be one scalar configuration struct.');
end
expectedTop = {'schemaVersion','profile','id','label','parameters', ...
    'fixed','aeroMap','suspension'};
validateKnownFields(setup,expectedTop,'setup');
requireFields(setup,expectedTop,'setup');
if ~isnumeric(setup.schemaVersion) || ~isscalar(setup.schemaVersion) || ...
        ~isreal(setup.schemaVersion) || ~isfinite(double(setup.schemaVersion)) || ...
        double(setup.schemaVersion) ~= 1
    error('vehicleConfig:badSchema','Only scalar schemaVersion 1 is supported.');
end
if ~isTextScalar(setup.profile) || lower(string(setup.profile)) ~= "qss"
    error('vehicleConfig:badSchema','buildCar supports the QSS profile only.');
end
if ~isTextScalar(setup.id) || ~isTextScalar(setup.label)
    error('vehicleConfig:badSchema','setup id and label must be scalar text.');
end

names = vehicleConfig.parameterNames();
if ~isstruct(setup.parameters) || ~isscalar(setup.parameters)
    error('vehicleConfig:badSchema','setup.parameters must be one scalar struct.');
end
validateKnownFields(setup.parameters,names,'setup.parameters');
requireFields(setup.parameters,names,'setup.parameters');
q = struct();
for k = 1:numel(names)
    name = names{k};
    value = setup.parameters.(name);
    if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ...
            ~isfinite(double(value))
        invalidParameter(['parameters.' name], ...
            'must be a finite, real numeric scalar');
    end
    q.(name) = double(value);
end
validateSourceConstraints(q);

fixedNames = {'gears','primary_reduction','torque_fn','Fx_parameters', ...
    'Fy_parameters','friction_scaling_factor','model_type','model_artifact', ...
    'model_uncertainty'};
mapNames = {'enabled','path','claScale','cdaScale','copOffset'};
suspensionNames = {'rollSplitPolicy','frontArbStiffness_NmPerRad', ...
    'rearArbStiffness_NmPerRad'};
validateNestedStruct(setup,'fixed',fixedNames);
validateNestedStruct(setup,'aeroMap',mapNames);
validateNestedStruct(setup,'suspension',suspensionNames);
fixed = setup.fixed;
map = setup.aeroMap;
suspension = setup.suspension;

if ~isnumeric(fixed.gears) || isempty(fixed.gears) || ~isvector(fixed.gears) || ...
        ~isreal(fixed.gears) || any(~isfinite(double(fixed.gears(:)))) || ...
        any(double(fixed.gears(:)) <= 0)
    invalidParameter('fixed.gears','must be a finite positive numeric vector');
end
validatePositiveScalar(fixed.primary_reduction,'fixed.primary_reduction');
validateFiniteNumericArray(fixed.torque_fn,'fixed.torque_fn');
if ~ismatrix(fixed.torque_fn) || size(fixed.torque_fn,1) ~= 2 || ...
        size(fixed.torque_fn,2) < 2
    invalidParameter('fixed.torque_fn', ...
        'must be a two-row torque curve with at least two samples');
end
validateFiniteNumericArray(fixed.Fx_parameters,'fixed.Fx_parameters');
validateFiniteNumericArray(fixed.Fy_parameters,'fixed.Fy_parameters');
validatePositiveScalar(fixed.friction_scaling_factor, ...
    'fixed.friction_scaling_factor');
if ~isTextScalar(fixed.model_type)
    invalidParameter('fixed.model_type','must be scalar text');
end
modelType = lower(string(fixed.model_type));
if ~any(modelType == ["legacy","lc0_nd","nondimensional"])
    invalidParameter('fixed.model_type','must be legacy or nondimensional');
end
if ~isTextScalar(fixed.model_artifact) || ~isTextScalar(fixed.model_uncertainty)
    error('vehicleConfig:badSchema', ...
        'fixed.model_artifact and fixed.model_uncertainty must be scalar text.');
end
if ~any(lower(string(fixed.model_uncertainty)) == ["low","nominal","high"])
    invalidParameter('fixed.model_uncertainty','must be low, nominal, or high');
end

validateFlag(map.enabled,'aeroMap.enabled');
if ~isTextScalar(map.path)
    error('vehicleConfig:badSchema','aeroMap.path must be scalar text.');
end
validatePositiveScalar(map.claScale,'aeroMap.claScale');
validatePositiveScalar(map.cdaScale,'aeroMap.cdaScale');
validateFiniteScalar(map.copOffset,'aeroMap.copOffset');
if ~isTextScalar(suspension.rollSplitPolicy) || ...
        lower(string(suspension.rollSplitPolicy)) ~= "explicit"
    error('vehicleConfig:badSchema', ...
        'QSS construction requires the explicit rollSplitPolicy.');
end
validateFiniteScalar(suspension.frontArbStiffness_NmPerRad, ...
    'suspension.frontArbStiffness_NmPerRad');
validateFiniteScalar(suspension.rearArbStiffness_NmPerRad, ...
    'suspension.rearArbStiffness_NmPerRad');
identity = setupAssetIdentity(setup);
end

function validateSourceConstraints(q)
positiveNames = {'wheelbase','wheel_radius','spring_rate_front_lb_in', ...
    'spring_rate_rear_lb_in','motion_ratio_front','motion_ratio_rear'};
for k = 1:numel(positiveNames)
    name = positiveNames{k};
    if q.(name) <= 0
        invalidParameter(['parameters.' name],'must be positive');
    end
end
if q.weight_dist <= 0 || q.weight_dist >= 1
    invalidParameter('parameters.weight_dist','must be strictly between zero and one');
end
end

function role = validateRole(role)
if ~isTextScalar(role)
    error('vehicleConfig:badRole','role must be lap or accel.');
end
role = lower(string(role));
if ~any(role == ["lap","accel"])
    error('vehicleConfig:badRole','role must be lap or accel.');
end
end

function assets = validateAssets(setup,identity,assets)
required = {'identity','map','camberModelData','camberRatioData'};
if ~isstruct(assets) || ~isscalar(assets) || ~all(isfield(assets,required)) || ...
        ~isstruct(assets.identity) || ~isscalar(assets.identity)
    assetMismatch();
end
identityFields = fieldnames(identity);
if ~all(isfield(assets.identity,identityFields))
    assetMismatch();
end
for k = 1:numel(identityFields)
    name = identityFields{k};
    if ~isequaln(assets.identity.(name),identity.(name))
        assetMismatch();
    end
end
if ~isstruct(assets.camberModelData) || ~isscalar(assets.camberModelData) || ...
        ~all(isfield(assets.camberModelData,{'betaL','betaR'})) || ...
        ~isstruct(assets.camberRatioData) || ~isscalar(assets.camberRatioData) || ...
        ~all(isfield(assets.camberRatioData, ...
            {'camber2indices','camber4indices','camber2ratio','camber4ratio'}))
    assetMismatch();
end
if setup.aeroMap.enabled && ~isa(assets.map,'AeroMap')
    assetMismatch();
end
if ~isempty(assets.map)
    if ~isa(assets.map,'AeroMap') || ...
            string(assets.map.sourcePath) ~= identity.aeroMapPath
        assetMismatch();
    end
end
modelType = lower(string(setup.fixed.model_type));
if any(modelType == ["lc0_nd","nondimensional"]) && ...
        (~isfield(assets,'tireModel') || ~isstruct(assets.tireModel) || ...
        ~isscalar(assets.tireModel))
    assetMismatch();
end
end

function tire = buildTire(fixed,q,assets)
modelType = lower(string(fixed.model_type));
if modelType == "legacy"
    tireParams = struct( ...
        'model_type',modelType, ...
        'model_artifact',fixed.model_artifact, ...
        'model_uncertainty',fixed.model_uncertainty, ...
        'p_i',q.p_i, ...
        'Fx_parameters',fixed.Fx_parameters, ...
        'Fy_parameters',fixed.Fy_parameters, ...
        'friction_scaling_factor',fixed.friction_scaling_factor, ...
        'camber_data',assets.camberRatioData);
    tire = buildTireModel(tireParams);
else
    tire = NondimensionalTire(assets.tireModel,q.p_i, ...
        fixed.friction_scaling_factor,fixed.model_uncertainty);
end
end

function config = makeRideHeightAero(q)
lbInToNpm = 175.126835;
config = struct( ...
    'enabled',true, ...
    'spring_rate_front_lb_in',q.spring_rate_front_lb_in, ...
    'spring_rate_rear_lb_in',q.spring_rate_rear_lb_in, ...
    'motion_ratio_front',q.motion_ratio_front, ...
    'motion_ratio_rear',q.motion_ratio_rear, ...
    'wheel_rate_front_Npm',q.spring_rate_front_lb_in*lbInToNpm*q.motion_ratio_front^2, ...
    'wheel_rate_rear_Npm',q.spring_rate_rear_lb_in*lbInToNpm*q.motion_ratio_rear^2, ...
    'static_front_ride_height_in',q.static_front_ride_height_in, ...
    'static_rear_ride_height_in',q.static_rear_ride_height_in, ...
    'map_reference_front_ride_height_in',q.map_reference_front_ride_height_in, ...
    'map_reference_rear_ride_height_in',q.map_reference_rear_ride_height_in);
end

function config = makeCamberKinematics(q,camberModelData)
config = struct( ...
    'roll_gradient_deg_per_g',q.roll_gradient_deg_per_g, ...
    'rear_roll_camber_outer_deg_per_deg',q.rear_roll_camber_outer_deg_per_deg, ...
    'rear_roll_camber_inner_deg_per_deg',q.rear_roll_camber_inner_deg_per_deg, ...
    'ride_camber_front_deg_per_in',q.ride_camber_front_deg_per_in, ...
    'ride_camber_rear_deg_per_in',q.ride_camber_rear_deg_per_in, ...
    'camberModelData',camberModelData);
end

function validateNestedStruct(setup,name,expected)
value = setup.(name);
if ~isstruct(value) || ~isscalar(value)
    error('vehicleConfig:badSchema','setup.%s must be one scalar struct.',name);
end
validateKnownFields(value,expected,['setup.' name]);
requireFields(value,expected,['setup.' name]);
end

function validateKnownFields(value,expected,label)
unknown = setdiff(fieldnames(value),expected,'stable');
if ~isempty(unknown)
    error('vehicleConfig:unknownField', ...
        '%s contains unsupported field %s.',label,unknown{1});
end
end

function requireFields(value,expected,label)
missing = setdiff(expected,fieldnames(value),'stable');
if ~isempty(missing)
    error('vehicleConfig:badSchema','%s is missing required field %s.', ...
        label,missing{1});
end
end

function tf = isTextScalar(value)
tf = ischar(value) && (isrow(value) || isempty(value));
if isstring(value)
    tf = isscalar(value) && ~ismissing(value);
end
end

function validateFlag(value,name)
if ~(isnumeric(value) || islogical(value)) || ~isscalar(value) || ...
        ~isreal(value) || ~isfinite(double(value)) || ...
        ~any(double(value) == [0 1])
    invalidParameter(name,'must be one finite zero-or-one flag');
end
end

function validatePositiveScalar(value,name)
validateFiniteScalar(value,name);
if double(value) <= 0
    invalidParameter(name,'must be positive');
end
end

function validateFiniteScalar(value,name)
if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ...
        ~isfinite(double(value))
    invalidParameter(name,'must be one finite real numeric scalar');
end
end

function validateFiniteNumericArray(value,name)
if ~isnumeric(value) || isempty(value) || ~isreal(value) || ...
        any(~isfinite(double(value(:))))
    invalidParameter(name,'must be a non-empty finite real numeric array');
end
end

function invalidParameter(name,requirement)
error('vehicleConfig:invalidParameter', ...
    'setup.%s %s.',name,requirement);
end

function assetMismatch()
error('vehicleConfig:assetMismatch', ...
    'The supplied asset context belongs to a different setup identity. Reload assets for this setup.');
end

function tf = isRampProfile(setup)
tf = isstruct(setup) && isscalar(setup) && isfield(setup,'profile') && ...
    isTextScalar(setup.profile) && lower(string(setup.profile)) == "ramp";
end

function [car,derived] = buildRampCar(setup,role,assets)
[p,fixed,mapSpec,suspension] = validateRampSetup(setup);
role = validateRole(role);
if role ~= "accel"
    error('vehicleConfig:badRole', ...
        'The ramp profile supports the accel role only.');
end

roll = rampSpeed.deriveRollStiffness(p.frontSpringRate_lb_in, ...
    p.rearSpringRate_lb_in,p.motion_ratio_front,p.motion_ratio_rear, ...
    p.track_width,suspension.frontArbStiffness_NmPerRad, ...
    suspension.rearArbStiffness_NmPerRad);
if isempty(assets) || (isstruct(assets) && isscalar(assets) && ...
        isempty(fieldnames(assets)))
    assets = struct();
elseif ~isstruct(assets) || ~isscalar(assets)
    assetMismatch();
end

map = [];
if mapSpec.enabled
    if isfield(assets,'map') && isa(assets.map,'AeroMap') && ...
            string(assets.map.sourcePath) == string(mapSpec.path)
        map = assets.map;
    else
        map = AeroMap(char(mapSpec.path));
    end
    if map.claScale ~= mapSpec.claScale || ...
            map.cdaScale ~= mapSpec.cdaScale || ...
            map.copOffset ~= mapSpec.copOffset
        map = map.withCorrections(mapSpec.claScale,mapSpec.cdaScale, ...
            mapSpec.copOffset);
    end
end

aq = p;
aq.weight_dist = p.rearWeightDistribution_percent/100;
aq.spring_rate_front_lb_in = p.frontSpringRate_lb_in;
aq.spring_rate_rear_lb_in = p.rearSpringRate_lb_in;
aq.static_front_ride_height_in = p.frontRideHeight_in;
aq.static_rear_ride_height_in = p.rearRideHeight_in;
aero = Aero(aq.accel_cda,aq.accel_cla,aq.aero_distribution, ...
    aq.accel_cla_p_deg_p,aq.accel_D_p_deg_p,map);
powertrain = Powertrain(aq.redline,aq.shift_point,fixed.gears, ...
    fixed.primary_reduction,fixed.torque_fn,aq.shift_time,aq.final_drive, ...
    aq.wheel_radius,aq.drivetrain_efficiency,aq.G_d1,aq.G_d2_overrun, ...
    aq.G_d2_driving,aq.brake_distribution,aq.max_braking_torque);
camberRatios = [];
if isfield(assets,'camberRatios')
    camberRatios = assets.camberRatios;
elseif isfield(assets,'camberRatioData')
    camberRatios = assets.camberRatioData;
end
tire = Tire2(aq.p_i,fixed.Fx_parameters,fixed.Fy_parameters, ...
    fixed.friction_scaling_factor,camberRatios);
rideHeightAero = makeRideHeightAero(aq);
camberKinematics = struct( ...
    'roll_gradient_deg_per_g',aq.roll_gradient_deg_per_g, ...
    'rear_roll_camber_outer_deg_per_deg',aq.rear_roll_camber_outer_deg_per_deg, ...
    'rear_roll_camber_inner_deg_per_deg',aq.rear_roll_camber_inner_deg_per_deg, ...
    'ride_camber_front_deg_per_in',aq.ride_camber_front_deg_per_in, ...
    'ride_camber_rear_deg_per_in',aq.ride_camber_rear_deg_per_in);
if isfield(assets,'camberModelData')
    camberKinematics.camberModelData = assets.camberModelData;
end

car = buildCarFromParameterSet(aq,p.driverWeight_kg,aero,powertrain,tire, ...
    rideHeightAero,camberKinematics,roll.R_sf);
car.k_rf = suspension.frontArbStiffness_NmPerRad;
car.k_rr = suspension.rearArbStiffness_NmPerRad;
car.k_f_b = roll.wheelRateFront_Npm;
car.k_r_b = roll.wheelRateRear_Npm;
car.k_f_r = roll.wheelRateFront_Npm;
car.k_r_r = roll.wheelRateRear_Npm;
car.rs_total = (roll.frontRollStiffness_NmPerRad + ...
    roll.rearRollStiffness_NmPerRad)*pi/180;

derived = roll;
derived.driverWeight_kg = p.driverWeight_kg;
derived.rearWeightDistribution_percent = p.rearWeightDistribution_percent;
derived.frontWeightDistribution_percent = ...
    100-p.rearWeightDistribution_percent;
derived.totalMass_kg = car.M;
derived.staticRearLoad_N = car.M*car.g* ...
    p.rearWeightDistribution_percent/100;
derived.staticFrontLoad_N = car.M*car.g* ...
    (1-p.rearWeightDistribution_percent/100);
if isfield(assets,'fingerprints')
    derived.assetFingerprints = assets.fingerprints;
end
end

function [p,fixed,map,suspension] = validateRampSetup(setup)
if ~isstruct(setup) || ~isscalar(setup)
    error('vehicleConfig:badSchema','setup must be one scalar configuration struct.');
end
expectedTop = {'schemaVersion','profile','id','label','parameters', ...
    'fixed','aeroMap','suspension'};
validateKnownFields(setup,expectedTop,'setup');
requireFields(setup,expectedTop,'setup');
if ~isnumeric(setup.schemaVersion) || ~isscalar(setup.schemaVersion) || ...
        ~isreal(setup.schemaVersion) || ~isfinite(double(setup.schemaVersion)) || ...
        double(setup.schemaVersion) ~= 1
    error('vehicleConfig:badSchema','Only scalar schemaVersion 1 is supported.');
end
if ~isTextScalar(setup.profile) || lower(string(setup.profile)) ~= "ramp"
    error('vehicleConfig:badSchema','buildCar supports ramp profile only on this path.');
end
if ~isTextScalar(setup.id) || ~isTextScalar(setup.label)
    error('vehicleConfig:badSchema','setup id and label must be scalar text.');
end

parameterNames = rampParameterNames();
if ~isstruct(setup.parameters) || ~isscalar(setup.parameters)
    error('vehicleConfig:badSchema','setup.parameters must be one scalar struct.');
end
validateKnownFields(setup.parameters,parameterNames,'setup.parameters');
requireFields(setup.parameters,parameterNames,'setup.parameters');
p = struct();
for k = 1:numel(parameterNames)
    name = parameterNames{k};
    value = setup.parameters.(name);
    if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ...
            ~isfinite(double(value))
        invalidParameter(['parameters.' name], ...
            'must be a finite, real numeric scalar');
    end
    p.(name) = double(value);
end
if p.mass <= 0 || p.wheelbase <= 0 || p.wheel_radius <= 0 || ...
        p.motion_ratio_front <= 0 || p.motion_ratio_rear <= 0
    invalidParameter('parameters','mass, geometry, and motion ratios must be positive');
end
if p.driverWeight_kg < 0 || p.rearWeightDistribution_percent <= 0 || ...
        p.rearWeightDistribution_percent >= 100
    invalidParameter('parameters', ...
        'driver weight and rear distribution must be within ramp setup bounds');
end

fixedNames = {'gears','primary_reduction','torque_fn','Fx_parameters', ...
    'Fy_parameters','friction_scaling_factor','model_type'};
validateNestedStruct(setup,'fixed',fixedNames);
fixed = setup.fixed;
if ~isnumeric(fixed.gears) || isempty(fixed.gears) || ~isvector(fixed.gears) || ...
        ~isreal(fixed.gears) || any(~isfinite(double(fixed.gears(:)))) || ...
        any(double(fixed.gears(:)) <= 0)
    invalidParameter('fixed.gears','must be a finite positive numeric vector');
end
validatePositiveScalar(fixed.primary_reduction,'fixed.primary_reduction');
validateFiniteNumericArray(fixed.torque_fn,'fixed.torque_fn');
validateFiniteNumericArray(fixed.Fx_parameters,'fixed.Fx_parameters');
validateFiniteNumericArray(fixed.Fy_parameters,'fixed.Fy_parameters');
validatePositiveScalar(fixed.friction_scaling_factor, ...
    'fixed.friction_scaling_factor');
if ~isTextScalar(fixed.model_type) || lower(string(fixed.model_type)) ~= "legacy"
    error('vehicleConfig:badSchema', ...
        'The ramp profile requires the legacy Tire2 model.');
end

mapNames = {'enabled','path','claScale','cdaScale','copOffset'};
suspensionNames = {'rollSplitPolicy','frontArbStiffness_NmPerRad', ...
    'rearArbStiffness_NmPerRad'};
validateNestedStruct(setup,'aeroMap',mapNames);
validateNestedStruct(setup,'suspension',suspensionNames);
map = setup.aeroMap;
suspension = setup.suspension;
validateFlag(map.enabled,'aeroMap.enabled');
if ~isTextScalar(map.path)
    error('vehicleConfig:badSchema','aeroMap.path must be scalar text.');
end
validatePositiveScalar(map.claScale,'aeroMap.claScale');
validatePositiveScalar(map.cdaScale,'aeroMap.cdaScale');
validateFiniteScalar(map.copOffset,'aeroMap.copOffset');
if logical(map.enabled) && strlength(string(map.path)) == 0
    error('vehicleConfig:badSchema','An enabled ramp aero map requires a path.');
end
if ~isTextScalar(suspension.rollSplitPolicy) || ...
        lower(string(suspension.rollSplitPolicy)) ~= "derived"
    error('vehicleConfig:badSchema', ...
        'Ramp construction requires the derived rollSplitPolicy.');
end
validateFiniteScalar(suspension.frontArbStiffness_NmPerRad, ...
    'suspension.frontArbStiffness_NmPerRad');
validateFiniteScalar(suspension.rearArbStiffness_NmPerRad, ...
    'suspension.rearArbStiffness_NmPerRad');
end