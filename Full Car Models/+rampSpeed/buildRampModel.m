function model = buildRampModel(car,setupSpec,profile,speeds)
%BUILDRAMPMODEL Build a serializable, data-only Ramp Speed model snapshot.
%   MODEL = rampSpeed.buildRampModel(CAR,SETUPSPEC,PROFILE) records the
%   setup, scalar vehicle data, and model provenance without retaining CAR
%   or any of its component objects.
%
%   MODEL = rampSpeed.buildRampModel(...,SPEEDS) also compiles and stores
%   the continuous powertrain envelope for each requested speed.

if nargin < 1 || isempty(car) || ~isobject(car) || ~isscalar(car)
    error('rampSpeed:invalidCar', ...
        'car must be one scalar solver-ready Car object.');
end
if nargin < 2 || isempty(setupSpec)
    setupSpec = struct();
end
if ~isstruct(setupSpec) || ~isscalar(setupSpec)
    error('rampSpeed:invalidSetup', ...
        'setupSpec must be a scalar struct.');
end
if nargin < 3
    profile = [];
end
if nargin < 4
    speeds = [];
end

resolvedProfile = rampSpeed.resolveSolverProfile(profile);
serializedProfile = toDataOnly( ...
    rampSpeed.serializeSolverProfile(resolvedProfile));

setupSnapshot = serializeSetupSpec(setupSpec);
model = struct();
model.schemaVersion = 1;
model.modelKind = 'rampSpeedLite';
model.modelVersion = 'ramp-speed-lite-1';
model.setupId = textField(setupSpec,'id');
model.label = textField(setupSpec,'label');
model.setupLabel = model.label;
model.baselineVersion = textField(setupSpec,'baselineVersion');
model.powertrainModel = textField(serializedProfile,'powertrainModel');
model.solverProfile = serializedProfile;
model.setupSpec = setupSnapshot;
model.vehicle = vehicleSnapshot(car);
model.driver = driverSnapshot(setupSpec,model.vehicle);
model.suspension = suspensionSnapshot(car,setupSpec);
model.aero = aeroSnapshot(car,setupSpec);
model.powertrain = powertrainSnapshot(car,model.powertrainModel);
model.speed_mps = zeros(0,1);

if ~isempty(speeds)
    if ~isnumeric(speeds) || ~isreal(speeds) || ...
            any(~isfinite(double(speeds(:)))) || any(double(speeds(:)) <= 0)
        error('rampSpeed:invalidSpeed', ...
            'speeds must contain finite positive real numeric values.');
    end
    model.speed_mps = unique(double(speeds(:)),'stable');
    [compiled,failures] = rampSpeed.compilePowertrainEnvelope( ...
        car,model.speed_mps,true);
    model.continuousEnvelope = toDataOnly(compiled);
    model.continuousEnvelopeErrors = toDataOnly(failures);
end

assertDataOnly(model,'model');
end

function snapshot = serializeSetupSpec(setupSpec)
snapshot = struct();
textNames = {'id','label','baselineVersion','source','aeroMapId', ...
    'aeroMapLabel','aeroMapRelativePath','aeroMapFingerprint'};
for i = 1:numel(textNames)
    name = textNames{i};
    value = textField(setupSpec,name);
    if ~isempty(value)
        snapshot.(name) = value;
    end
end
numericNames = {'schemaVersion','rearArbStiffness_NmPerRad', ...
    'frontSpringRate_lb_in','rearSpringRate_lb_in','frontRideHeight_in', ...
    'rearRideHeight_in','driverWeight_kg', ...
    'rearWeightDistribution_percent'};
for i = 1:numel(numericNames)
    name = numericNames{i};
    [value,found] = numericField(setupSpec,name,true);
    if found
        snapshot.(name) = value;
    end
end
[isBaseline,found] = logicalField(setupSpec,'isBaseline');
if found
    snapshot.isBaseline = isBaseline;
end
if isfield(setupSpec,'assetFingerprints') && ...
        isstruct(setupSpec.assetFingerprints) && ...
        isscalar(setupSpec.assetFingerprints)
    fingerprints = struct();
    fingerprintNames = {'aeroMap','camberRatios','camberModels'};
    for i = 1:numel(fingerprintNames)
        name = fingerprintNames{i};
        value = textField(setupSpec.assetFingerprints,name);
        if ~isempty(value)
            fingerprints.(name) = value;
        end
    end
    if ~isempty(fieldnames(fingerprints))
        snapshot.assetFingerprints = fingerprints;
    end
end
end

function vehicle = vehicleSnapshot(car)
vehicle = struct();
propertyMap = {
    'mass_kg','M'
    'wheelbase_m','W_b'
    'cgToFrontAxle_m','l_f'
    'cgToRearAxle_m','l_r'
    'frontTrack_m','t_f'
    'rearTrack_m','t_r'
    'wheelRadius_m','R'
    'cgHeight_m','h_g'
    'frontRollCenterHeight_m','h_rf'
    'rearRollCenterHeight_m','h_rr'
    'yawInertia_kg_m2','I_zz'
    'gravity_mps2','g'
    'rollingResistanceCoefficient','Crr'
    'wheelInertia_kg_m2','I_wheel'
    'drivelineInertia_kg_m2','I_driveline'};
for i = 1:size(propertyMap,1)
    [value,found] = numericMember(car,propertyMap{i,2},true);
    if found
        vehicle.(propertyMap{i,1}) = value;
    end
end
end

function driver = driverSnapshot(setupSpec,vehicle)
driver = struct();
[weight,hasWeight] = numericField(setupSpec,'driverWeight_kg',true);
if hasWeight
    % Keep the setup field name and the compact-model convenience name.
    driver.driverWeight_kg = weight;
    driver.weight_kg = weight;
end
[rearDistribution,hasDistribution] = numericField( ...
    setupSpec,'rearWeightDistribution_percent',true);
if hasDistribution
    driver.rearWeightDistribution_percent = rearDistribution;
    driver.frontWeightDistribution_percent = 100 - rearDistribution;
end
if hasWeight && isfield(vehicle,'mass_kg')
    driver.vehicleMass_kg = vehicle.mass_kg;
    driver.massWithoutDriver_kg = vehicle.mass_kg - weight;
end
end

function suspension = suspensionSnapshot(car,setupSpec)
suspension = struct();
setupNames = {'rearArbStiffness_NmPerRad','frontSpringRate_lb_in', ...
    'rearSpringRate_lb_in','frontRideHeight_in','rearRideHeight_in'};
for i = 1:numel(setupNames)
    name = setupNames{i};
    [value,found] = numericField(setupSpec,name,true);
    if found
        suspension.(name) = value;
    end
end

rideHeightAero = memberValue(car,'rideHeightAero');
rideNames = {'enabled','spring_rate_front_lb_in', ...
    'spring_rate_rear_lb_in','motion_ratio_front','motion_ratio_rear', ...
    'wheel_rate_front_Npm','wheel_rate_rear_Npm', ...
    'static_front_ride_height_in','static_rear_ride_height_in', ...
    'map_reference_front_ride_height_in','map_reference_rear_ride_height_in'};
rideSnapshot = copyPrimitiveFields(rideHeightAero,rideNames);
if ~isempty(fieldnames(rideSnapshot))
    suspension.rideHeightAero = rideSnapshot;
end

[frontArb,hasFrontArb] = numericMember(car,'k_rf',true);
[rearArb,hasRearArb] = numericMember(car,'k_rr',true);
if hasFrontArb
    suspension.frontArbStiffness_NmPerRad = frontArb;
end
if hasRearArb && ~isfield(suspension,'rearArbStiffness_NmPerRad')
    suspension.rearArbStiffness_NmPerRad = rearArb;
end

[rollSplit,hasRollSplit] = numericMember(car,'R_sf',true);
if hasRollSplit
    suspension.R_sf = rollSplit;
end
[wheelRateFront,hasWheelRateFront] = numericMember( ...
    rideHeightAero,'wheel_rate_front_Npm',true);
[wheelRateRear,hasWheelRateRear] = numericMember( ...
    rideHeightAero,'wheel_rate_rear_Npm',true);
if hasWheelRateFront
    suspension.wheelRateFront_Npm = wheelRateFront;
end
if hasWheelRateRear
    suspension.wheelRateRear_Npm = wheelRateRear;
end

[frontSpring,hasFrontSpring] = numericField( ...
    setupSpec,'frontSpringRate_lb_in',true);
[rearSpring,hasRearSpring] = numericField( ...
    setupSpec,'rearSpringRate_lb_in',true);
[motionRatioFront,hasMotionRatioFront] = numericMember( ...
    rideHeightAero,'motion_ratio_front',true);
[motionRatioRear,hasMotionRatioRear] = numericMember( ...
    rideHeightAero,'motion_ratio_rear',true);
[track,hasTrack] = numericMember(car,'t_f',true);
if hasFrontSpring && hasRearSpring && hasMotionRatioFront && ...
        hasMotionRatioRear && hasTrack && hasFrontArb && hasRearArb
    try
        derived = rampSpeed.deriveRollStiffness(frontSpring,rearSpring, ...
            motionRatioFront,motionRatioRear,track,frontArb,rearArb);
        derivedNames = fieldnames(derived);
        for i = 1:numel(derivedNames)
            suspension.(derivedNames{i}) = double(derived.(derivedNames{i}));
        end
        % Preserve the exact roll split used by the generated Car.
        if hasRollSplit
            suspension.R_sf = rollSplit;
        end
    catch
        % The setup and raw Car scalar remain useful if legacy optional
        % suspension data is incomplete or cannot be re-derived.
    end
end
end

function aero = aeroSnapshot(car,setupSpec)
aero = struct();
mapId = textField(setupSpec,'aeroMapId');
if ~isempty(mapId)
    aero.aeroMapId = mapId;
end

catalogLabel = '';
relativePath = '';
mapPath = '';
if ~isempty(mapId)
    try
        catalog = rampSpeed.aeroMapCatalog();
        for i = 1:numel(catalog)
            if strcmp(textValue(catalog(i).id),mapId)
                catalogLabel = textValue(catalog(i).label);
                relativePath = textValue(catalog(i).relativePath);
                mapPath = textValue(catalog(i).path);
                break
            end
        end
    catch
        % The setup's map ID is still retained when a catalog is unavailable.
    end
end
if ~isempty(catalogLabel)
    aero.aeroMapLabel = catalogLabel;
end
if ~isempty(relativePath)
    aero.aeroMapRelativePath = relativePath;
end

aeroObject = memberValue(car,'aero');
aeroNames = {'cda','cla','D_f','D_r','cla_p_deg_p','D_p_deg_p','rho'};
aeroSnapshotValues = copyPrimitiveFields(aeroObject,aeroNames);
aeroFields = fieldnames(aeroSnapshotValues);
for i = 1:numel(aeroFields)
    aero.(aeroFields{i}) = aeroSnapshotValues.(aeroFields{i});
end

mapObject = memberValue(aeroObject,'map');
if isobject(mapObject) && isscalar(mapObject)
    sourcePath = textMember(mapObject,'sourcePath');
    if isempty(mapPath)
        mapPath = sourcePath;
    end
    [frontRange,hasFrontRange] = numericMember(mapObject,'frontRangeIn',false);
    [rearRange,hasRearRange] = numericMember(mapObject,'rearRangeIn',false);
    if hasFrontRange && numel(frontRange) == 2
        aero.aeroMapFrontRange_in = frontRange(:).';
    end
    if hasRearRange && numel(rearRange) == 2
        aero.aeroMapRearRange_in = rearRange(:).';
    end
    interpolationMode = textMember(mapObject,'interpolationMode');
    if ~isempty(interpolationMode)
        aero.aeroMapInterpolationMode = interpolationMode;
    end
    for name = {'claScale','cdaScale','copOffset'}
        [value,found] = numericMember(mapObject,name{1},true);
        if found
            aero.(['aeroMap' upper(name{1}(1)) name{1}(2:end)]) = value;
        end
    end
end

fingerprint = textField(setupSpec,'aeroMapFingerprint');
if isempty(fingerprint) && isfield(setupSpec,'assetFingerprints') && ...
        isstruct(setupSpec.assetFingerprints) && ...
        isscalar(setupSpec.assetFingerprints)
    fingerprint = textField(setupSpec.assetFingerprints,'aeroMap');
end
if isempty(fingerprint) && ~isempty(mapPath)
    try
        fingerprint = textValue(rampSpeed.fingerprintFile(mapPath));
    catch
        % Fingerprinting is optional; retain map identity and path provenance.
    end
end
if ~isempty(fingerprint)
    aero.aeroMapFingerprint = fingerprint;
end

rideHeightAero = memberValue(car,'rideHeightAero');
[enabled,hasEnabled] = logicalMember(rideHeightAero,'enabled');
if hasEnabled
    aero.rideHeightAeroEnabled = enabled;
end
end

function powertrain = powertrainSnapshot(car,powertrainModel)
powertrain = struct('powertrainModel',powertrainModel);
powertrainObject = memberValue(car,'powertrain');
names = {'redline','shift_point','gears','primary_reduction', ...
    'torque_fn','shift_time','final_drive','wheel_radius', ...
    'drivetrain_efficiency','brake_distribution','G_d1', ...
    'G_d2_overrun','G_d2_driving','max_braking_torque', ...
    'switch_gear_velocities'};
for i = 1:numel(names)
    [value,found] = numericMember(powertrainObject,names{i},false);
    if found
        powertrain.(names{i}) = value;
    end
end
if isfield(powertrain,'gears') && isfield(powertrain,'final_drive') && ...
        isfield(powertrain,'primary_reduction')
    powertrain.totalGearReductions = powertrain.gears .* ...
        powertrain.final_drive .* powertrain.primary_reduction;
end
end

function output = copyPrimitiveFields(source,names)
output = struct();
for i = 1:numel(names)
    name = names{i};
    [value,found] = numericMember(source,name,false);
    if found
        output.(name) = value;
        continue
    end
    [value,found] = logicalMember(source,name);
    if found
        output.(name) = value;
    end
end
end

function value = memberValue(owner,name)
value = [];
if isstruct(owner) && isscalar(owner) && isfield(owner,name)
    value = owner.(name);
elseif isobject(owner) && isscalar(owner) && isprop(owner,name)
    try
        value = owner.(name);
    catch
        value = [];
    end
end
end

function [value,found] = numericField(source,name,scalarOnly)
[value,found] = numericMember(source,name,scalarOnly);
end

function [value,found] = numericMember(owner,name,scalarOnly)
value = [];
found = false;
raw = memberValue(owner,name);
if isempty(raw) || (~isnumeric(raw) && ~islogical(raw)) || ~isreal(raw)
    return
end
raw = double(raw);
if any(~isfinite(raw(:))) || (scalarOnly && ~isscalar(raw))
    return
end
value = raw;
found = true;
end

function [value,found] = logicalField(source,name)
[value,found] = logicalMember(source,name);
end

function [value,found] = logicalMember(owner,name)
value = false;
found = false;
raw = memberValue(owner,name);
if (islogical(raw) || isnumeric(raw)) && isreal(raw) && ...
        isscalar(raw) && isfinite(double(raw))
    value = logical(raw);
    found = true;
end
end

function value = textField(source,name)
value = '';
if isstruct(source) && isscalar(source) && isfield(source,name)
    value = textValue(source.(name));
end
end

function value = textMember(owner,name)
value = textValue(memberValue(owner,name));
end

function value = textValue(raw)
value = '';
if ischar(raw) && (isrow(raw) || isempty(raw))
    value = raw;
elseif isstring(raw) && isscalar(raw) && ~ismissing(raw)
    value = char(raw);
end
end

function output = toDataOnly(value)
if isstruct(value)
    output = value;
    names = fieldnames(value);
    for element = 1:numel(value)
        for i = 1:numel(names)
            name = names{i};
            output(element).(name) = toDataOnly(value(element).(name));
        end
    end
elseif iscell(value)
    output = cell(size(value));
    for i = 1:numel(value)
        output{i} = toDataOnly(value{i});
    end
elseif isstring(value)
    if isscalar(value)
        if ismissing(value)
            output = '';
        else
            output = char(value);
        end
    else
        output = cellstr(value);
    end
elseif isnumeric(value) || islogical(value) || ischar(value)
    output = value;
else
    error('rampSpeed:nonSerializableModelField', ...
        'Model metadata cannot contain values of class %s.',class(value));
end
end

function assertDataOnly(value,path)
if isstruct(value)
    names = fieldnames(value);
    for element = 1:numel(value)
        for i = 1:numel(names)
            assertDataOnly(value(element).(names{i}), ...
                [path '.' names{i}]);
        end
    end
elseif iscell(value)
    for i = 1:numel(value)
        assertDataOnly(value{i},sprintf('%s{%d}',path,i));
    end
elseif isnumeric(value) || islogical(value) || ischar(value)
    return
else
    error('rampSpeed:nonSerializableModelField', ...
        'The %s field has unsupported class %s.',path,class(value));
end
end
