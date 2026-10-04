function assets = loadAssets(setup)
%LOADASSETS Load immutable QSS map, camber, and optional LC0 tire assets once.
% Call once for a batch and pass the returned context to each buildCar call.
identity = setupAssetIdentity(setup);
if identity.profile ~= "qss"
    error('vehicleConfig:badSchema', ...
        'loadAssets currently supports the QSS profile only.');
end

map = [];
mapFingerprint = "";
if identity.aeroMapEnabled
    if strlength(identity.aeroMapPath) == 0 || ~isfile(identity.aeroMapPath)
        error('vehicleConfig:assetNotFound', ...
            'Enabled QSS aero map does not exist: %s',identity.aeroMapPath);
    end
    map = AeroMap(char(identity.aeroMapPath));
    map = map.withCorrections(identity.claScale,identity.cdaScale,identity.copOffset);
    mapFingerprint = rampSpeed.fingerprintFile(identity.aeroMapPath);
end

if ~isfile(identity.camberModelPath)
    error('vehicleConfig:assetNotFound', ...
        'Camber model data does not exist: %s',identity.camberModelPath);
end
if ~isfile(identity.camberRatioPath)
    error('vehicleConfig:assetNotFound', ...
        'Camber ratio data does not exist: %s',identity.camberRatioPath);
end
camberModelData = load(char(identity.camberModelPath),'betaL','betaR');
camberRatioData = load(char(identity.camberRatioPath), ...
    'camber2indices','camber4indices','camber2ratio','camber4ratio');
if ~all(isfield(camberModelData,{'betaL','betaR'})) || ...
        ~all(isfield(camberRatioData, ...
            {'camber2indices','camber4indices','camber2ratio','camber4ratio'}))
    error('vehicleConfig:badAsset', ...
        'Camber assets are missing fields required by the QSS constructors.');
end

model = [];
tireFingerprint = "";
if strlength(identity.tireModelPath) > 0
    if ~isfile(identity.tireModelPath)
        error('buildTireModel:missingArtifact', ...
            'Nondimensional tire artifact is missing: %s',identity.tireModelPath);
    end
    loaded = load(char(identity.tireModelPath),'model');
    if ~isfield(loaded,'model')
        error('buildTireModel:badArtifact', ...
            'Nondimensional artifact must contain a variable named model.');
    end
    model = loaded.model;
    tireFingerprint = rampSpeed.fingerprintFile(identity.tireModelPath);
end

assets = struct();
assets.projectRoot = string(fileparts(fileparts(char(identity.camberModelPath))));
assets.map = map;
assets.mapPath = identity.aeroMapPath;
assets.camberModelPath = identity.camberModelPath;
assets.camberRatioPath = identity.camberRatioPath;
assets.camberModelData = camberModelData;
assets.camberRatioData = camberRatioData;
assets.tireModelPath = identity.tireModelPath;
assets.tireModel = model;
assets.identity = identity;
assets.fingerprints = struct( ...
    'aeroMap',mapFingerprint, ...
    'camberRatios',rampSpeed.fingerprintFile(identity.camberRatioPath), ...
    'camberModels',rampSpeed.fingerprintFile(identity.camberModelPath), ...
    'tireModel',tireFingerprint);
end
