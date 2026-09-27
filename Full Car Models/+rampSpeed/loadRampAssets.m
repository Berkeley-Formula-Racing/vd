function assets = loadRampAssets(config,setup)
%LOADRAMPASSETS Resolve and load the immutable assets used by one setup.
validateattributes(config,{'struct'},{'scalar'},mfilename,'config');
validateattributes(setup,{'struct'},{'scalar'},mfilename,'setup');
if ~isfield(setup,'aeroMapId')
    error('rampSpeed:invalidSetup','Setup must contain aeroMapId.');
end

root = fileparts(fileparts(mfilename('fullpath')));
componentRoot = fullfile(root,'carComponents');
catalog = rampSpeed.aeroMapCatalog();
ids = string({catalog.id});
index = find(ids == string(setup.aeroMapId),1);
if isempty(index)
    error('rampSpeed:aeroMapNotFound', ...
        'Aero-map ID %s is not present in rampSpeed.aeroMapCatalog.', ...
        string(setup.aeroMapId));
end

mapPath = string(catalog(index).path);
camberRatioPath = string(fullfile(componentRoot,'camberratiossmoothed.mat'));
camberModelPath = string(fullfile(componentRoot,'camber_models_fast.mat'));
map = rampSpeed.validateAeroMapFile(mapPath);
ratioData = load(char(camberRatioPath), ...
    'camber2indices','camber4indices','camber2ratio','camber4ratio');
modelData = load(char(camberModelPath),'betaL','betaR');

assets = struct();
assets.projectRoot = string(root);
assets.map = map;
assets.mapPath = mapPath;
assets.mapId = string(catalog(index).id);
assets.mapLabel = string(catalog(index).label);
assets.mapRelativePath = string(catalog(index).relativePath);
assets.camberRatioPath = camberRatioPath;
assets.camberModelPath = camberModelPath;
assets.camberRatios = ratioData;
assets.camberModelData = modelData;
assets.fingerprints = struct( ...
    'aeroMap',rampSpeed.fingerprintFile(mapPath), ...
    'camberRatios',rampSpeed.fingerprintFile(camberRatioPath), ...
    'camberModels',rampSpeed.fingerprintFile(camberModelPath));
end
