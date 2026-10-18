function cases = buildCases(setups,assets)
%BUILDCASES Build named lap and acceleration Cars for sampled setups.
if nargin < 2
    assets = [];
end
if ~isstruct(setups)
    error('vehicleConfig:badSchema', ...
        'setups must be a struct array produced from vehicle configurations.');
end

template = struct('id',"",'label',"",'designRow',NaN, ...
    'setup',struct(),'lapCar',[],'accelCar',[]);
cases = repmat(template,numel(setups),1);

cachedIdentities = cell(0,1);
cachedAssets = cell(0,1);
if ~isempty(assets)
    if ~isstruct(assets) || ~isscalar(assets) || ~isfield(assets,'identity')
        error('vehicleConfig:assetMismatch', ...
            'A shared asset context must include its source identity.');
    end
    cachedIdentities{end+1,1} = assets.identity;
    cachedAssets{end+1,1} = assets;
end

for i = 1:numel(setups)
    identity = setupAssetIdentity(setups(i));
    cacheIndex = find(cellfun( ...
        @(candidate) isequaln(candidate,identity),cachedIdentities),1);
    if isempty(cacheIndex)
        context = vehicleConfig.loadAssets(setups(i));
        cachedIdentities{end+1,1} = identity;
        cachedAssets{end+1,1} = context;
    else
        context = cachedAssets{cacheIndex};
    end

    cases(i).id = "case-" + compose('%06d',i);
    cases(i).label = setups(i).label;
    cases(i).designRow = i;
    cases(i).setup = setups(i);
    cases(i).lapCar = vehicleConfig.buildCar(setups(i),"lap",context);
    cases(i).accelCar = vehicleConfig.buildCar(setups(i),"accel",context);
end
end
