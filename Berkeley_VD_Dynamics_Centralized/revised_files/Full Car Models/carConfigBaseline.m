function [car,config] = carConfigBaseline()
%CARCONFIGBASELINE Return the stable single-car ramp-speed baseline.
%   [CAR,CONFIG] = CARCONFIGBASELINE creates exactly one acceleration-oriented
%   Car object. CONFIG is the frozen, serializable source for future setup
%   duplication and editing; it reads the same master defaults as carConfig.

root = fileparts(mfilename('fullpath'));
addpath(fullfile(root,'config'));
config = vdRampConfig();
defaultSetup = config.defaultSetup;
assets = rampSpeed.loadRampAssets(config,defaultSetup);
[car,~] = buildSingleRampCar(config.fixedParameters,defaultSetup, ...
    assets.map,assets);
end
