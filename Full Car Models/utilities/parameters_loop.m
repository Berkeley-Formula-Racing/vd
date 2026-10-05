function [car_cell,sampleTable] = parameters_loop(cP,aP,eP,DTp,Bp,tP,samplingType,numSamples)
%PARAMETERS_LOOP Preserve the legacy cell adapter for QSS case generation.
%   Sampling and car construction use the named setup core. The public
%   output remains an N-by-2 cell array with lap then acceleration cars.

if nargin < 7 || isempty(samplingType), samplingType = "FullFactorial"; end
if nargin < 8, numSamples = []; end

baseSetup = vehicleConfig.fromLegacyGroups(cP,aP,eP,DTp,Bp,tP);
[setups,sampleTable] = ...
    vehicleConfig.sampleSetups(baseSetup,samplingType,numSamples);

% Keep the direct legacy wrapper's public error identifier and check order.
if baseSetup.aeroMap.enabled && ~readableMap(baseSetup.aeroMap.path)
    error('parameters_loop:aeroMapNotFound', ...
        'use_aeromap is true but aeroParams.map_path does not name a readable CSV file.')
end

cases = vehicleConfig.buildCases(setups);
car_cell = vehicleConfig.toLegacyCells(cases);
end

function tf = readableMap(path)
tf = false;
if isempty(path) || ~(ischar(path) || (isstring(path) && isscalar(path)))
    return
end
if ~isfile(path)
    return
end
[fid,~] = fopen(char(path),'r');
if fid < 0
    return
end
fclose(fid);
tf = true;
end
