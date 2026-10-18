function [output,eventParams,designTable,baselineTable] = carConfig( ...
    samplingType,numSamples,tireModel,uncertaintyMode,outputFormat)
%CARCONFIG Build legacy car cells or named QSS cases.
%   By default, output is an N-by-2 cell array with lap cars in column 1
%   and acceleration cars in column 2. Pass "cases" as the fifth input to
%   return an N-by-1 named case struct array instead.
%
%   [cars,eventParams,X,B] = carConfig();
%   [cars,eventParams,X] = carConfig("Explicit",overrideTable);
%   [cases,eventParams,X,B] = carConfig( ...
%       "Explicit",overrideTable,"legacy","nominal","cases");

if nargin < 1 || isempty(samplingType), samplingType = "FullFactorial"; end
if nargin < 2, numSamples = []; end
if nargin < 3 || isempty(tireModel), tireModel = "legacy"; end
if nargin < 4 || isempty(uncertaintyMode), uncertaintyMode = "nominal"; end
if nargin < 5 || isempty(outputFormat), outputFormat = "cells"; end
outputFormat = validateOutputFormat(outputFormat);

[baseSetup,eventParams] = vehicleConfig.defaultSetup(tireModel,uncertaintyMode);
[baselineSetups,baselineTable] = ...
    vehicleConfig.sampleSetups(baseSetup,"FullFactorial",[]);

if baseSetup.aeroMap.enabled && ~readableMap(baseSetup.aeroMap.path)
    error('parameters_loop:aeroMapNotFound', ...
        'use_aeromap is true but aeroParams.map_path does not name a readable CSV file.')
end

request = lower(string(samplingType));
switch request
    case {"fullfactorial","grid"}
        setups = baselineSetups;
        designTable = baselineTable;

    case "explicit"
        [setups,designTable] = ...
            vehicleConfig.sampleSetups(baseSetup,"Explicit",numSamples);

    case {"lhs","random"}
        warning('carConfig:legacyDOE', ...
            ['carConfig legacy DOE sampling is deprecated; use Explicit ' ...
             'with a resolved DOE design table.']);
        study = DOEStudyConfig();
        resolved = doeResolveStudy(study,baselineTable);
        [~,overrides] = ...
            doeInitialDesign(resolved,numSamples,study.randomSeed,request);
        [setups,designTable] = ...
            vehicleConfig.sampleSetups(baseSetup,"Explicit",overrides);

    otherwise
        error('carConfig:badRequest', ...
            'samplingType must be FullFactorial, Explicit, LHS, or Random.')
end

cases = vehicleConfig.buildCases(setups);
if outputFormat == "cases"
    output = cases;
else
    output = vehicleConfig.toLegacyCells(cases);
end
end

function outputFormat = validateOutputFormat(outputFormat)
if ~(ischar(outputFormat) || ...
        (isstring(outputFormat) && isscalar(outputFormat))) || ...
        (ischar(outputFormat) && ~isrow(outputFormat))
    error('carConfig:badOutputFormat', ...
        'outputFormat must be scalar text "cells" or "cases".');
end
outputFormat = lower(string(outputFormat));
if ~any(outputFormat == ["cells","cases"])
    error('carConfig:badOutputFormat', ...
        'outputFormat must be "cells" or "cases".');
end
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
