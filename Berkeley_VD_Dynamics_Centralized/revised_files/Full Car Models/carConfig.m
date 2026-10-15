function [carCell,eventParams,designTable,baselineTable] = carConfig(samplingType,numSamples,tireModel,uncertaintyMode)
% outputs carCell:     one car per grid combination or DOE sample,
%                      column 1 the lap car, column 2 the accel car
%        eventParams:  everything the dynamic events need that is NOT a
%                      property of the car -- course geometry, lap counts,
%                      driver calibration and the scoring reference times.
%                      Events2 requires it; ask for it whenever you build one:
%                          [carCell,eventParams] = carConfig();
%                          ev = Events2(carCell{1,1},carCell{1,2},eventParams);
%                      Asking for only carCell still works, so nothing that
%                      does not build an Events2 needs changing.
%        designTable:  exact sampled parameter values, one row per car.
%        baselineTable: calibrated scalar parameter values, one row.
%
% Default is the calibrated baseline/full-factorial configuration. For DOE:
%   [cars,eventParams,X] = carConfig("Explicit",overrideTable);
%   [cars,eventParams,X] = carConfig("FullFactorial",[],"lc0_nd");
%   [cars,eventParams,X] = carConfig("FullFactorial",[],"lc0_nd","low");
if nargin < 1 || isempty(samplingType), samplingType = "FullFactorial"; end
if nargin < 2, numSamples = []; end
if nargin < 3 || isempty(tireModel), tireModel = "legacy"; end
if nargin < 4 || isempty(uncertaintyMode), uncertaintyMode = "nominal"; end

% Shared defaults. Study overrides still belong in the sampling table.
addpath(fullfile(fileparts(mfilename('fullpath')),'config'));
cfg = vdVehicleConfig();
carParams = cfg.car;
aeroParams = cfg.aero;
eParams = cfg.engine;
eParams.torque_fn = eParams.torque_fn();
DTparams = cfg.drivetrain;
Bparams = cfg.brakes;
tireParams = vdTireConfig(tireModel,uncertaintyMode);
tireParams.model_type = lower(string(tireModel));
tireParams.model_uncertainty = lower(string(uncertaintyMode));
eventParams = cfg.events;

%% calibrated scalar baseline and requested cases
[baselineCars,baselineTable] = parameters_loop(carParams,aeroParams,eParams, ...
    DTparams,Bparams,tireParams,"FullFactorial",[]);

request = lower(string(samplingType));
switch request
    case {"fullfactorial","grid"}
        carCell = baselineCars;
        designTable = baselineTable;

    case "explicit"
        [carCell,designTable] = parameters_loop(carParams,aeroParams,eParams, ...
            DTparams,Bparams,tireParams,"Explicit",numSamples);

    case {"lhs","random"}
        warning('carConfig:legacyDOE', ...
            ['carConfig legacy DOE sampling is deprecated; use Explicit ' ...
             'with a resolved DOE design table.']);
        study = DOEStudyConfig();
        resolved = doeResolveStudy(study,baselineTable);
        [~,overrides] = doeInitialDesign(resolved,numSamples,study.randomSeed,request);
        [carCell,designTable] = parameters_loop(carParams,aeroParams,eParams, ...
            DTparams,Bparams,tireParams,"Explicit",overrides);

    otherwise
        error('carConfig:badRequest', ...
            'samplingType must be FullFactorial, Explicit, LHS, or Random.')
end
