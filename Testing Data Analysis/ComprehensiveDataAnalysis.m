clear
close all
clc
set(0,'DefaultTextInterpreter','none')

%% Import Data
filename = "C:\Users\johny\Downloads\Logs\0045.csv";
opt = detectImportOptions(filename);
opt.VariableUnitsLine = 16;
T = readtable(filename, opt);
T = T(2:end,:); % remove first row

%% Inputs
% time selection
timeRange = [0,Inf];

% plot selection
overviewPlot = 0;
bodyMovementPlot = 0;
tireTemperaturePlot = 0;
damperVelocityHistogram = 0;
aeroPlot = 0;
powerPlot = 0;

%% Car Data
C = carConfig();
car = C{1,1};
car.MR_F = 0.74;
car.MR_R = 1;
car.k = 200*4.45*39.37; % N/m

%% ADL vs Telemetry Unit
if ~any(strcmp(T.Properties.VariableNames,'STEERINGANGLE'))
    T.WheelSpdFL = T.WheelSpeedFrontLeft;
    T.WheelSpdFR = T.WheelSpeedFrontRight;
    T.WheelSpdRL = T.WheelSpeedRearLeft;
    T.WheelSpdRR = T.WheelSpeedRearRight;
    T.BrakePres_F = 0*T.Time;
    T.BrakePres_R = 0*T.Time;
    T.TPS = 0*T.Time;
    T.RPM = T.EngineSpeed;
end

%% Filter Data
T.STEERINGANGLE = T.STEERINGANGLE * 20/45;
variablesToFilter = {'SHOCKFL','SHOCKFR','SHOCKRL','SHOCKRR','STEERINGANGLE', 'WheelSpeedRearLeft','WheelSpeedRearRight','WheelSpeedFrontLeft','WheelSpeedFrontRight'};
meanRangeSeconds = 0.25;
meanTimestep = mean(diff(T.Time));
order = meanRangeSeconds/meanTimestep;
T(:,variablesToFilter) = array2table(movmean(table2array(T(:,variablesToFilter)),order,1));

%% Zero Sensors
t0 = 10;
T.SHOCKFL = T.SHOCKFL - mean(T.SHOCKFL(T.Time<t0));
T.SHOCKFR = T.SHOCKFR - mean(T.SHOCKFR(T.Time<t0));
T.SHOCKRL = T.SHOCKRL - mean(T.SHOCKRL(T.Time<t0));
T.SHOCKRR = T.SHOCKRR - mean(T.SHOCKRR(T.Time<t0));

T.wheelPosFL = T.SHOCKFL./car.MR_F;
T.wheelPosFR = T.SHOCKFR./car.MR_F;
T.wheelPosRL = T.SHOCKRL./car.MR_R;
T.wheelPosRR = T.SHOCKRR./car.MR_R;
T.Properties.VariableUnits(strcmp(T.Properties.VariableNames, 'wheelPosFL')) = {'mm'};
T.Properties.VariableUnits(strcmp(T.Properties.VariableNames, 'wheelPosFR')) = {'mm'};
T.Properties.VariableUnits(strcmp(T.Properties.VariableNames, 'wheelPosRL')) = {'mm'};
T.Properties.VariableUnits(strcmp(T.Properties.VariableNames, 'wheelPosRR')) = {'mm'};

%% Calculate New Channels
T.SuspHeave = (T.wheelPosFL + T.wheelPosFR + T.wheelPosRL + T.wheelPosRR)/4;
T.FrontHeave = (T.wheelPosFL + T.wheelPosFR)/2;
T.RearHeave = (T.wheelPosRL + T.wheelPosRR)/2;
T.SuspPitch = rad2deg((T.RearHeave-T.FrontHeave)/(car.W_b*1000));
T.SuspRoll = rad2deg((T.wheelPosFL + T.wheelPosRL - T.wheelPosFR - T.wheelPosRR)/(2*car.t_f*1000));
T = setUnits(T, {'SuspHeave', 'FrontHeave', 'RearHeave'}, 'mm');
T = setUnits(T, {'SuspPitch', 'SuspRoll'},'deg');

T.damperVelFL = [0; diff(T.wheelPosFL)./diff(T.Time)];
T.damperVelFR = [0; diff(T.wheelPosFR)./diff(T.Time)];
T.damperVelRL = [0; diff(T.wheelPosRL)./diff(T.Time)];
T.damperVelRR = [0; diff(T.wheelPosRR)./diff(T.Time)];
T = setUnits(T, {'damperVelFL', 'damperVelFR', 'damperVelRL', 'damperVelRR'}, 'mm/s');

T.Speed = (T.WheelSpdFL + T.WheelSpdFR + T.WheelSpdRL + T.WheelSpdRR)/4 * 1000/3600;
T = setUnits(T, 'Speed', 'm/s');

T.TheoreticalDrag = (car.aero.cda * car.aero.rho/2)* T.Speed.^2;
T.TheoreticalLift = (car.aero.cla * car.aero.rho/2) * T.Speed.^2;
T = setUnits(T, {'TheoreticalDrag', 'TheoreticalLift'}, 'N');
