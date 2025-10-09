clear
close all
clc
set(0,'DefaultTextInterpreter','none')
%% Import Data

% filename = '20220402-0011602.csv';
% filename = 'braketests-alameda-09.csv';

filename = 'autocross_3.csv';

%filename = 'straighlinepullsm400.csv';

% filename = 'full_skidpad1.csv';
% filename = 'full_skidpad2.csv';
% filename = 'skidpad2_50hz.csv'; % shows accel granularity - also good tire temp

%filename = 'skippad_alex_5runs_with_start.csv';

% filename = 'damper_testing_rebounds_50hz.csv';
% filename = 'AllDay_DamperTuning.csv';

% include units
opt = detectImportOptions(filename);
opt.VariableUnitsLine = 16;
T = readtable(filename, opt);

T = T(2:end,:); % remove first row

%% Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% time selection
timeRange = [0,Inf]; %time range to be plotted

% plot selection
overviewPlot = 0;
lateralPlot = 0;
longitudinalPlot = 0;
bodyMovementPlot = 0;
tireTemperaturePlot = 0;
damperVelocityHistogram = 0;

understeerPlot = 0;
understeerPlotRadiusDependent = 0;

GGPlot = 0;
GGVPlot = 0;

aeroPlot = 0;

powerPlot = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Car Data

C = carConfig();
car = C{1,1};

car.MR_F = 0.74;
car.MR_R = 1;
car.k = 200*4.45*39.37; % N/m ---add k front and k rear

%% ADL vs Telemetry Unit

if ~any(strcmp(T.Properties.VariableNames,'StAngle')) % if data is from ADL
    %T.StAngle = T.SteeringAngle;

    T.WheelSpdFL = T.GroundSpeedLeft; % Km/h (?)
    T.WheelSpdFR = T.GroundSpeedRight;
    T.WheelSpdRL = T.DriveSpeedLeft;
    T.WheelSpdRR = T.DriveSpeedRight;
    
    T.AccelX = 0*T.Time;
    T.AccelY = 0*T.Time;
    T.AccelZ = 0*T.Time;
    
    T.BrakePres_F = 0*T.Time;
    T.BrakePres_R = 0*T.Time;
    
    T.TPS = 0*T.Time;
    T.RPM = T.EngineRPM;
end

%% Filter Data

T.StAngle = T.StAngle * 20/45; % max steering angle is 20 deg, sensors read -45 -> 45 deg
T.AccelX = -T.AccelX; % long acceleration reverse (+ should be increasing forward speed)

% moving mean filter on selected data values
variablesToFilter = {'SuspPosFL','SuspPosFR','SuspPosRL','SuspPosRR',...
    'StAngle', 'WheelSpdRL','WheelSpdRR','WheelSpdFL','WheelSpdFR', 'AccelX', 'AccelY', 'AccelZ'};
meanRangeSeconds = 0.25; % length of averaging window (s)
meanTimestep = mean(diff(T.Time));
order = meanRangeSeconds/meanTimestep;
T(:,variablesToFilter) = ...
    array2table(movmean(table2array(T(:,variablesToFilter)),order,1));

%% Zero Sensors

% zero out wheel positions to their averages for the first 10 seconds
t0 = 10;
T.SuspPosFL = T.SuspPosFL - mean(T.SuspPosFL(T.Time<t0));
T.SuspPosFR = T.SuspPosFR - mean(T.SuspPosFR(T.Time<t0));
T.SuspPosRL = T.SuspPosRL - mean(T.SuspPosRL(T.Time<t0));
T.SuspPosRR = T.SuspPosRR - mean(T.SuspPosRR(T.Time<t0));

T.wheelPosFL = T.SuspPosFL./car.MR_F;
T.wheelPosFR = T.SuspPosFR./car.MR_F;
T.wheelPosRL = T.SuspPosRL./car.MR_R;
T.wheelPosRR = T.SuspPosRR./car.MR_R;
T = setUnits(T, {'wheelPosFL', 'wheelPosFR', 'wheelPosRL', 'wheelPosRR'},...
    getUnits(T, 'SuspPosRR'));

%% Calculate New Channels

% heave, roll, pitch, front heave, rear heave, using shock pot values
T.SuspHeave = (T.wheelPosFL + T.wheelPosFR + T.wheelPosRL + T.wheelPosRR)/4;
T.FrontHeave = (T.wheelPosFL + T.wheelPosFR)/2;
T.RearHeave = (T.wheelPosRL + T.wheelPosRR)/2;
T.SuspPitch = rad2deg((T.RearHeave-T.FrontHeave)/(car.W_b*1000));
T.SuspRoll = rad2deg((T.wheelPosFL + T.wheelPosRL - T.wheelPosFR - T.wheelPosRR)/(2*car.t_f*1000));
T = setUnits(T, {'SuspHeave', 'FrontHeave', 'RearHeave'},...
    getUnits(T, 'wheelPosFL'));
T = setUnits(T, {'SuspPitch', 'SuspRoll'},'deg');

% damperVelocities
T.damperVelFL = [0; diff(T.wheelPosFL)./diff(T.Time)];
T.damperVelFR = [0; diff(T.wheelPosFR)./diff(T.Time)];
T.damperVelRL = [0; diff(T.wheelPosRL)./diff(T.Time)];
T.damperVelRR = [0; diff(T.wheelPosRR)./diff(T.Time)];
T = setUnits(T, {'damperVelFL', 'damperVelFR', 'damperVelRL', 'damperVelRR'},...
    [getUnits(T, 'wheelPosRR') '/s']);

% speed
T = setUnits(T, {'WheelSpdFL', 'WheelSpdFR', 'WheelSpdRL', 'WheelSpdRR'}, 'km/h');
T.Speed = (T.WheelSpdFL + T.WheelSpdFR + T.WheelSpdRL + T.WheelSpdRR)/4 * 1000/3600; % km/h -> m/s
T = setUnits(T, 'Speed', 'm/s');
% corner radius channel using r= v^2/Ay
T.cornerRadius = T.Speed.^2./(T.AccelY*9.8);
T = setUnits(T, 'cornerRadius', 'm');

%% Use Aero parameters to calculate expected Drag and Lift

T.TheoreticalDrag = (car.aero.cda * car.aero.rho/2)* T.Speed.^2;
T.TheoreticalLift = (car.aero.cla * car.aero.rho/2) * T.Speed.^2;
T = setUnits(T, {'TheoreticalDrag', 'TheoreticalDrag'} , 'N');

%% time selection plots

if overviewPlot
    figure
    subplot(4,1,1)
    plotLine(T,timeRange,'WheelSpdFL')
    hold on
    plotLine(T,timeRange,'WheelSpdFR')
    plotLine(T,timeRange,'WheelSpdRL')
    plotLine(T,timeRange,'WheelSpdRR')
    legend('Interpreter','none')
    grid

    subplot(4,1,2)
    plotLine(T,timeRange,'BrakePres_F')
    hold on
    plotLine(T,timeRange,'BrakePres_R')
    legend('Interpreter','none')
    grid
    
    subplot(4,1,3)
    plotLine(T,timeRange,'TPS')
    legend('Interpreter','none')
    grid

    
    subplot(4,1,4)
    plotLine(T,timeRange,'StAngle')
    legend('Interpreter','none')
    grid
    
    sgtitle('Lockup Plot')
end

if lateralPlot
    figure
    subplot(5,1,1)
    plotLine(T,timeRange,'StAngle')
    legend('Interpreter','none')
    grid
    
    subplot(5,1,2)
    plotLine(T,timeRange,'AccelY')
    legend('Interpreter','none')
    grid
    
    subplot(5,1,3)
    plotLine(T,timeRange,'SuspRoll')
    legend('Interpreter','none')
    grid

    subplot(5,1,4)
    plotLine(T,timeRange,'Speed')
    yyaxis right
    plotLine(T,timeRange,'cornerRadius')
    legend('Interpreter','none')
    grid

    subplot(5,1,5)
    plotLine(T,timeRange,'AccelY')
    legend('Interpreter','none')
    grid

    sgtitle('Lateral Plot')
end

if longitudinalPlot
    figure
    subplot(5,1,1)
    plotLine(T,timeRange,'BrakePres_F')
    hold on
    plotLine(T,timeRange,'BrakePres_R')
    legend('Interpreter','none')
    grid
    
    subplot(5,1,2)
    plotLine(T,timeRange,'TPS')
    legend('Interpreter','none')
    grid
    
    subplot(5,1,3)
    plotLine(T,timeRange,'SuspPitch')
    legend('Interpreter','none')
    grid

    subplot(5,1,4)
    plotLine(T,timeRange,'Speed')
    yyaxis right
    ylim([-20,20])
    plotLine(T,timeRange,'RPM')
    legend('Interpreter','none')
    grid

    subplot(5,1,5)
    plotLine(T,timeRange,'AccelX')
    legend('Interpreter','none')
    grid

    sgtitle('Longitudinal Plot')
end

if bodyMovementPlot
    figure
    subplot(3,1,1)
    plotLine(T,timeRange,'SuspHeave')
    hold on
    yyaxis right
    plotLine(T,timeRange,'SuspPitch')
    plotLine(T,timeRange,'SuspRoll')
    legend('Interpreter','none')
    grid

    subplot(3,1,2)
    plotLine(T,timeRange,'FrontHeave')
    hold on
    plotLine(T,timeRange,'RearHeave')
    legend('Interpreter','none')
    grid

    subplot(3,1,3)
    scatter([1 2], [1 2])
    plotLine(T,timeRange,'wheelPosFL')
    hold on
    plotLine(T,timeRange,'wheelPosFR')
    plotLine(T,timeRange,'wheelPosRL')
    plotLine(T,timeRange,'wheelPosRR')
    legend('Interpreter','none')
    grid

    sgtitle('Body Movement')
end

if tireTemperaturePlot
    figure
    subplot(4,1,1)
    hold on
    plotLine(T,timeRange,'TireTmpLF2')
    plotLine(T,timeRange,'TireTmpLF4')
    plotLine(T,timeRange,'TireTmpLF6')
    legend('Interpreter','none')
    grid

    subplot(4,1,2)
    hold on
    plotLine(T,timeRange,'TireTmpRF2')
    plotLine(T,timeRange,'TireTmpRF4')
    plotLine(T,timeRange,'TireTmpRF6')
    legend('Interpreter','none')
    grid

    subplot(4,1,3)
    hold on
    plotLine(T,timeRange,'TireTmpLR2')
    plotLine(T,timeRange,'TireTmpLR4')
    plotLine(T,timeRange,'TireTmpLR6')
    legend('Interpreter','none')
    grid

    subplot(4,1,4)
    hold on
    plotLine(T,timeRange,'TireTmpRR2')
    plotLine(T,timeRange,'TireTmpRR4')
    plotLine(T,timeRange,'TireTmpRR6')
    legend('Interpreter','none')
    grid

    sgtitle('Tire Temperature')
end

%% damper velocity histogram
    
if damperVelocityHistogram   
    select = (T.Speed > 0.1);

    figure
    subplot(3,2,1)
    histogram(T.damperVelFL(select)), xlim([-100,100])
    title('damperVelFL'), xlabel('velocity (mm/s)')
    subplot(3,2,2)
    histogram(T.damperVelFR(select)), xlim([-100,100])
    title('damperVelFR'), xlabel('velocity (mm/s)')
    subplot(3,2,3)
    histogram(T.damperVelRL(select)), xlim([-100,100])
    title('damperVelRL'), xlabel('velocity (mm/s)')
    subplot(3,2,4)
    histogram(T.damperVelRR(select)), xlim([-100,100])
    title('damperVelRR'), xlabel('velocity (mm/s)')

    subplot(3,1,3);
    hold on
    plotLine(T, timeRange, 'damperVelFL')
    plotLine(T, timeRange, 'damperVelFR')
    plotLine(T, timeRange, 'damperVelRL')
    plotLine(T, timeRange, 'damperVelRR')
    title('Damper velocity vs time'), ylabel('velocity (mm/s)'), grid
end

%% Understeer gradient fitting
if understeerPlot
    cornerRadiusBounds = [5,7]; % m

    t01 = 20; % s
    minAy = 0.2; % g
    
    % moving mean filter on data
    meanRangeSeconds = 1; % s
    meanTimestep = mean(diff(T.Time));
    order = meanRangeSeconds/meanTimestep;
    
    delta = abs(movmean(T.StAngle, order));
    latAccel = abs(movmean(T.AccelY, order));
    cornerRadius = abs(movmean(T.cornerRadius, order));

    % select only some of data
    select1 = (T.Time>t01 & T.Time<(T.Time(end)-t01)) & (abs(latAccel) > minAy);

    delta = delta(select1);
    latAccel = latAccel(select1);
    cornerRadius = cornerRadius(select1);

    % exlcude data with corner radii out of bounds   
    exclude =  (cornerRadius < cornerRadiusBounds(1) | ...
        cornerRadius > cornerRadiusBounds(2));
    
    [f1, gof1] = fit(latAccel, delta,'poly1', 'Exclude', exclude);
    b1 = confint(f1);

    figure
    subplot(2,1,1)
    plot(f1, latAccel, delta, exclude),...
        title(['Steering Angle vs A_y | Understeer Gradient = '...
        num2str(f1.p1) ' (' num2str(b1(1,1)) ', ' num2str(b1(2,1)) ')'],...
        'Interpreter','tex')
    legend('Location','NorthWest')
    xlim([0 2]), ylim([-20, 20])
    xlabel('Lateral Acceleration (G)'), ylabel('Steering angle (deg)'), grid
    
    % corner radius histogram
    subplot(2,1,2)
    histogram(T.cornerRadius(abs(T.cornerRadius)<20),50)
    xline(cornerRadiusBounds), xline(-cornerRadiusBounds)
    xlabel(['Corner Radius (' getUnits(T, 'cornerRadius') ')'])
    title('Corner Radius')
end

%% Understeer gradient fitting (Corner Radius Dependent)
if understeerPlotRadiusDependent

    % disregard first and last <t01> seconds 
    t01 = 20; % s
    minAy = 0.2; % G
    
    % moving mean filter on data
    meanRangeSeconds = 1; % s
    meanTimestep = mean(diff(T.Time));
    order = meanRangeSeconds/meanTimestep;

    delta = abs(movmean(T.StAngle, order));
    latAccel = abs(movmean(T.AccelY, order));
    cornerRadius = abs(movmean(T.cornerRadius, order)); % average corner radius over 3x the time

    select1 = (T.Time>t01 & T.Time<(T.Time(end)-t01)) & (abs(latAccel) > minAy);

    delta = delta(select1);
    latAccel = latAccel(select1);
    cornerRadius = cornerRadius(select1);

    % break up data by corner radius, fit understeer gradient to each group
    numZones = 4;

    cornerRadiusBounds = linspace(3,30, numZones+1);
    fits = cell(2, 10);

    % 1-UG, 2- 95% interval, 3- 95% interval, 4- r^2
    understeerGradient = zeros(4,numZones); 

    figure
    for i = 1:numZones
        bounds = cornerRadiusBounds([i,(i+1)]);
        exclude =  (cornerRadius < bounds(1) | ...
            cornerRadius > bounds(2));
    
        [f1, gof1] = fit(latAccel, delta,'poly1', 'Exclude', exclude);
        b1 = confint(f1);

        understeerGradient(1,i) = f1.p1;
        understeerGradient(2,i) = b1(1,1);
        understeerGradient(3,i) = b1(2,1);
        understeerGradient(4,i) = gof1.rsquare;

        understeerGradient(5,i) = mean(bounds);
    
        subplot(3,numZones,i)
        scatter(latAccel(~exclude), delta(~exclude))
        hold on
        plot(f1)

        legend('off')
        xlabel('A_y (G)', 'Interpreter', 'tex'), ylabel('\delta (deg)', 'Interpreter', 'tex'), grid
        title([num2str(round(bounds(1))) ' < R < ' num2str(round(bounds(2))) ' UG = ' num2str(f1.p1)])

        ylim([-20, 20]);
        xlim([0,2]);
    end
    

    subplot(3,1,2)
    hold on
    errorbar(understeerGradient(5,:), understeerGradient(1,:), ...
        understeerGradient(2,:)-understeerGradient(1,:), '-s')

    ylabel('UG (deg/G)')
    xlabel(['Corner Radius (' getUnits(T, 'cornerRadius') ')'])

    title('Understeer Gradient vs Corner Radius (95% Confidence Bounds)')

    
    % corner radius histogram
    subplot(3,1,3)
    histogram(T.cornerRadius(abs(T.cornerRadius)<20),50)
    xline(cornerRadiusBounds), xline(-cornerRadiusBounds)
    xlabel(['Corner Radius (' getUnits(T, 'cornerRadius') ')'])
    title('Corner Radius')

    sgtitle('Understeer Gradient')
end

%% GG plot
if GGPlot
    figure
    scatter(T,'AccelY','AccelX'), title('GG Plot')
    xlabel('Lateral Gs'), ylabel('Longitudinal Gs')
    grid
    xlim([-2,2]), ylim([-2,2])
    hold on
    xline(0),yline(0)
end

if GGVPlot
    figure
    scatter3(T.AccelY, T.AccelX, T.Speed, '.'), title('GG Plot')
    xlabel('Lateral Gs'), ylabel('Longitudinal Gs'), zlabel('Speed (m/s)')
    xlim([-2,2]), ylim([-2,2])
    hold on
    xline(0),yline(0)
end

%% Aero, plots downforce vs speed, fits CLA
if aeroPlot
    include = ~or(isnan(T.Speed), isnan(T.SuspHeave));
    exclude =  T.Speed(include)<1; % m/s

    T.Speed = T.Speed * 22/9;
        
    % calculate downforce
    totalHeaveStiffness = 2*(350*4.45*39.37)*car.MR_F^2 + ... 
                          2*(250*4.45*39.37)*car.MR_R^2; % N/m
    
    % moving mean filter on data
    meanRangeSeconds = 1; % s
    meanTimestep = mean(diff(T.Time));
    order = meanRangeSeconds/meanTimestep;
    
    FzTotalSprings = T.SuspHeave(include)/1000*totalHeaveStiffness; % N
    FzTotalSprings = abs(movmean(FzTotalSprings, order));

    ft = fittype(@(a, b, x) -(a*car.aero.rho/2)*x.^2 + b); 
    f2 = fit(T.Speed(include), FzTotalSprings, ft, 'Exclude', exclude, 'start', [0 0]);
    totalHeaveStiffness = 2*car.k/car.MR_F^2 + 2*car.k/car.MR_R^2;
    CLA = round(f2.a, 2);

    figure
    if any(strcmp(T.Properties.VariableNames,'SuspForceFL'))
        subplot(2,1,1);
    end
    plot(f2,T.Speed(include), FzTotalSprings, exclude)
    title(['Speed Vs Downforce (suspension compression) | CLA = ' num2str(CLA) 'm^2'])
    xlabel('Speed (m/s)');
    ylabel('Force (N)');
    %ylim([0, Inf])
    grid    


    if any(strcmp(T.Properties.VariableNames,'SuspForceFL')) %% TODO: add motion ratio calculation
        FzTotalLoadCells = (T.SuspForceFL + T.SuspForceFR + T.SuspForceRL + T.SuspForceRR)/4 * 9.8; % N
        FzTotalLoadCells = abs(movmean(FzTotalLoadCells, order));

        ft = fittype(@(a, b, x) -(a*car.aero.rho/2)*x.^2 + b); 
        f2 = fit(T.Speed(include), FzTotalLoadCells, ft, 'Exclude', exclude, 'start', [0 0]);
        totalHeaveStiffness = 2*car.k/car.MR_F^2 + 2*car.k/car.MR_R^2;
        CLA = round(f2.a, 2);
    
        subplot(2,1,2);
        plot(f2,T.Speed(include), FzTotalLoadCells, exclude)
        title(['Speed Vs Downforce (load cells) | CLA = ' num2str(CLA) 'm^2'])
        xlabel('Speed (m/s)');
        ylabel('Force (N)');
        %ylim([0, Inf])
        grid
    end
    
    
end

% CLA = 2*F/(rho v^2)
% F = (CLA*rho/2)*v^2

%% Power Estimation
if powerPlot
    
    accelWheelCalc = [0; diff(T.Speed) ./ diff(T.Time)];
    select = (accelWheelCalc > 0.1) & (T.Time > 0) & (T.Time < Inf);
    power = ((accelWheelCalc .* car.M) + T.TheoreticalDrag)...
        .* T.Speed * 2.6; % W
%     power = (((abs(T.AccelX) * 9.8) .* car.M) + T.TheoreticalDrag)...
%         .* T.Speed; % W
    engineRPM = T.RPM;

    power = power(select);
    engineRPM = engineRPM(select);

    numZones = 30;
    percentilePower = 99.5;
    movmeanSmoothingOrder = 7;

    RPMbounds = linspace(0, max(engineRPM), numZones+1);

    fittedCurve = zeros(2,numZones);

    for i = 1:numZones
        selectRPM = (engineRPM > RPMbounds(i)) & (engineRPM < RPMbounds(i+1));
        fittedCurve(1, i) = mean(RPMbounds([i,i+1]));
        fittedCurve(2, i) = prctile(power(selectRPM), percentilePower);
    end

    fittedCurve(isnan(fittedCurve)) = 0;

    fittedCurve(2,:) = movmean(fittedCurve(2,:), movmeanSmoothingOrder);

    ft = fittype(@(a, b, c, d, x) a*x + b*x.^2 + c*x.^3 + d*x.^4); 
    powerFittedCurve = fit(fittedCurve(1,:)', fittedCurve(2,:)' * 0.00134102, ft, 'start', [0.0026 -1.79*10^(-6) 4.99*10^(-10) -3.047*10^(-14)]);

    figure

    select = power*0.00134102<50;
    
    hold on
    scatter(engineRPM(select), power(select) * 0.00134102, 'DisplayName', 'Instantaneous Power = (A_x \cdot M_{car} + F_{drag}) \cdot V_x');
    plot(fittedCurve(1,:), powerFittedCurve(fittedCurve(1,:)), 'LineWidth', 3, ...
         'DisplayName', [num2str(percentilePower) ' Percentile Power (WOT)']);

     
%     hold on
%     plot(fittedCurve(1,:), fittedCurve(2,:) * 0.00134102, 'LineWidth', 3, ...
%         'DisplayName', [num2str(percentilePower) ' Percentile Power (WOT)']);
     ylabel('Power (hp)');
     ylim([0 80])
% 
    % plot torque
    %torque = fittedCurve(2,:) ./ (fittedCurve(1,:)/60*2*pi); % Nm
    torque = powerFittedCurve(fittedCurve(1,:))' ./ (fittedCurve(1,:)/5252); % Ft-Lbs


    % * 0.73756
    yyaxis right
    plot(fittedCurve(1,:), torque, 'LineWidth', 3, ...
        'DisplayName', [num2str(percentilePower) ' Percentile Torque (WOT)']);
    legend();
    xlabel('Engine RPM');
    ylabel('Torque (Ft-Lbs)')
    ylim([0 80])
    xlim([2500 Inf])
    grid
    title('Engine Power Estimation');
end

%% FUNctions :D 
% this has to have been greg

function varargout = plotLine(T, timeRange, name, varargin)
    y = T.(name);
    y = y(T.Time>timeRange(1)&T.Time<timeRange(2));
    x = T.Time(T.Time>timeRange(1)&T.Time<timeRange(2));

    [varargout{1:nargout}] = plot(x,y,...
        'displayName', name, varargin{:});
    xlabel('time (s)')
    ylabel(getUnits(T, name));
end

function unit = getUnits(T, name)
    i = strcmp(T.Properties.VariableNames, name);
    unit = T.Properties.VariableUnits{i};
end

function T = setUnits(T, names, unit)
    if iscell(names)
        for j = 1:numel(names)
            i = strcmp(T.Properties.VariableNames, names{j});
            T.Properties.VariableUnits{i} = unit;
        end
    else
        i = strcmp(T.Properties.VariableNames, names);
        T.Properties.VariableUnits{i} = unit;
    end
end

