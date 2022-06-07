% Lapsim2
% Main steady-state lapsim script. This is traditional lapsim. It finds the
% g-g diagram over different velocities (max lateral and longitudinal
% acceleration at a given velocity), then uses those to predict performance
% in the dynamic events.
% HOW TO USE:
% 1) carConfig.m: define the car you want to test
% 2) run this script2
% 3) results are in the comp object, which is stored in the corresponding
% car object. Index into carCell to get the car you want, then open the
% comp object.
clear
setup_paths
carCell = carConfig(); %generate all cars to sim over
numCars = size(carCell,1);
time = struct();time.prev = 0; time.curr = 0;
tic
% Set numWorkers to number of cores for better performance
numWorkers = 0;
if numWorkers ~= 0
    disp('The parallel toolbox takes a few minutes to start.')
    disp('Set numWorkers to 0 for single-car runs')
end

for i = 1:numCars
    car = carCell{i,1};
    accelCar = carCell{i,2};
    fprintf("car %d of %d - starting g-g\n",[i numCars]);
    paramArr = gg2(car,numWorkers);
    fprintf("car %d of %d - g-g complete\n",[i numCars]);
    time.curr = floor(toc);
    fprintf("Stage Time: %d s; Total time elapsed: %d s\n",[time.curr-time.prev time.curr]);
    time.prev = time.curr;
    car = makeGG(paramArr,car); %post-processes gg data and stores in car
    comp = Events2(car,accelCar); 
    comp.calcTimes();       %run events and calc points
    car.comp = comp;        %store in array
    carCell{i,1} = car; %put updated car back into array. Matlab is pass by value, not pass by reference
    fprintf("car %d of %d - points calculated\n",[i numCars]);
    time.curr = floor(toc);
    fprintf("Stage Time: %d s; Total time elapsed: %d s\n",[time.curr-time.prev time.curr]);
end
fprintf("done\n");

%% Saving
save('FinalDriveSweep6-6-2022-num2.mat','carCell');

%% Points Plotting

disp("car 1 points: " + num2str(carCell{1,1}.comp.points.total));

% options
display_point_values_above_bar_flag = true;

label_cars_automatically_flag = true;

%automatic car labeling
automatic_label_name = 'Final Drive Ratio';
%automatic_label = @(car) (1/2+car.powertrain.G_d2_driving)/(1/2-car.powertrain.G_d2_driving);%TBR
automatic_label = @(car) car.powertrain.final_drive;%Car mass
%automatic_label = @(car) car.tire.gamma;%Car cda

% 1 to select, 0 to exclude
selected_categories = find([ ... 
     1 ... %Accel
     1 ... %Autocross
     1 ... %Endurance
     1 ... %Skidpad
     1 ... %Total
]);

plot_lapsim_points(carCell, display_point_values_above_bar_flag, true,...
    [], automatic_label_name, automatic_label, selected_categories);
%% Car Plotting

% select desired car object
desiredCarIndex = 1;
car = carCell{desiredCarIndex,1};

% set desired plots to 1
plot1 = 0; % velocity-dependent g-g diagram scatter plot
plot2 = 0; % velocity-dependent g-g diagram surface
plot3 = 0; % max accel for given velocity and lateral g w/ scattered interpolant
plot4 = 0; % max braking for given velocity and lateral w/ scattered interpolant
plot5 = 0; % 2D g-g diagram for velocity specified below (gg_vel)

g_g_vel = [12 14 26]; % can input vector to overlay different velocities

plot_choice = [plot1 plot2 plot3 plot4 plot5];
plotter(car,g_g_vel,plot_choice);

%% Event Plotting

% select desired comp object
comp = carCell{1,1}.comp;

% set desired plots to 1
plot1 = 0; % autocross track distance vs curvature
plot2 = 0; % endurance track distance vs curvature
plot3 = 0; % max possible velocity for given radius
plot4 = 0; % max possible long accel for given velocity
plot5 = 0; % accel event longitudinal velocity vs time
plot6 = 0; % accel event longitudinal accel vs time
plot7 = 0; % autocross gear shifts
plot8 = 0; % autocross slip angle vs distance

plot_choice = [plot1 plot2 plot3 plot4 plot5 plot6 plot7 plot8];
event_plotter(comp,plot_choice);

