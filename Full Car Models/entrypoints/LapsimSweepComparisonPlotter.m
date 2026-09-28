%% setup
entrypointDir = fileparts(mfilename('fullpath'));
run(fullfile(entrypointDir, 'bootstrap.m'));

%% settings
%file_names = {'Lapsim18R25B_non_turbo.mat','Lapsim18R25B_turbo.mat'};
%file_names = {'Lapsim_B22_no_undertray.mat','Lapsim_B22_with_undertray.mat'};
file_names = {'FinalDriveSweep6-6-2022.mat'};

% options
display_point_values_above_bar_flag = true;

label_cars_automatically_flag = true;

%manual car labeling
manual_car_labels = {'NA', 'Turbo'};
%manual_car_labels = {'No Undertray', 'Undertray'};

%automatic car labeling
automatic_label_name = 'Torque Bias Ratio';
%automatic_label = @(car) (1/2+car.powertrain.G_d2_driving)/(1/2-car.powertrain.G_d2_driving);%TBR
%automatic_label = @(car) car.M;%Car mass
automatic_label = @(car) car.powertrain.final_drive;%Car mass

% 1 to select, 0 to exclude
selected_categories = find([ ...
     1 ... %Accel
     1 ... %Autocross
     1 ... %Endurance
     1 ... %Skidpad
     1 ... %Total
]);

%combine cars from each file into one large car cell array
carCellComplete = {};
num_files = numel(file_names);
for i = 1:num_files
    carCellComplete = [carCellComplete; load(file_names{i}).carCell];
end

plot_lapsim_points(carCellComplete, display_point_values_above_bar_flag, label_cars_automatically_flag,...
    manual_car_labels, automatic_label_name, automatic_label, selected_categories);
