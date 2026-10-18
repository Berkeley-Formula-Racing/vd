%% car_setup_example
% Edit a vehicle configuration and build named cases without running a solver.
entrypointDir = fileparts(mfilename('fullpath'));
run(fullfile(entrypointDir,'bootstrap.m'));
[setup,eventParams] = vehicleConfig.defaultSetup();
variant = setup;
variant.parameters.wheelbase = 1.60;
variant.parameters.final_drive = 3.2;
car = vehicleConfig.buildCar(variant,"lap");
Q = table([155;170],'VariableNames',{'mass'});
[cases,eventParams,designTable] = ...
    carConfig("Explicit",Q,"legacy","nominal","cases");
lapCar = cases(2).lapCar;
