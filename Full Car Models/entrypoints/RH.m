options = struct();
options.numWorkers = 15;
options.saveRawResults = true;
options.outputPath = "C:\Sim Results\ride_height_stiffness_225.mat";

result = run_ride_height_stiffness_sweep(options);