clear; clc;

carCell = carConfig_pm(); 


fprintf('\n--- factorial results ---\n');
for i = 1:numel(carCell)
    fprintf('case %d: lap = %.2f s\n', i, carCell{i}.lap_time);
end


figure; hold on; grid on;
for i = 1:min(5, numel(carCell))
    s = carCell{i}.spec;    
    if i == 1
        trk = track_pm.loadFromMichiganMat('michigantrack2024.mat');
    end
    plot(trk.s, carCell{i}.v_profile, 'DisplayName', ...
        sprintf('case %d (%.2fs)', i, carCell{i}.lap_time));
end
xlabel('s [m]'); ylabel('v [m/s]');
title('Speed profiles (first few cases)');
legend show;
