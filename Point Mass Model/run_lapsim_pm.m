clear; clc;

carCell = carConfig_pm(); 

for i = 1:numel(carCell)
    carCell{i} = events_pm(carCell{i});
end

