function carCell = toLegacyCells(cases)
%TOLEGACYCELLS Pack named cases as lap/acceleration Cars without rebuilding.
carCell = cell(numel(cases),2);
for i = 1:numel(cases)
    carCell{i,1} = cases(i).lapCar;
    carCell{i,2} = cases(i).accelCar;
end
end
