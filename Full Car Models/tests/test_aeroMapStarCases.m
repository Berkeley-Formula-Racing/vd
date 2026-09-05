function tests = test_aeroMapStarCases
tests = functiontests(localfunctions);
end

function testBuildsFiveIndependentStarsAroundMapReference(testCase)
modelRoot = fileparts(which('carConfig'));
cd(modelRoot); setup_paths
[baseCell,~] = carConfig();

[carCell,plan] = aeroMapStarCases(baseCell);

verifySize(testCase,carCell,[21 2]);
verifyEqual(testCase,height(plan),21);
verifyEqual(testCase,plan.swept_parameter(1),"baseline");
verifyEqual(testCase,plan.front_static_ride_height_in(1),0,'AbsTol',1e-12);
verifyEqual(testCase,plan.rear_static_ride_height_in(1),0,'AbsTol',1e-12);
verifyEqual(testCase,plan.cla_scale(1),1,'AbsTol',1e-12);
verifyEqual(testCase,plan.cda_scale(1),1,'AbsTol',1e-12);
verifyEqual(testCase,plan.cop_offset(1),0,'AbsTol',1e-12);

parameters = ["FrontRideHeightIn" "RearRideHeightIn" "ClAScale" "CdAScale" "CoPOffset"];
for p = parameters
    verifyEqual(testCase,nnz(plan.swept_parameter == p),4);
end

i = find(plan.swept_parameter == "FrontRideHeightIn" & ...
    abs(plan.swept_value-0.25) < 1e-12,1);
car = carCell{i,1};
verifyEqual(testCase,car.rideHeightAero.static_front_ride_height_in,0.25,'AbsTol',1e-12);
verifyEqual(testCase,car.rideHeightAero.static_rear_ride_height_in,0,'AbsTol',1e-12);
verifyEqual(testCase,car.aero.map.claScale,1,'AbsTol',1e-12);
verifyEqual(testCase,car.aero.map.cdaScale,1,'AbsTol',1e-12);
verifyEqual(testCase,car.aero.map.copOffset,0,'AbsTol',1e-12);
end

function testSensitivityPlanLabelsMapInputsAndUnits(testCase)
modelRoot = fileparts(which('carConfig'));
cd(modelRoot); setup_paths
[baseCell,~] = carConfig();
[~,plan] = aeroMapStarCases(baseCell);

info = aeroMapSensitivityPlan(plan);

verifyEqual(testCase,info.baseIdx,1);
verifyEqual(testCase,string({info.parameters.name}), ...
    ["FrontRideHeightIn" "RearRideHeightIn" "ClAScale" "CdAScale" "CoPOffset"]);
verifyEqual(testCase,string({info.parameters.unit}), ...
    ["in" "in" "scale" "scale" "front-balance fraction"]);
end
