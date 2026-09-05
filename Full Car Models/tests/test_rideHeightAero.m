function tests = test_rideHeightAero
tests = functiontests(localfunctions);
end

function testConfiguredCarReportsCoupledRideHeightAndAeroLoads(testCase)
[cars,~] = carConfig();
car = cars{1,1};

state = [0,0,20,0,0,0,0,0,0];
m = car.metrics(state);

verifyTrue(testCase,car.rideHeightAero.enabled);
verifyGreaterThan(testCase,m.aero_iterations,0);
verifyLessThan(testCase,m.front_ride_height_in, ...
    car.rideHeightAero.static_front_ride_height_in);
verifyLessThan(testCase,m.rear_ride_height_in, ...
    car.rideHeightAero.static_rear_ride_height_in);
verifyEqual(testCase,m.aero_downforce_front_N + m.aero_downforce_rear_N, ...
    m.downforce,'AbsTol',1e-9);
verifyEqual(testCase,m.front_ride_height_offset_in, ...
    m.front_ride_height_in-car.rideHeightAero.map_reference_front_ride_height_in, ...
    'AbsTol',1e-12);
verifyEqual(testCase,mean([m.ride_height_FL_in m.ride_height_FR_in]), ...
    m.front_ride_height_in,'AbsTol',1e-12);
verifyEqual(testCase,mean([m.ride_height_RL_in m.ride_height_RR_in]), ...
    m.rear_ride_height_in,'AbsTol',1e-12);
end

function testAccelerationCarRetainsStaticAeroFallback(testCase)
[cars,~] = carConfig();
accelCar = cars{1,2};

m = accelCar.metrics([0,0,20,0,0,0,0,0,0]);

verifyFalse(testCase,accelCar.aero.hasMap());
verifyEqual(testCase,m.downforce,accelCar.aero.lift(20),'AbsTol',1e-12);
verifyEqual(testCase,m.drag,accelCar.aero.drag(20),'AbsTol',1e-12);
verifyTrue(testCase,isnan(m.front_ride_height_in));
verifyTrue(testCase,isnan(m.rear_ride_height_in));
end
