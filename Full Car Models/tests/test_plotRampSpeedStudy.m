function tests = test_plotRampSpeedStudy
tests = functiontests(localfunctions);
end

function testOverlaysBalanceMetricsForEveryCar(testCase)
% Would fail if a car is omitted, F/R aero loads are not separated, or the
% total handling metric is not shown as a zero-centred K_linear trace.
R(1) = mockRamp([5;10;15], [100;400;900], [0.53;0.54;0.55], ...
    [0.48;0.48;0.48], [0.8;0.2;-0.4]);
R(2) = mockRamp([5;10;15], [90;350;760], [0.49;0.50;0.51], ...
    [0.52;0.52;0.52], [1.1;0.5;-0.1]);

shown = ["aero_front_load" "aero_rear_load" "aero_balance" "handling_balance"];
fig = plotRampSpeedStudy(R,["baseline" "comparison"], ...
    struct('visible','off','outputs',shown));
cleaner = onCleanup(@() close(fig)); %#ok<NASGU>

ax = flipud(findall(fig,'Type','axes'));
verifyEqual(testCase,numel(ax),4);
verifyEqual(testCase,string(ax(1).YLabel.String),"front aero load (N)");
verifyEqual(testCase,string(ax(2).YLabel.String),"rear aero load (N)");
verifyEqual(testCase,string(ax(3).YLabel.String),"front aero balance (%)");
verifyEqual(testCase,string(ax(4).YLabel.String),"K_{linear} (deg/g)");
verifyGreaterThanOrEqual(testCase,numel(findall(ax(1),'Type','line')),4);
verifyGreaterThanOrEqual(testCase,numel(findall(ax(4),'Type','line')),2);
verifyEqual(testCase,numel(findall(ax(4),'Type','constantline')),1);
end

function testPlotsSelectedFrontRearRideHeightAndCamberOutputs(testCase)
R = mockRamp([5;10],[100;400],[0.53;0.54],[0.48;0.48],[0.8;0.2]);
shown = ["aero_front_load" "aero_rear_load" "front_ride_height" "rear_camber" "pitch_angle"];

fig = plotRampSpeedStudy(R,'baseline',struct('visible','off','outputs',shown));
cleaner = onCleanup(@() close(fig)); %#ok<NASGU>

ax = flipud(findall(fig,'Type','axes'));
verifyEqual(testCase,numel(ax),5);
verifyEqual(testCase,string(ax(1).YLabel.String),"front aero load (N)");
verifyEqual(testCase,string(ax(2).YLabel.String),"rear aero load (N)");
verifyEqual(testCase,string(ax(3).YLabel.String),"front ride height (in)");
verifyEqual(testCase,string(ax(4).YLabel.String),"rear camber magnitude (deg)");
verifyEqual(testCase,string(ax(5).YLabel.String),"pitch angle (deg)");
end

function testPlotsSelectedFrontAndRearShockTravelOutputs(testCase)
R = mockRamp([5;10],[100;400],[0.53;0.54],[0.48;0.48],[0.8;0.2]);

fig = plotRampSpeedStudy(R,'baseline',struct('visible','off', ...
    'outputs',["front_shock_travel" "rear_shock_travel"]));
cleaner = onCleanup(@() close(fig)); %#ok<NASGU>

ax = flipud(findall(fig,'Type','axes'));
verifyEqual(testCase,numel(ax),2);
verifyEqual(testCase,string(ax(1).YLabel.String),"front shock compression (in)");
verifyEqual(testCase,string(ax(2).YLabel.String),"rear shock compression (in)");
end

function testAcceptsOneCharacterVectorLabelForOneCar(testCase)
R = mockRamp([5;10],[100;400],[0.53;0.54],[0.48;0.48],[0.8;0.2]);

fig = plotRampSpeedStudy(R,'baseline',struct('visible','off'));
cleaner = onCleanup(@() close(fig)); %#ok<NASGU>

verifyTrue(testCase,isgraphics(fig));
end

function R = mockRamp(v,df,cop,lltd,k)
R = struct();
R.perSpeed = table(v,df,cop,lltd,k,df.*cop,df.*(1-cop), ...
    -0.01*v,-0.008*v,1+0.01*v,0.8+0.01*v,-0.1*v,0.01*v,0.008*v, ...
    'VariableNames',{'vCar','downforce','aero_balance','mech_balance','K_linear', ...
    'aero_downforce_front_N','aero_downforce_rear_N', ...
    'front_ride_height_in','rear_ride_height_in', ...
    'front_camber_deg','rear_camber_deg','pitch_angle_deg', ...
    'front_shock_travel_in','rear_shock_travel_in'});
R.settings = struct('mode','coast','nRamp',12);
end
