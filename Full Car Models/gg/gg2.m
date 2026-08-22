function paramArr = gg2(car,numWorkers)
% creates velocity-dependent g-g diagram 
% describes max lateral acceleration, max longitudinal acceleration for
%   certain velocity
% inputs car: car object, numWorkers: number of workers for parallelization
% outputs paramArr: array of ParamSets (contains info about each
%   optimization point). 
%   rows: different longitudinal velocities, columns: different lat accels

minV = 5; maxV = car.max_vel;
longVinterval = 1; % m/s
latAgrid = 20;

% maxV is car.max_vel, not max_vel-0.5. That 0.5 standoff dates to the 2019
% first revision with no explanation, and the obvious reason for it -- that the
% lateral solve goes marginal where the car needs all its grip just to hold top
% speed -- is not what happens. The full row converges at exactly max_vel on
% every car on file and at four different final drives spanning max_vel 22.3 to
% 39.1, worst lateral residual 9.3e-9. The guard that was actually needed is
% the -0.001 already inside car.max_vel, which keeps the engine off redline.
%
% This DOES move lap times on a car that can reach the band. Today's car cannot
% -- its best longitudinal acceleration crosses zero at 30.3 m/s, well under
% the old top row -- so autocross and endurance are unchanged to 1e-5 s. On a
% 40/9 final drive with half the drag it is worth +3.05 s of endurance, because
% the old grid froze a POSITIVE acceleration at the top row and pinned peak
% speed at max_vel instead of letting it settle at its true terminal speed.
% That is the car this exists for.
%
% Snap the last colon node onto maxV rather than appending it, when the two are
% within a quarter interval. Appending leaves a final gap equal to the
% fractional part of maxV -- 0.0596 m/s on the stock final drive, and 1e-4 at
% some others -- and two velocity rows that close are coincident to the
% scatteredInterpolant the lap sim triangulates, which unique() will not
% collapse because they are not equal. Snapping bounds the last interval to
% [0.25, 1.25] m/s and costs nothing.
longVelArr = minV:longVinterval:maxV;
if isempty(longVelArr)
    longVelArr = maxV;                       % max_vel below minV: one row
elseif maxV - longVelArr(end) > 0.25*longVinterval
    longVelArr(end+1) = maxV;
else
    longVelArr(end) = maxV;
end
paramArr(numel(longVelArr),latAgrid) = ParamSet();

% iterate through velocities
parfor (c1 = 1:numel(longVelArr),numWorkers) %parfor
    longVel = longVelArr(c1);
    [maxLatx,maxLatLatAccel,maxLatLongAccel,maxLatx0] = max_lat_accel(longVel,car);
    latAccelArr = linspace(0.1,maxLatLatAccel-0.1,latAgrid);
    row = ParamSet();
    row(numel(latAccelArr)) = ParamSet();
    % iterate through lateral accelerations
    for c2 = 1:numel(latAccelArr)
        latAccel = latAccelArr(c2);
        %disp([longVel, latAccel]);
        %disp('here1')
        [xAccel,longAccel,longAccelx0] = max_long_accel_cornering(longVel,latAccel,car);
        %disp('here2')
        [xBraking,longDecel,brakingDecelx0] = max_braking_decel_cornering(longVel,latAccel,car);
        %disp('here3');
        carParams = ParamSet(car,longVel); 
        carParams = carParams.setMaxLatParams(maxLatx,maxLatLatAccel,maxLatLongAccel,maxLatx0);
        carParams = carParams.setMaxAccelParams(xAccel,longAccel,latAccel,longAccelx0);
        carParams = carParams.setMaxDecelParams(xBraking,longDecel,latAccel,brakingDecelx0);
        row(c2) = carParams;
    end
    paramArr(c1,:) = row;
end