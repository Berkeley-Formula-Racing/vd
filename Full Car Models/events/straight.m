function [time_vec,ending_vel,long_accel_vector,long_vel_vector] = straight(long_vel,length,...
    long_vel_interp,long_accel_interp,max_vel,accelCar)
% inputs: length of straight and starting velocity
% outputs: time and ending velocity 

distance = [0 linspace(0,length,10000)];
long_accel_vector = zeros(1,10000);
long_vel_vector = zeros(1,10000);

obj = accelCar.powertrain;
current_gear = 1;
shift_time_cumulative = 0;
start_shifting = false;

% Clutch launch. The capability lookup ties engine rpm rigidly to wheel speed,
% so at a standstill the engine is at ~0 rpm making ~0 torque and the car can
% barely start: the looked-up acceleration near v = 0 is essentially zero, and
% once rolling resistance is included it goes NEGATIVE, which makes the sqrt
% below imaginary and NaNs the whole run. A real car slips the clutch to hold
% the engine in its power band off the line. Model that minimally: below the
% speed where the drivetrain first overcomes resistance, launch at the
% capability available AT that speed rather than at the idle-rpm value. The
% dead zone is a fraction of a metre, so this is a faithful stand-in for the
% clutch, not a new degree of freedom -- and it is the ONLY place the launch
% is modelled, since Track_Solver's rollout comes through here too.
iPos = find(long_accel_interp > 0,1);
if isempty(iPos)
    vLaunch = 0; aLaunch = 0;      % nothing productive anywhere: leave as-is
else
    vLaunch = long_vel_interp(iPos);
    aLaunch = long_accel_interp(iPos);
end

time = zeros(1,10000);
for i = 1:10000 
    
    if long_vel>obj.switch_gear_velocities(current_gear)
        start_shifting = true;
    end
    
    if start_shifting ...
            && shift_time_cumulative < obj.shift_time
        long_accel = -accelCar.aero.drag(long_vel)/accelCar.M;
        shift_time_cumulative = shift_time_cumulative+time(i-1);
    else
        % basic kinematics equations
        long_accel = lininterp1(long_vel_interp,long_accel_interp,long_vel);
        % clutch-limited launch: hold the productive capability through the
        % dead zone the no-clutch lookup shows near a standstill
        if long_vel < vLaunch
            long_accel = aLaunch;
        end
    end

    if long_vel == max_vel
        long_accel = 0;
    end
     
    % reset shift time counter
    if shift_time_cumulative>obj.shift_time
        start_shifting = false;
        shift_time_cumulative = 0;
        current_gear = current_gear+1;
    end
    
    long_accel_vector(i) = long_accel;
    long_vel_initial = long_vel;
    long_vel_vector(i) = long_vel_initial;
    % max(...,0): a net-negative step (heavy braking, or resistance exceeding
    % drive) can drive the argument below zero and the sqrt imaginary. The
    % launch floor above keeps this from stalling the car at rest.
    long_vel = sqrt(max(long_vel^2+2*long_accel*(distance(i+1)-distance(i)),0));
    % limit top speed to max velocity. Latent on the current car -- the 75 m
    % accel run tops out at 28.31 m/s against a max_vel of 33.06 -- but without
    % it the long_vel == max_vel guard above can never fire, since only this
    % min() can produce that exact equality.
    long_vel = min(long_vel,max_vel);
    time(i) = 2*(distance(i+1)-distance(i))/(long_vel+long_vel_initial);
    
end

time_vec = cumsum(time);

% velocity at end of straight
ending_vel = long_vel;

end

