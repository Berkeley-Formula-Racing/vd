
steering_angle = [14.1 13.6 13.4 13.8 14.8 14.7 15.3 14.9 16.2];
lat_accel = [0.26 0.36 0.4 0.47 0.57 0.67 0.56 0.76 0.81];
groundspeed_left = [12.3 14.1 14.8 16.0 18.4 20.5 21.8 23.2 24.1];
groundspeed_right = [10.9 12.5 13.1 14.2 16.2 18.0 19.2 20.5 21.2];
velocity = (groundspeed_left+groundspeed_right)/2*0.44704; % m/s
pseudo_lat_accel = -velocity.^2/8.625/9.81;
scatter(pseudo_lat_accel,steering_angle)
