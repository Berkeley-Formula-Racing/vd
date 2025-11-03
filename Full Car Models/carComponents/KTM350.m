function torque_fn = KTM350()

% Preliminary Spark Tune (from Feb 2019)
points = [...
%0	0
6000	22.5535
6500	24.2606
7000	25.6219
7500	26.7192
8000	27.1255
8500	26.9020
9000	25.8048
9500	25.6219
10000	24.7076
10500	22.249
11000	20.8674
11500	20.461
12000	19.5263
12500	18.0837];

RPM = transpose(points(:,1));
TQ = transpose(points(:,2));

torque_fn = [RPM; TQ];

plot(RPM, TQ);