function [cuoSteerDeg,steerFromYawDeg] = cuoSteerFromYaw(steerAvgDeg,yawRate,vCar,wheelbase)
%CUOSTEERFROMYAW Road-wheel steer residual relative to yaw kinematics.
% Positive CUOsteerFromYaw is understeer: the car needs more road-wheel
% steer than a neutral bicycle model at the same yaw rate and speed.

steerFromYawDeg = atan2(wheelbase.*yawRate,vCar)*180/pi;
cuoSteerDeg = steerAvgDeg - steerFromYawDeg;
end
