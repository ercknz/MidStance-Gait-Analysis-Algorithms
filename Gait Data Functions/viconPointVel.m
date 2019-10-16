function [velocity, speed] = viconPointVel(pointXYZ, times)
%% viconPointVel Function
% This function takes in a the XYZ cordinates of a marker and generates the
% XYZ velocities and the speed. 
% 
% Function by erick nunez

%% Calculates velocity and speed
u = gradient(pointXYZ(:,1),times);
v = gradient(pointXYZ(:,2),times);
w = gradient(pointXYZ(:,3),times);
velocity = [u,v,w];
speed = sqrt(u.^2 + v.^2 + w.^2);

end