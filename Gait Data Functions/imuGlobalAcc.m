function globalAccXYZ = imuGlobalAcc(anglesXYZ, accXYZ)
%% Global accelerations
% Function finds global accelerations give angles and local accelerations
% of IMUs. Function assumed the following orientation of the data:
% +x: direction of walking
% +y: medial to Left Side
% +z: normal to ground (vertical)
% 
% Global Accelerations are founds as such:
% Ax = Direction of walking.
% Ay = Lateral direction.
% Az = Vertical direction.
% A = Ry*a - G
%
% function by erick nunez

%% variables
g = 9.81;
ax = accXYZ(:,1);       ay = accXYZ(:,2);       az = accXYZ(:,3);
alpha = (pi/180)*anglesXYZ(:,1); 
beta  = (pi/180)*anglesXYZ(:,2);  
theta = (pi/180)*anglesXYZ(:,3);

%% global accelerations calculated
Ax =  ax.*cos(beta) + az.*sin(beta);
Ay =  ay;
Az = -ax.*sin(beta) + az.*cos(beta) - g;

%% Output
globalAccXYZ = [Ax Ay Az];
end