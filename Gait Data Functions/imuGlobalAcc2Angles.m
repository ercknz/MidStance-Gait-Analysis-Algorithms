function globalAccXYZ = imuGlobalAcc2Angles(anglesXYZ, accXYZ)
%% Global accelerations
% Function finds global accelerations give angles and local accelerations
% of IMUs. Function assumed the following orientation of the data:
% +x: direction of walking
% +y: medial to Left Side
% +z: normal to ground (vertical)
% 
% Global Accelerations are founds as such:
% Ax = direction of walking.
% Ay = lateral direction. 
% Az = Vertical direction.
% A = Ry*Rx*a - g
%
% function by erick nunez

%% variables
g = 9.81;
ax = accXYZ(:,1);       ay = accXYZ(:,2);       az = accXYZ(:,3);
alpha = (pi/180)*anglesXYZ(:,1); 
beta  = (pi/180)*anglesXYZ(:,2);  
theta = (pi/180)*anglesXYZ(:,3);

%% global accelerations calculated
Ax =  ax.*cos(beta) + ay.*sin(alpha).*sin(beta) + az.*cos(alpha).*sin(beta);
Ay =                  ay.*cos(alpha)            - az.*sin(alpha);
Az = -ax.*sin(beta) + ay.*cos(beta).*sin(alpha) + az.*cos(alpha).*cos(alpha) - g;

%% Output
globalAccXYZ = [Ax Ay Az];

end