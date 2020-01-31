function globalAccXYZ = imuGlobalAcc3Angles(anglesXYZ, accXYZ)
%% Global accelerations
% Function finds global accelerations give angles and local accelerations
% of IMUs. Function assumed the following orientation of the data:
% +x: direction of walking
% +y: medial to Left Side
% +z: normal to ground (vertical)
% 
% Global Accelerations are founds as such:
% Ax = ax*Cos(thetaY) + az*Sin(thetaY) -> direction of walking.
% Ay = ay
% Az = -ax*Sin(thetaY) + az*Cos(thetaY) - gravity -> Vertical direction.
%
% function by erick nunez

%% variables
g = 9.81;
ax = accXYZ(:,1);       ay = accXYZ(:,2);       az = accXYZ(:,3);
alpha = (pi/180)*anglesXYZ(:,1); 
beta  = (pi/180)*anglesXYZ(:,2);  
theta = (pi/180)*anglesXYZ(:,3);

%% global accelerations calculated
Ax =  ax.*(cos(beta).*cos(theta) + sin(alpha).*sin(beta).*sin(theta)) + ay.*(sin(alpha).*sin(beta).*cos(theta) - cos(beta).*sin(theta)) + az.*cos(alpha).*sin(beta);
Ay =  ax.*cos(alpha).*sin(theta)                                      + ay.*cos(alpha).*cos(theta)                                      - az.*sin(alpha);
Az =  ax.*(cos(beta).*sin(alpha).*sin(theta) - sin(beta).*cos(theta)) + ay.*(sin(beta).*sin(theta) + cos(beta).*sin(alpha).*cos(theta)) + az.*cos(alpha).*cos(beta) - g;

%% Output
globalAccXYZ = [Ax Ay Az];

end