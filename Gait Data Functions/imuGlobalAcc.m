function globalAccXYZ = imuGlobalAcc(anglesXYZ, accXYZ)
%% Global accelerations
% Function finds global accelerations give angles and local accelerations
% of IMUs. Function assumed the following orientation of the data:
% +x: direction of walking
% +y: medial-to-lateral
% +z: lower-to-upper body
% 
% Global Accelerations are founds as such:
% Ax = axCos(theta)-azSin(theta) -> direction of walking.
% Ay is ignored at the moment
% Az = axSin(theta)+azCos(theta)-gravity -> Vertical direction.
%
% function by erick nunez

%% variables
gravity = 9.81;

%% global accelerations calculated
globalAccXYZ(:,1) = accXYZ(:,1) .* cosd(anglesXYZ(:,2)) - accXYZ(:,3) .* sind(anglesXYZ(:,2));
globalAccXYZ(:,2) = nan(size(accXYZ(:,2)));
globalAccXYZ(:,3) = accXYZ(:,1) .* sind(anglesXYZ(:,2)) + accXYZ(:,3) .* cosd(anglesXYZ(:,2)) - gravity;

end