function compAngles = compFusionAngles(highPass, accXYZ, gyroXYZ, magXYZ, sampleFreq)
%% imuFusionAngles
% Finding angles and fuses data from IMUs. Function takes in [samples,3]
% matrices which contain XYZ data from an IMU. Then finds the angles using
% the below formulas. highPass is the weight for the complimentary
% filter. This assumes the following orientation:
% +x: direction of walking
% +y: medial to Left Side
% +z: normal to ground (vertical)
% 
% roll = atan(acc(y)/acc(z))
% pitch = atan(-acc(x)/sqrt(acc(y)^2+acc(z)^2))
% 
% M(x) = mag(x)*cos(pitch) + mag(z)*sin(pitch)
% M(y) = mag(x)*sin(roll)*sin(pitch) + mag(y)*cos(roll) - mag(z)*sin(roll)*cos(pitch)
% yaw = atan(M(y)/M(x)) 
% accAng = [roll, pitch, yaw]
%                              
% 
% compAng(i) = highPass*(compAng(i-1)+gyro*dt)+lowPass*accAng
% 
% Function by erick nunez

%% Sets up variables and checks them
dt = 1/sampleFreq;
if highPass <= 1 && highPass >= 0
    lowPass = 1 - highPass;
else
    error('highPass value invalid')
end

%% Calculates Accel angles based on Accel/mag readings
roll = atan2d(accXYZ(:,2),accXYZ(:,3));
pitch = atan2d(-accXYZ(:,1),sqrt((accXYZ(:,2).^2)+(accXYZ(:,3).^2)));

Mx = magXYZ(:,1).*cosd(pitch) + magXYZ(:,3).*sind(pitch);

My = magXYZ(:,1).*sind(roll).*sind(pitch) ...
    + magXYZ(:,2).*cosd(roll) ...
    - magXYZ(:,3).*sind(roll).*cosd(pitch);

yaw = atan2d(My,Mx);

accAng = [roll, pitch, yaw];

%% Calculates angles using complimentary filter
compAngles(1,:) = highPass .* gyroXYZ(1,:) .* dt + lowPass .* accAng(1,:);
for i = 2:length(gyroXYZ)
    compAngles(i,:) = highPass .* (compAngles(i-1,:) + gyroXYZ(i,:) .* dt) + lowPass .* accAng(i,:);
end

end