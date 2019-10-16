function compAngles = compFusionAngles(highPass, accXYZ, gyroXYZ, magXYZ, sampleFreq)
%% imuFusionAngles
% Finding angles and fuses data from IMUs. Function takes in [samples,3]
% matrices which contain XYZ data from an IMU. Then finds the angles using
% the below formulas. highPass is the weight for the complimentary
% filter. This assumes the following orientation:
% +x: direction of walking
% +y: medial-to-lateral
% +z: lower-to-upper body
% 
% accAng(x) = atan(acc(y)/sqrt(acc(x)^2+acc(z)^2))
% accAng(y) = atan(acc(x)/sqrt(acc(y)^2+acc(z)^2))
% 
% magAng(x) = mag(x)*cos(accAng(y)) ...
%            + mag(y)*sin(accAng(x))*sin(accAng(y)) ...
%            + mag(z)*cos(accAng(x))*sin(accAng(y))
% magAng(y) = mag(y)*cos(accAng(x)) ...
%            - mag(z)*sin(accAng(x))
% accAng(z) = atan(-magAng(y)/magAng(x)) 
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
accAng(:,1) = atan2d(accXYZ(:,2),sqrt((accXYZ(:,1).^2)+(accXYZ(:,3).^2)));
accAng(:,2) = atan2d(accXYZ(:,1),sqrt((accXYZ(:,2).^2)+(accXYZ(:,3).^2)));
Mx = magXYZ(:,1) .* cosd(accAng(:,2)) ...
    + magXYZ(:,2) .* sind(accAng(:,1)) .* sind(accAng(:,2)) ...
    - magXYZ(:,3) .* cosd(accAng(:,1)) .* sind(accAng(:,2));
My = magXYZ(:,2) .* cosd(accAng(:,1)) ...
    - magXYZ(:,3) .* sind(accAng(:,1));
accAng(:,3) = atan2d(-My,Mx);

%% Calculates angles using complimentary filter
compAngles(1,:) = highPass .* gyroXYZ(1,:) .* dt + lowPass .* accAng(1,:);
for i = 2:length(gyroXYZ)
    compAngles(i,:) = highPass .* (compAngles(i-1,:) + gyroXYZ(i,:) .* dt) + lowPass .* accAng(i,:);
end

end