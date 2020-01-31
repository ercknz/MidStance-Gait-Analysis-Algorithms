function [waist, rightAnkle, leftAnkle] = gaitSpeed(measuredWaist, waist, rightAnkle, leftAnkle)
%% GaitVelocity
% This functions takes in the measured waist data, the separated right
% steps, and the separated  left steps. Then, based on the leading foot, it
% pulls the waist accelerations per cycle (HS-HS). Next it calculates the
% velocities for the waist and ankles. Data is plotted and saved locally.
% Function assumed the following orientation of the data:
% +x: direction of walking
% +y: medial to Left Side
% +z: normal to ground (vertical)
%
% Function by erick nunez

%% Uses leading Foot to find waist data.
if rightAnkle.steps.indexes{1}(1) < leftAnkle.steps.indexes{1}(1)
    waist.Step = imuDataResampling(rightAnkle.steps.indexes, ...
        measuredWaist.times, waist.globalAcc);
else
    waist.Step = imuDataResampling(leftAnkle.steps.indexes, ...
        measuredWaist.times, waist.globalAcc);
end

%% Finds Velocities
waist.Step.vel      =  calculateInt(     waist.Step.meanAcc,      mean(diff(waist.Step.avgTimes)));
rightAnkle.Step.vel =  calculateInt(rightAnkle.Step.meanAcc, mean(diff(rightAnkle.Step.avgTimes)));
leftAnkle.Step.vel  =  calculateInt( leftAnkle.Step.meanAcc,  mean(diff(leftAnkle.Step.avgTimes)));

%% Finds displacement
waist.Step.disp      =  calculateInt(     waist.Step.vel,      mean(diff(waist.Step.avgTimes)));
rightAnkle.Step.disp =  calculateInt(rightAnkle.Step.vel, mean(diff(rightAnkle.Step.avgTimes)));
leftAnkle.Step.disp  =  calculateInt( leftAnkle.Step.vel,  mean(diff(leftAnkle.Step.avgTimes)));

%% total Distance
waist.Step.distance      =   sqrt(     waist.Step.disp(end,1)^2 +      waist.Step.disp(end,3)^2);
rightAnkle.Step.distance =   sqrt(rightAnkle.Step.disp(end,1)^2 + rightAnkle.Step.disp(end,3)^2);
leftAnkle.Step.distance  =   sqrt( leftAnkle.Step.disp(end,1)^2 +  leftAnkle.Step.disp(end,3)^2);
% waist.Step.distance      =   sqrt(     waist.Step.disp(end,1)^2 +      waist.Step.disp(end,2)^2 +      waist.Step.disp(end,3)^2);
% rightAnkle.Step.distance =   sqrt(rightAnkle.Step.disp(end,1)^2 + rightAnkle.Step.disp(end,2)^2 + rightAnkle.Step.disp(end,3)^2);
% leftAnkle.Step.distance  =   sqrt( leftAnkle.Step.disp(end,1)^2 +  leftAnkle.Step.disp(end,2)^2 +  leftAnkle.Step.disp(end,3)^2);

%% Walking speed
waist.avgSpeed      =      waist.Step.distance/      waist.Step.avgTimes(end);
rightAnkle.avgSpeed = rightAnkle.Step.distance/ rightAnkle.Step.avgTimes(end);
leftAnkle.avgSpeed  =  leftAnkle.Step.distance/  leftAnkle.Step.avgTimes(end);

end

%% Additional Functions
function dataOut = calculateInt(rate, dt)
dataOut = zeros(size(rate));
for i = 2:length(dataOut)
    dataOut(i,:) = dataOut(i-1,:) + (rate(i,:) + rate(i-1,:)).*(dt/2);
end
end

% function dataOut = resetDrift(dataIn, times)
% tStart = times(1); tEnd = times(end);
% dataOut = nan(size(dataIn));
% for i = 1:length(dataIn)
%    dataOut(i,:) = ((tEnd - times(i))/(tEnd - tStart)) .* dataIn(i,:);
% end
% end