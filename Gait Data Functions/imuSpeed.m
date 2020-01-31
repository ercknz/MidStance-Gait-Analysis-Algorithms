function out = imuSpeed(imu, otherIndexes)
%% imuSpeed
% Function assumed the following orientation of the data:
% +x: direction of walking
% +y: medial to Left Side
% +z: normal to ground (vertical)
%
% Function by erick nunez

%% Variables to use
if nargin == 1
    indexes = imu.steps.indexes;
elseif nargin == 2
    indexes = otherIndexes;
end
acc = cell(size(indexes));        times = {cell(size(indexes))};
vel = cell(size(indexes));        speed = cell(size(indexes));
meanSpeed = [];


%% Pull Accelerations
for i = 1:length(indexes)
    times{i} = imu.times(indexes{i});
    acc{i} = imu.globalAcc(indexes{i},:);
end

%% Find Step Velocities and speeds
for i = 1:length(indexes)
    vel{i} = calculateInt(acc{i}, mean(diff(times{i})));
    speedProfile = sqrt(vel{i}(:,1).^2 + vel{i}(:,2).^2 + vel{i}(:,3).^2);
    speed{i} = speedProfile;
    meanSpeeds(i) = mean(speedProfile);
end

%% Save outputs
out.speed = speed;
out.meanSpeeds = meanSpeeds;
out.vel = vel;

end

%% Additional Functions
function dataOut = calculateInt(rate, dt)
dataOut = zeros(size(rate));
for i = 2:length(dataOut)
    dataOut(i,:) = dataOut(i-1,:) + (rate(i,:) + rate(i-1,:)).*dt/2;
end
end