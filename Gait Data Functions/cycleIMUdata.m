function steps = cycleIMUdata(steps, accXYZ, times)
%% IMU Data Resampling
% This function takes in IMU data and calculated angles. Then using the
% steps indexes to resample the data. This is used to normalized the data
% based on heel strike to heel strike. StepsCell should contain the step
% indexes.  Gyro data and angles are optional arguments.
%
% The function assumes the following orientation of the data:
% +x: direction of walking
% +y: medial-to-lateral
% +z: lower-to-upper body
%
% Function by Erick Nunez

%% Variables to use
indexes = steps.indexes;
steps.Acc = cell(size(indexes));        steps.Times = {cell(size(indexes))};
steps.Vel = cell(size(indexes));        steps.Displace = cell(size(indexes));
steps.Distance = nan(size(indexes));   steps.Speed = nan(size(indexes));

%% Pull Accelerations
for i = 1:length(indexes)
    steps.Times{i} = times(indexes{i});
    steps.Acc{i} = accXYZ(indexes{i},:);
end

%% Find Step Velocities and displacements
for i = 1:length(indexes)
    steps.Vel{i} =       calculateInt(steps.Acc{i}, mean(diff(steps.Times{i})));
    steps.Displace{i} =  calculateInt(steps.Vel{i}, mean(diff(steps.Times{i})));
end

%% Finds HS to HS length
for i = 1:length(indexes)
    steps.Distance(i) = sqrt(steps.Displace{i}(end,1)^2 + steps.Displace{i}(end,3)^2);
end

%% Finds Speeds
for i = 1:length(indexes)
   steps.Speed(i) = steps.Distance(i)/(steps.Times{i}(end)-steps.Times{i}(1)); 
end


end

%% Functions
function dataOut = calculateInt(rate, dt)
dataOut = zeros(size(rate));
for i = 2:length(dataOut)
    dataOut(i,:) = dataOut(i-1,:) + (rate(i,:) + rate(i-1,:)).*dt/2;
end
end