function steps = cycleIMUdata(steps, accXYZ, times)
%% cycleIMUdata
% This function takes in IMU data and the seperated steps, either HS to HS
% or MS to MS, and find the velocity, distance, displacement, and speed per
% step.
%
% The function assumes the following orientation of the data:
% +x: direction of walking
% +y: medial to Left Side
% +z: normal to ground (vertical)
%
% Function by Erick Nunez

%% Variables to use
indexes = steps.indexes;
steps.Acc = cell(size(indexes));        steps.Times = {cell(size(indexes))};
steps.Vel = cell(size(indexes));        steps.Displace = cell(size(indexes));
steps.Distance = nan(size(indexes));    steps.Speed = nan(size(indexes));
steps.Duration = nan(size(indexes));

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

%% Finds segment length
for i = 1:length(indexes)
    %     steps.Distance(i) = sqrt(steps.Displace{i}(end,1)^2 + steps.Displace{i}(end,3)^2);
    steps.Distance(i) = sqrt(steps.Displace{i}(end,1)^2 + steps.Displace{i}(end,2)^2 + steps.Displace{i}(end,3)^2);
end

%% Finds Speeds
for i = 1:length(indexes)
    steps.Duration(i) = steps.Times{i}(end)-steps.Times{i}(1);
    steps.Speed(i) = steps.Distance(i)/steps.Duration(i);
end

end

%% Additional Functions
function dataOut = calculateInt(rate, dt)
dataOut = zeros(size(rate));
for i = 2:length(dataOut)
    dataOut(i,:) = dataOut(i-1,:) + (rate(i,:) + rate(i-1,:)).*dt/2;
end
end
