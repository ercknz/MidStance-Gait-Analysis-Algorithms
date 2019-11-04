function steps = findMSIndexes(anglesXYZ, times, sigma)
%% Finds MS Indexes
% This function is used to find the Mid Stance point for each step based on 
% the angles. The function assumes the following orientation of
% the data:
% +x: direction of walking
% +y: medial to Left Side
% +z: normal to ground (vertical)
% 
% Function by Erick Nunez

%% Variables to be used
indexes = {};       removed = 0;

%% Finds Mins
MSindex = find(islocalmin(abs(anglesXYZ(:,2)), 'MinSeparation', 0.2, 'SamplePoints',times));

%% Find rate of change of angles to find slope
angRateY = gradient(anglesXYZ(:,2), mean(diff(times)));
for i = 1:(length(MSindex) - removed)
    if angRateY(MSindex(i - removed)) < 0
       MSindex(i- removed) = [];
       removed = removed + 1;
    end
end

%% separates steps
for i =2:length(MSindex)
    indexes{end+1} = MSindex(i-1):MSindex(i)-1;
end

%% Finds mean and std for steps
for i = 1:length(indexes)
    frameLengths(i) = length(indexes{i});
end
meanFrames = mean(frameLengths);    stdFrames = std(frameLengths);

%% Removes outliers in steps
removed = 0;
for i = 1:length(indexes)
    if (length(indexes{i-removed}) < meanFrames - stdFrames*sigma) || (length(indexes{i-removed}) > meanFrames + stdFrames*sigma)
        indexes(i-removed) = [];
        removed = removed + 1;
    end
end

%% Save data for output
steps.removed = removed;            steps.indexes = indexes;
steps.meanFrames = meanFrames;      steps.stdFrames = stdFrames;
steps.MSindexes = MSindex;

