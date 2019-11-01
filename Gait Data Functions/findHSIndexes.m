function steps = findHSIndexes(gyroXYZ, times, sigma)
%% Finds Step Indexes
% This function is used to find the maxs in the gyro data. Function expects
% gyro data and will use it to find local max to separate into Heel Strike
% - Heel Strike segments. The function assumes the following orientation of
% the data:
% +x: direction of walking
% +y: medial-to-lateral
% +z: lower-to-upper body
% 
% Function by Erick Nunez

%% Variables to be used
HSindex = [];       indexes = {};

%% Finds Mins
maxIndex = find(islocalmax(gyroXYZ(:,2), 'MinSeparation', 0.4, 'SamplePoints',times));
maxValues = gyroXYZ(maxIndex,2);

%% Finds Heel Strikes
for i = 2:length(maxIndex)-1
   if maxValues(i) > maxValues(i-1) && maxValues(i) > maxValues(i+1)
      HSindex(end+1) = maxIndex(i); 
   end
end

%% Sepatates steps
for i =2:length(HSindex)
    indexes{end+1} = HSindex(i-1):HSindex(i)-1;
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
steps.HSindexes = HSindex;

end
