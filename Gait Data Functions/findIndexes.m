function steps = findIndexes(gyroXYZ, times, sigma, removeY)
%% Finds All Indexes
% This function is used to find the Mid Stance point for each step based on
% the gyroscope data. It find the local maxs first which are the HS and TO
% indexes. Then it find the midpoint between them which ends up being
% midstance.
% +x: direction of walking
% +y: medial to Left Side
% +z: normal to ground (vertical)
%
% Function by Erick Nunez

%% Variables to be used
indexes = {};

%% Finds Mins
idx = find(islocalmin(gyroXYZ(:,2), 'MinSeparation', 0.4, 'SamplePoints',times));
meanMins = mean(gyroXYZ(idx,2));

%% Finds MidSwing(MSw)
removed = 0;
for i = 1:length(idx)
    if gyroXYZ(idx(i-removed),2) >  meanMins
        idx(i-removed) = [];
        removed = removed + 1;
    end
end
MSwIdx = idx;

%% Finds ToeOff(TO)
TOIdx = nan(1,length(MSwIdx));
for i = 1:length(MSwIdx)
    for j = 2:30
        if MSwIdx(i)-1-j>0
            if gyroXYZ(MSwIdx(i)-1-j,2) < gyroXYZ(MSwIdx(i)-j,2) && gyroXYZ(MSwIdx(i)+1-j,2) < gyroXYZ(MSwIdx(i)-j,2)
                TOIdx(i) = MSwIdx(i)-j;
                break
            end
        end
    end
end
TOIdx = TOIdx(~isnan(TOIdx));

%% Finds HeelStrike(HS)
HSIdx = nan(1,length(MSwIdx));
for i = 1:length(MSwIdx)
    for j = 2:30
        if MSwIdx(i)+1+j<=length(gyroXYZ(:,2))
            if gyroXYZ(MSwIdx(i)-1+j,2) < gyroXYZ(MSwIdx(i)+j,2) && gyroXYZ(MSwIdx(i)+1+j,2) < gyroXYZ(MSwIdx(i)+j,2)
                HSIdx(i) = MSwIdx(i)+j;
                break
            end
        end
    end
end
HSIdx = HSIdx(~isnan(HSIdx));

%% Finds MidStance before ToeOff(MStTOIdx)
MStTOIdx = nan(1,length(TOIdx));
for i = 1:length(TOIdx)
    for j = 2:30
        if TOIdx(i)-1-j > 0
            if gyroXYZ(TOIdx(i)-1-j,2) > gyroXYZ(TOIdx(i)-j,2) && gyroXYZ(TOIdx(i)+1-j,2) > gyroXYZ(TOIdx(i)-j,2)
                MStTOIdx(i) = TOIdx(i)-j;
                break
            end
        end
    end
end
MStTOIdx = MStTOIdx(~isnan(MStTOIdx));

%% Finds MidStance(MSt) after HeelStrike(MStHSIdx)
MStHSIdx = nan(1,length(HSIdx));
for i = 1:length(HSIdx)
    for j = 2:30
        if HSIdx(i)+1+j <= length(gyroXYZ(:,2))
            if gyroXYZ(HSIdx(i)-1+j,2) > gyroXYZ(HSIdx(i)+j,2) && gyroXYZ(HSIdx(i)+1+j,2) > gyroXYZ(HSIdx(i)+j,2)
                MStHSIdx(i) = HSIdx(i)+j;
                break
            end
        end
    end
end
MStHSIdx = MStHSIdx(~isnan(MStHSIdx));

%% Combining Estimated MidStances(MStIdx)
idx = unique([MStHSIdx, MStTOIdx],'sorted');
MStIdx = [];    i = 1;
while i < length(idx)
    if abs(idx(i) - idx(i+1)) <= 15
        MStIdx = [MStIdx,round((idx(i)+idx(i+1))/2)];
        i = i + 2;
    else
        MStIdx = [MStIdx,idx(i)];
        i = i + 1;
    end
end
if abs(idx(end) - idx(end-1)) > 15
    MStIdx = [MStIdx,idx(end)];
end

%% Separates steps
for i =1:length(MStIdx)-1
    idx = MStIdx(i):MStIdx(i+1)-1;
    if ~isempty(intersect(idx,MSwIdx))
        indexes{end+1} = idx;
    end
end

%% Finds mean and std for steps
for i = 1:length(indexes)
    frameLengths(i) = length(indexes{i});
end
meanFrames = round(mean(frameLengths));    stdFrames = round(std(frameLengths));

%% Removes outliers in steps
removed = 0;
for i = 1:length(indexes)
    if (length(indexes{i-removed}) < meanFrames - stdFrames*sigma) || (length(indexes{i-removed}) > meanFrames + stdFrames*sigma)
        indexes(i-removed) = [];
        removed = removed + 1;
    end
end

%% Removes Outliers in Y axis using resampled data
if removeY
    % Resamples to normalize MSt to MSt
    ReSampGyroSteps = cell(1,length(indexes));
    for i = 1:length(indexes)
        ReSampGyroSteps{i} = resample(gyroXYZ(indexes{i},2), meanFrames, length(indexes{i}));
    end
    % Finds Average Profile
    [meanGyro, stdGyro] = findProfile(meanFrames, ReSampGyroSteps);
    % Removes Y outliers
    removed = 0;
    for i = 1:length(ReSampGyroSteps)
        for j = 1:meanFrames
            if (ReSampGyroSteps{i-removed}(j) < (meanGyro(j) - sigma * stdGyro(j))) || (ReSampGyroSteps{i-removed}(j) > (meanGyro(j) + sigma * stdGyro(j)))
                indexes(i-removed) = [];
                ReSampGyroSteps(i-removed) = [];
                removed = removed + 1;
                break
            end
        end
    end
end

%% Save data for output
steps.removed = removed;            steps.indexes = indexes;
steps.meanFrames = meanFrames;      steps.stdFrames = stdFrames;
steps.MStIdx = MStIdx;              steps.MSwIdx = MSwIdx;
steps.HSIdx = HSIdx;                steps.TOIdx = TOIdx;

end

%% Find Profile function
% This function finds the mean profile for the data and standard deviation
% of the XYZ data provided.
%
% Function by erick nunez
function [meanProfile, stdProfile] = findProfile(frames, dataXYZcell)
for i = 1:frames
    samples = [];
    for j = 1:length(dataXYZcell)
        samples(end+1,:) = dataXYZcell{j}(i,:);
    end
    meanProfile(i,:) = mean(samples);   stdProfile(i,:) = std(samples);
end
end

