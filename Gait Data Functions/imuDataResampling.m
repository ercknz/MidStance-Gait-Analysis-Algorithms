function avgStep = imuDataResampling(indexes, times, accXYZ, gyroXYZ, anglesXYZ)
%% IMU Data Resampling
% This function takes in IMU data and calculated angles. Then using the
% steps indexes to resample the data. This is used to normalized the data
% based on heel strike to heel strike. StepsCell should contain the step
% indexes.  Gyro data and angles are optional arguments.
%
% The function assumes the following orientation of the data:
% +x: direction of walking
% +y: medial to Left Side
% +z: normal to ground (vertical)
%
% Function by Erick Nunez

%% Variables to be used
stepsAcc = cell(size(indexes));        stepsTimes = {cell(size(indexes))};
if nargin == 4 || nargin == 5
    stepsGyro = cell(size(indexes));
end
if nargin == 5
    stepsAngle = cell(size(indexes));
end

%% Seperates data into Steps
for i = 1:length(indexes)
    stepsTimes{i} = times(indexes{i});
    stepsAcc{i} = accXYZ(indexes{i},:);
    if nargin == 4 || nargin == 5
        stepsGyro{i} = gyroXYZ(indexes{i},:);
    end
    if nargin == 5
        stepsAngle{i} = anglesXYZ(indexes{i},:);
    end
end

%% Finds mean and std of step frame length and mean times
for i = 1:length(indexes)
    frameLengths(i) = length(indexes{i});
    stepDurations(i) = stepsTimes{i}(end) - stepsTimes{i}(1);
end
meanFrames = round(mean(frameLengths));     stdFrames = std(frameLengths);
meanDuration = mean(stepDurations);         stdDurations = std(stepDurations);

%% Resamples the IMU data
for i = 1:length(indexes)
    stepsAcc{i} = resample(stepsAcc{i}, meanFrames, length(indexes{i}));
    if nargin == 4 || nargin == 5
        stepsGyro{i} = resample(stepsGyro{i}, meanFrames, length(indexes{i}));
    end
    if nargin == 5
        stepsAngle{i} = resample(stepsAngle{i}, meanFrames, length(indexes{i}));
    end
end
avgStep.avgTimes = linspace(0,meanDuration,meanFrames);

%% Finds average step profile
[avgStep.meanAcc, avgStep.stdAcc] = findProfile(meanFrames, stepsAcc);
if nargin == 4 || nargin == 5
    [avgStep.meanGyro, avgStep.stdGyro] = findProfile(meanFrames, stepsGyro);
end
if nargin == 5
    [avgStep.meanAngle, avgStep.stdAngle] = findProfile(meanFrames, stepsAngle);
end
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