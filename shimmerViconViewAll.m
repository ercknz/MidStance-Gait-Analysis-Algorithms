%% View All data of Shimmer-Vicon test
% Script by Erick Nunez

%% clean up work space
clear; clc; close all;
%% Subject and session variables
looking = true;
% sensorLoc = {'leftAnkle','rightAnkle','waist'}; 
sensorLoc = {'leftAnkle','leftWrist','rightAnkle','rightWrist','waist'};

%% Load Mat-file
load('../../../Google Drive/School/PhD - Biomedical/Rutgers - Rotation/Gait Data/Shimmer-Vicon Test/shimmerViconTestRawData.mat');

%% Test times, frequency and period
raw.freqC = 20; raw.filtOrder = 2;
for i=1:length(sensorLoc)
    raw.(sensorLoc{i}).unixTimes = raw.(sensorLoc{i}).IMUdata(:,1);
    raw.(sensorLoc{i}).startTime = raw.(sensorLoc{i}).IMUdata(1,1);
    raw.(sensorLoc{i}).times = raw.(sensorLoc{i}).unixTimes - raw.(sensorLoc{i}).startTime;
    raw.(sensorLoc{i}).period = 0;
    for j=2:length(raw.(sensorLoc{i}).unixTimes)
        raw.(sensorLoc{i}).period = raw.(sensorLoc{i}).period + (raw.(sensorLoc{i}).unixTimes(j) - raw.(sensorLoc{i}).unixTimes(j-1));
    end
    raw.(sensorLoc{i}).period  = raw.(sensorLoc{i}).period/(j-1)/1000;
    raw.(sensorLoc{i}).freqS  = 1/raw.(sensorLoc{i}).period;
end

%% Butterworth low pass filter
raw.freqS = 0; raw.filtA = 0; raw.filtB = 0;
for i=1:length(sensorLoc)
    raw.freqS = raw.freqS + raw.(sensorLoc{i}).freqS/length(sensorLoc);
end
[raw.filtB, raw.filtA] = butter(raw.filtOrder,raw.freqC/(raw.freqS/2));

%% Pulls IMU data and applies filter
for i = 1:length(sensorLoc)
    % IMU Data for test day
    raw.(sensorLoc{i}).acc = raw.(sensorLoc{i}).IMUdata(:,2:4);
    raw.(sensorLoc{i}).gyro = raw.(sensorLoc{i}).IMUdata(:,6:8);
    raw.(sensorLoc{i}).mag = raw.(sensorLoc{i}).IMUdata(:,9:11);
    for j = 1:3
        % Filter trial IMU Data for subject
        raw.(sensorLoc{i}).fAcc(:,j) = filtfilt(raw.filtB, raw.filtA, raw.(sensorLoc{i}).acc(:,j));
        raw.(sensorLoc{i}).fGyro(:,j) = filtfilt(raw.filtB, raw.filtA, raw.(sensorLoc{i}).gyro(:,j));
        raw.(sensorLoc{i}).fMag(:,j) = filtfilt(raw.filtB, raw.filtA, raw.(sensorLoc{i}).mag(:,j));
    end
end

%% View raw.
if looking
    for i = 1:length(sensorLoc)
        figure(i)
        subplot(3,1,1)
        plot(raw.(sensorLoc{i}).times, raw.(sensorLoc{i}).fAcc);
        grid on; xlabel('seconds'); ylabel('m/s^2'); legend('X','Y','Z');
        title(['Filtered Raw ',sensorLoc{i},' Shimmer-Vicon Accel'])
        subplot(3,1,2)
        plot(raw.(sensorLoc{i}).times, raw.(sensorLoc{i}).fGyro);
        grid on; xlabel('seconds'); ylabel('degrees/sec'); legend('X','Y','Z');
        title(['Filtered Raw ',sensorLoc{i},' Shimmer-Vicon Gyro'])
        subplot(3,1,3)
        plot(raw.(sensorLoc{i}).times, raw.(sensorLoc{i}).fMag);
        grid on; xlabel('seconds'); ylabel('localFlux'); legend('X','Y','Z');
        title(['Filtered Raw ',sensorLoc{i},' Shimmer-Vicon Mag'])
    end
end