%% IMU syncing between sensors
% Script by Erick Nunez

%% clean up work space
clear; clc; close all;
%% Subject and session variables
% subjects 1, 2, 4, 5, 6, 7, 12, 13 have good data. 
subj.Num = 14;
looking = true;
sensorLoc = {'leftAnkle','rightAnkle','waist'}; 
% sensorLoc = {'leftAnkle','leftWrist','rightAnkle','rightWrist','waist'};
session.freqC = 20; session.filtOrder = 2;

%% Load Mat-file
load(['../../../Google Drive/School/PhD - Biomedical/Rutgers - Rotation/Gait Data/subjectMatLabData/subject',num2str(subj.Num),'RawDataLogs.mat']);

%% Session times
for i=1:length(sensorLoc)
    session.(sensorLoc{i}).startTime = raw.(sensorLoc{i})(1,1);
    session.(sensorLoc{i}).unixTimes = raw.(sensorLoc{i})(:,1);
    session.(sensorLoc{i}).times = session.(sensorLoc{i}).unixTimes - session.(sensorLoc{i}).startTime;
end

%% Sampling frequency and period
for i=1:length(sensorLoc)
    session.(sensorLoc{i}).period = 0;
    for j=2:length(session.(sensorLoc{i}).unixTimes)
        session.(sensorLoc{i}).period = session.(sensorLoc{i}).period + (session.(sensorLoc{i}).unixTimes(j) - session.(sensorLoc{i}).unixTimes(j-1));
    end
    session.(sensorLoc{i}).period  = session.(sensorLoc{i}).period/(j-1)/1000;
    session.(sensorLoc{i}).freqS  = 1/session.(sensorLoc{i}).period;
end

%% Butterworth low pass filter
session.freqS = 0; session.filtA = 0; session.filtB = 0;
for i=1:length(sensorLoc)
    session.freqS = session.freqS + session.(sensorLoc{i}).freqS/length(sensorLoc);
end
[session.filtB, session.filtA] = butter(session.filtOrder,session.freqC/(session.freqS/2));

%% Pulls session IMU data and finds resultant
for i = 1:length(sensorLoc)
    % Baseline IMU Data for session
    session.(sensorLoc{i}).acc = raw.(sensorLoc{i})(:,2:4);
    session.(sensorLoc{i}).gyro = raw.(sensorLoc{i})(:,6:8);
    session.(sensorLoc{i}).mag = raw.(sensorLoc{i})(:,9:11);
    session.(sensorLoc{i}).rAcc = [];
    % resultant Baseline IMU Data for session
    for j = 1:length(session.(sensorLoc{i}).acc)
        session.(sensorLoc{i}).rAcc(end+1) = sqrt(session.(sensorLoc{i}).acc(j,1)^2 + session.(sensorLoc{i}).acc(j,2)^2 + session.(sensorLoc{i}).acc(j,3)^2);
    end
end

%% Filters subject IMU data
for i = 1:length(sensorLoc)
    for j = 1:3
        % Filter trial IMU Data for subject
        session.(sensorLoc{i}).fAcc(:,j) = filtfilt(session.filtB, session.filtA, session.(sensorLoc{i}).acc(:,j));
        session.(sensorLoc{i}).fGyro(:,j) = filtfilt(session.filtB, session.filtA, session.(sensorLoc{i}).gyro(:,j));
        session.(sensorLoc{i}).fMag(:,j) = filtfilt(session.filtB, session.filtA, session.(sensorLoc{i}).mag(:,j));
    end
    session.(sensorLoc{i}).fRAcc = [];
    for j = 1:length(session.(sensorLoc{i}).acc)
        session.(sensorLoc{i}).fRAcc(end+1) = sqrt(session.(sensorLoc{i}).fAcc(j,1)^2 + session.(sensorLoc{i}).fAcc(j,2)^2 + session.(sensorLoc{i}).fAcc(j,3)^2);
    end
end

%% Algorithm to find sync between sensors
% for i = 1:length(sensorLoc)-1
%     data = findLongestZero(session.(sensorLoc{i}).fRAcc,session.(sensorLoc{i+1}).fRAcc); 
% end


%% View raw.
if looking
    for i = 1:length(sensorLoc)
        figure(151)
        subplot(4,1,i)
        plot(session.(sensorLoc{i}).times, session.(sensorLoc{i}).acc);
        grid on; xlabel('seconds'); ylabel('m/s^2'); legend('X','Y','Z');
        title(['Raw ',(sensorLoc{i}),' Accel'])
    end
    subplot(4,1,4)
    for i = 1:length(sensorLoc)
        plot(session.(sensorLoc{i}).times, session.(sensorLoc{i}).rAcc); hold on;
    end
    grid on; xlabel('seconds'); ylabel('m/s^2'); legend('waist','rightAnkle','leftAnkle');
    title(['Raw ',(sensorLoc{i}),' Resultant Accel'])
    
    for i = 1:length(sensorLoc)
        figure(111)
        subplot(4,1,i)
        plot(session.(sensorLoc{i}).times, session.(sensorLoc{i}).fAcc);
        grid on; xlabel('seconds'); ylabel('m/s^2'); legend('X','Y','Z');
        title(['Raw ',(sensorLoc{i}),' Filtered Accel'])
    end
    subplot(4,1,4)
    for i = 1:length(sensorLoc)
        plot(session.(sensorLoc{i}).times, session.(sensorLoc{i}).fRAcc); hold on;
    end
    grid on; xlabel('seconds'); ylabel('m/s^2'); legend('waist','rightAnkle','leftAnkle');
    title(['Raw ',(sensorLoc{i}),' Filtered Resultant Accel'])
end
