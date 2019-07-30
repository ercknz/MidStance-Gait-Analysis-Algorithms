%% IMU Data Rotation
% Script by Erick Nunez

%% clean up work space
clear;
%% Subject and session variables
subj.Num = 15;
sensorLoc = {'leftAnkle','leftWrist','waist'};
dataIndex = [2,4,6,8,9,11]';

%% Load Mat-file
load(['../../../Google Drive/School/PhD - Biomedical/Rutgers - Rotation/Gait Data/subjectMatLabData/subject',num2str(subj.Num),'RawDataLogs.mat']);

%% rotate data
for i = 1:length(sensorLoc)
    for index = dataIndex
        newData = raw.(sensorLoc{i})(:,index) * -1;
        raw.(sensorLoc{i})(:,index) = newData;
    end
end

%% Save modifed mat-file
save(['../../../Google Drive/School/PhD - Biomedical/Rutgers - Rotation/Gait Data/subjectMatLabData-Rot/subject',num2str(subj.Num),'RawDataLogs-Rot.mat'],'raw','logFile','subj');