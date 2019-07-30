function imuDataImport(subject)
%% imuDataImport Function
% This function is used to import data into the matlab workspace from the
% csv files depending on the subject. Data should be in the specific
% folders and subject data should be updated in the
% "rewritten-completeDataLog" file. 
% 
% function by erick nunez

%% Import Data
subj.Num = subject;
logFile = readmatrix('rewritten-CompleteDataLog.xlsx');
switch subj.Num
    case {1,2,3}
        % For subjects 1, 2 and 3 use 2.13.19 data:
        data.leftAnkle = readmatrix('2.13.19/2019-02-13LeftAnkle/default_exp_Session1_LeftAnkle_Calibrated_SD.csv');
        data.leftWrist = readmatrix('2.13.19/2019-02-13LeftWrist/default_exp_Session2_LeftWrist_Calibrated_SD.csv');
        data.rightAnkle = readmatrix('2.13.19/2019-02-13RightAnkle/default_exp_Session2_RightAnkle_Calibrated_SD.csv');
        data.rightWrist = readmatrix('2.13.19/2019-02-13RightWrist/default_exp_Session2_RightWrist_Calibrated_SD.csv');
        data.waist = readmatrix('2.13.19/2019-02-13Waist/default_exp_Session2_Waist_Calibrated_SD.csv');
    case {4,5,6}
        % For subjects 4, 5 and 6 use 2.14.19 data:
        data.leftAnkle = readmatrix('2.14.19/2019-02-14_LeftAnkle/default_exp_Session1_LeftAnkle_Calibrated_SD.csv');
        data.leftWrist = readmatrix('2.14.19/2019-02-14_LeftWrist/default_exp_Session1_LeftWrist_Calibrated_SD.csv');
        data.rightAnkle = readmatrix('2.14.19/2019-02-14_RightAnkle/default_exp_Session1_RightAnkle_Calibrated_SD.csv');
        data.rightWrist = readmatrix('2.14.19/2019-02-14_RightWrist/default_exp_Session1_RightWrist_Calibrated_SD.csv');
        data.waist = readmatrix('2.14.19/2019-02-14_Waist/default_exp_Session1_Waist_Calibrated_SD.csv');
    case {7,8}
        % For subjects 7 and 8 use 2.18.19 data:
        data.leftAnkle = readmatrix('2.18.19/2019-02-18_LeftAnkle/default_exp_Session1_LeftAnkle_Calibrated_SD.csv');
        data.leftWrist = readmatrix('2.18.19/2019-02-18_LeftWrist/default_exp_Session1_LeftWrist_Calibrated_SD.csv');
        data.rightAnkle = readmatrix('2.18.19/2019-02-18_RightAnkle/default_exp_Session1_RightAnkle_Calibrated_SD.csv');
        data.rightWrist = readmatrix('2.18.19/2019-02-18_RightWrist/default_exp_Session1_RightWrist_Calibrated_SD.csv');
        data.waist = readmatrix('2.18.19/2019-02-18_Waist/default_exp_Session1_Waist_Calibrated_SD.csv');
    case {9,10,11}
        % For subjects 9, 10 and 11 use 2.19.19
        data.leftAnkle = readmatrix('2.19.19/2019-02-19_LeftAnkle/default_exp_Session1_LeftAnkle_Calibrated_SD.csv');
        data.leftWrist = readmatrix('2.19.19/2019-02-19_LeftWrist/default_exp_Session1_LeftWrist_Calibrated_SD.csv');
        data.rightAnkle = readmatrix('2.19.19/2019-02-19_RightAnkle/default_exp_Session1_RightAnkle_Calibrated_SD.csv');
        data.rightWrist = readmatrix('2.19.19/2019-02-19_RightWrist/default_exp_Session1_RightWrist_Calibrated_SD.csv');
        data.waist = readmatrix('2.19.19/2019-02-19_Waist/default_exp_Session1_Waist_Calibrated_SD.csv');
    case {12,13}
        % For subjects 12 and 13 use 3.27.19
        data.leftAnkle = readmatrix('3.27.19/LeftAnkle/default_exp_Session2_LeftAnkle_Calibrated_SD.csv');
        data.leftWrist = readmatrix('3.27.19/LeftWrist/default_exp_Session2_LeftWrist_Calibrated_SD.csv');
        data.rightAnkle = readmatrix('3.27.19/RightAnkle/default_exp_Session2_RightAnkle_Calibrated_SD.csv');
        data.rightWrist = readmatrix('3.27.19/RightWrist/default_exp_Session1_RightWrist_Calibrated_SD.csv');
        data.waist = readmatrix('3.27.19/Waist/default_exp_Session2_Waist_Calibrated_SD.csv');
    case {14,15}
        % For subjects 14 and 15 use 3.29.19
        data.leftAnkle = readmatrix('3.29.19/LeftAnkle/default_exp_Session1_LeftAnkle_Calibrated_SD.csv');
        data.leftWrist = readmatrix('3.29.19/LeftWrist/default_exp_Session1_LeftWrist_Calibrated_SD.csv');
        data.rightAnkle = readmatrix('3.29.19/RightAnkle/default_exp_Session1_RightAnkle_Calibrated_SD.csv');
        data.rightWrist = readmatrix('3.29.19/RightWrist/default_exp_Session1_RightWrist_Calibrated_SD.csv');
        data.waist = readmatrix('3.29.19/Waist/default_exp_Session1_Waist_Calibrated_SD.csv');
    case {100}
        % For subjects 100 use 4.15.19
        data.leftAnkle = readmatrix('4.15.19/LeftAnkle/default_exp_Session2_LeftAnkle_Calibrated_SD.csv');
        data.leftWrist = readmatrix('4.15.19/LeftWrist/default_exp_Session2_LeftWrist_Calibrated_SD.csv');
        data.rightAnkle = readmatrix('4.15.19/RightAnkle/default_exp_Session2_RightAnkle_Calibrated_SD.csv');
        data.rightWrist = readmatrix('4.15.19/RightWrist/default_exp_Session2_RightWrist_Calibrated_SD.csv');
        data.waist = readmatrix('4.15.19/Waist/default_exp_Session2_Waist_Calibrated_SD.csv');
    otherwise
        disp('no gait data found for subject')
end

%% clean data
sensorLoc = {'leftAnkle','leftWrist','rightAnkle','rightWrist','waist'};
for i = 1:length(sensorLoc)
    for j = 1:11
        raw.(sensorLoc{i})(:,j) = data.(sensorLoc{i})(~isnan(data.(sensorLoc{i})(:,j)),j);
    end
end

%% Save mat-file
save(['subjectMatLabData/subject',num2str(subj.Num),'RawDataLogs.mat'],'raw','logFile','subj');
end

