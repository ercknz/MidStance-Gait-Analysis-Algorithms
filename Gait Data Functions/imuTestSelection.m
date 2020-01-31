function path = imuTestSelection(testNum)
%% imuTestSelection
% This function is used to setup the paths for the IMU csv files. The
% file can vary depending on the test.
%
% function by erick nunez

%% Import Data
switch testNum
    case {101,102,103,104,105}
        % For tests 1-5 for subject 1
        path.RightAnkle = 'Shimmer IMU Files/2019-07-31_08.40.50_testright_SD_Session1/testright_Session1_Shimmer_CB32_Calibrated_SD.csv';
        path.LeftAnkle = 'Shimmer IMU Files/2019-07-31_08.40.30_Testwaist_SD_Session1/Testwaist_Session1_Shimmer_7696_Calibrated_SD.csv';
        path.RightWrist = 'Shimmer IMU Files/2019-07-31_08.40.12_DefaultTrial_SD_Session1/DefaultTrial_Session1_Shimmer_768F_Calibrated_SD.csv';
        path.LeftWrist = 'Shimmer IMU Files/2019-07-31_08.37.25_Testleft_SD_Session1/Testleft_Session1_Shimmer_7757_Calibrated_SD.csv';
        path.Waist = 'Shimmer IMU Files/2019-07-31_08.38.51_default_exp_SD_Session1/default_exp_Session1_Shimmer_7763_Calibrated_SD.csv';
    case {106,107,108,109,110,111}
        % For tests 6-11 for subject 1
        path.RightAnkle = 'Shimmer IMU Files/2019-11-14_10.19.09_testright_SD_Session1/testright_Session1_Shimmer_CB32_Calibrated_SD.csv';
        path.LeftAnkle = 'Shimmer IMU Files/2019-11-14_10.20.50_Testwaist_SD_Session1/Testwaist_Session1_Shimmer_7696_Calibrated_SD.csv';
        path.RightWrist = 'Shimmer IMU Files/2019-11-14_10.21.14_DefaultTrial_SD_Session1/DefaultTrial_Session1_Shimmer_768F_Calibrated_SD.csv';
        path.LeftWrist = 'Shimmer IMU Files/2019-11-14_10.18.38_Testleft_SD_Session1/Testleft_Session1_Shimmer_7757_Calibrated_SD.csv';
        path.Waist = 'Shimmer IMU Files/2019-11-14_10.21.39_DefaultTrial_SD_Session1/DefaultTrial_Session1_Shimmer_7763_Calibrated_SD.csv';
    case {201,202,203,204,205,206,207,208,209,210,211,212}
        % For tests 1-12 for subject 2
        path.RightAnkle = 'Shimmer IMU Files/2019-12-12_10.01.39_testright_SD_Session1/testright_Session1_Shimmer_CB32_Calibrated_SD.csv';
        path.LeftAnkle = 'Shimmer IMU Files/2019-12-12_10.00.39_Testwaist_SD_Session1/Testwaist_Session1_Shimmer_7696_Calibrated_SD.csv';
        path.RightWrist = 'Shimmer IMU Files/2019-12-12_10.01.22_DefaultTrial_SD_Session1/DefaultTrial_Session1_Shimmer_768F_Calibrated_SD.csv';
        path.LeftWrist = 'Shimmer IMU Files/2019-12-12_10.00.59_Testleft_SD_Session1/Testleft_Session1_Shimmer_7757_Calibrated_SD.csv';
        path.Waist = 'Shimmer IMU Files/2019-12-12_10.01.54_DefaultTrial_SD_Session1/DefaultTrial_Session1_Shimmer_7763_Calibrated_SD.csv';
    otherwise
        disp('no gait data found for subject')
end

end

