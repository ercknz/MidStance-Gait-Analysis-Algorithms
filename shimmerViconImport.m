function shimmerViconImport()
%% shimmerViconImport Function
% This function is used to import the data from the shimmer-vicon test and
% generate mat files to work within matlab.
%
% function by erick nunez

%% Import Data
% For subjects 1, 2 and 3 use 2.13.19 data:
data.rightAnkle = readmatrix('Shimmer-Vicon Test/2019-07-19_08.51.36_DefaultTrial_SD_Session1-77DC/DefaultTrial_Session1_Shimmer_77DC_Calibrated_SD.csv');
data.leftAnkle = readmatrix('Shimmer-Vicon Test/2019-07-19_08.53.18_default_exp_SD_Session1-7763/default_exp_Session1_Shimmer_7763_Calibrated_SD.csv');
data.rightWrist = readmatrix('Shimmer-Vicon Test/2019-07-19_08.55.29_default_exp_SD_Session1-CB32/testright_Session1_Shimmer_CB32_Calibrated_SD.csv');
data.waist = readmatrix('Shimmer-Vicon Test/2019-07-19_08.57.06_default_exp_SD_Session1-7757/Testleft_Session1_Shimmer_7757_Calibrated_SD.csv');
data.leftWrist = readmatrix('Shimmer-Vicon Test/2019-07-19_08.59.01_default_exp_SD_Session1-7696/Testwaist_Session1_Shimmer_7696_Calibrated_SD.csv');

%% clean data
sensorLoc = {'leftAnkle','leftWrist','rightAnkle','rightWrist','waist'};
for i = 1:length(sensorLoc)
    for j = 1:11
        raw.(sensorLoc{i}).IMUdata(:,j) = data.(sensorLoc{i})(~isnan(data.(sensorLoc{i})(:,j)),j);
    end
end

%% Save mat-file
save('Shimmer-Vicon Test/shimmerViconTestRawData.mat','raw');
end

