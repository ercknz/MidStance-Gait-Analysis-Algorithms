%% Shimmer Vicon test data import using classes
% This script use the classes used to organize the data collected from the
% shimmer-vicon test. This script is used to test the data stored in the
% classes and plot them.
%
% script by erick nunez

%% clean up
clear; clc; close all;

%% File Paths
addpath('Gait Data Classes');
addpath('Gait Data Functions');
path.RightAnkle = 'Shimmer IMU Files/2019-07-31_08.40.50_testright_SD_Session1/testright_Session1_Shimmer_CB32_Calibrated_SD.csv';
path.LeftAnkle = 'Shimmer IMU Files/2019-07-31_08.40.30_Testwaist_SD_Session1/Testwaist_Session1_Shimmer_7696_Calibrated_SD.csv';
path.RightWrist = 'Shimmer IMU Files/2019-07-31_08.40.12_DefaultTrial_SD_Session1/DefaultTrial_Session1_Shimmer_768F_Calibrated_SD.csv';
path.LeftWrist = 'Shimmer IMU Files/2019-07-31_08.37.25_Testleft_SD_Session1/Testleft_Session1_Shimmer_7757_Calibrated_SD.csv';
path.Waist = 'Shimmer IMU Files/2019-07-31_08.38.51_default_exp_SD_Session1/default_exp_Session1_Shimmer_7763_Calibrated_SD.csv';
logFileName = 'shimmerViconTest2Times.xlsx';

%% Variables
% 5 shimmer-vicon tests performed.
subjNum = 1;    method = 2;
compHighPass = 0.98;
accRange = 2:4;         gyroRange = 9:11;       magRange = 12:14;
cutOffFreq = 15;        filterOrder = 4;
sigma = 3;

%% Create measured data object and pull data
mData = shimmerGaitData(subjNum,logFileName,true);
mData.addSensor('rightAnkle',path.RightAnkle);
mData.addSensor('leftAnkle',path.LeftAnkle);
mData.addSensor('waist',path.Waist);
mData.getAcc(accRange);
mData.getGyro(gyroRange);
mData.getMag(magRange);

%% Filter data
mData.getSamplingFreq();
mData.createFilter(cutOffFreq, filterOrder);
mData.filterData();

%% Rotate sensors
% Rotate sensors to match the following orientation:
% +x: direction of walking
% +y: medial-to-lateral
% +z: lower-to-upper body
mData.waist.rotateIMU('y', 90);
mData.waist.rotateIMU('x',-90);
mData.rightAnkle.rotateIMU('x',-90);
mData.leftAnkle.rotateIMU('y', 180);
mData.leftAnkle.rotateIMU('x',-90);

%% Finds Angles using complimentary filter
waist.preAngles = compFusionAngles(compHighPass, mData.waist.fPreAcc, mData.waist.fPreGyro, mData.waist.fPreMag, mData.waist.sampleFreq);
waist.angles = compFusionAngles(compHighPass, mData.waist.fAcc, mData.waist.fGyro, mData.waist.fMag, mData.waist.sampleFreq);

rAnkle.preAngles = compFusionAngles(compHighPass, mData.rightAnkle.fPreAcc, mData.rightAnkle.fPreGyro, mData.rightAnkle.fPreMag, mData.rightAnkle.sampleFreq);
rAnkle.angles = compFusionAngles(compHighPass, mData.rightAnkle.fAcc, mData.rightAnkle.fGyro, mData.rightAnkle.fMag, mData.rightAnkle.sampleFreq);

lAnkle.preAngles = compFusionAngles(compHighPass, mData.leftAnkle.fPreAcc, mData.leftAnkle.fPreGyro, mData.leftAnkle.fPreMag, mData.leftAnkle.sampleFreq);
lAnkle.angles = compFusionAngles(compHighPass, mData.leftAnkle.fAcc, mData.leftAnkle.fGyro, mData.leftAnkle.fMag, mData.leftAnkle.sampleFreq);

%% Calculate Global Accelerations and removes baseline
% +x: direction of walking
% +y: medial-to-lateral
% +z: lower-to-upper body
waist.preGlobalAcc = imuGlobalAcc(waist.preAngles, mData.waist.fPreAcc);
waist.globalAcc = imuGlobalAcc(waist.angles, mData.waist.fAcc);
waist.globalAcc = waist.globalAcc - mean(waist.preGlobalAcc);

rAnkle.preGlobalAcc = imuGlobalAcc(rAnkle.preAngles, mData.rightAnkle.fPreAcc);
rAnkle.globalAcc = imuGlobalAcc(rAnkle.angles, mData.rightAnkle.fAcc);
rAnkle.globalAcc = rAnkle.globalAcc - mean(rAnkle.preGlobalAcc);

lAnkle.preGlobalAcc = imuGlobalAcc(lAnkle.preAngles, mData.leftAnkle.fPreAcc);
lAnkle.globalAcc = imuGlobalAcc(lAnkle.angles, mData.leftAnkle.fAcc);
lAnkle.globalAcc = lAnkle.globalAcc - mean(lAnkle.preGlobalAcc);

%% Segmentation of Data
rAnkle.steps = findStepIndexes(mData.rightAnkle.fGyro, mData.rightAnkle.times, sigma);
lAnkle.steps = findStepIndexes(mData.leftAnkle.fGyro, mData.leftAnkle.times, sigma);

%% Method #1
if method == 1
    % Data Resampling
    rAnkle.Step = imuDataResampling(rAnkle.steps.indexes, mData.rightAnkle.times, rAnkle.globalAcc, mData.rightAnkle.fGyro, rAnkle.angles);
    lAnkle.Step = imuDataResampling(lAnkle.steps.indexes, mData.leftAnkle.times, lAnkle.globalAcc, mData.leftAnkle.fGyro, lAnkle.angles);

    % Walking Speed
    [waist, rAnkle, lAnkle] = gaitSpeed(subjNum, mData.waist, waist, rAnkle, lAnkle);

    % Saves the workspace
    % save(['Shimmer Workspaces/shimmerViconIMU',num2str(subjNum),'.mat'])

    % Plot the calculated data
    plotStepIMUdata(subjNum, 'RightAnkle', rAnkle.steps.HSindexes, ...
        rAnkle.globalAcc, mData.rightAnkle.fGyro, rAnkle.angles, mData.rightAnkle.times, ...
        rAnkle.steps.indexes, rAnkle.Step.meanAcc, rAnkle.Step.meanGyro, rAnkle.Step.meanAngle, ...
        rAnkle.Step.stdAcc, rAnkle.Step.stdGyro, rAnkle.Step.stdAngle, rAnkle.Step.avgTimes)
    plotStepIMUdata(subjNum, 'LeftAnkle', lAnkle.steps.HSindexes, ...
        lAnkle.globalAcc, mData.leftAnkle.fGyro, lAnkle.angles, mData.leftAnkle.times, ...
        lAnkle.steps.indexes, lAnkle.Step.meanAcc, lAnkle.Step.meanGyro, lAnkle.Step.meanAngle, ...
        lAnkle.Step.stdAcc, lAnkle.Step.stdGyro, lAnkle.Step.stdAngle, lAnkle.Step.avgTimes)
end

%% Method #2
if method == 2
    % Speeds per step
    rAnkle.steps = cycleIMUdata(rAnkle.steps, rAnkle.globalAcc, mData.rightAnkle.times);
    lAnkle.steps = cycleIMUdata(lAnkle.steps,lAnkle.globalAcc, mData.leftAnkle.times);
    if rAnkle.steps.indexes{1}(1) < lAnkle.steps.indexes{1}(1)
        waist.steps = cycleIMUdata(rAnkle.steps, waist.globalAcc, mData.waist.times);
    else
        waist.steps = cycleIMUdata(lAnkle.steps, waist.globalAcc, mData.waist.times);
    end

    % Saves data
    writematrix([mean(waist.steps.Speed), mean(rAnkle.steps.Speed), mean(lAnkle.steps.Speed)], 'Shimmer Vicon Results/shimmerViconVel.xlsx', 'Range', ['K',num2str(subjNum+2)])
end

