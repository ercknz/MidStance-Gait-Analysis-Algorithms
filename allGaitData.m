%% All Gait data 
% This fucntion uses imuFusionGlobal.m and segmentGaitData.m to find the
% global accelerations and split the data into gait cycles. Then it says
% the data or plot into the appropriate folder or worksheet. Plots are
% saved in folder IMU Data Images and data is saved in
% cycleSeperationData.xlsx. Both of which should be in the same folder as
% the functions. 
%
% Script by erick nunez

%% clean up work space
clear; clc; close all;

%% Setting up varibles to iterate over.
% subjects 1, 2, 4, 5, 6, 7, 12, 13 have good gait data.
goodSubj = [1, 2, 4, 5];
timesSigma = [2.5, 3];
yAxis = [true, false];

%% Run Main functions
for i = goodSubj
    for j = timesSigma
        for k = yAxis
            imuFusionGlobal(i,j,k,false)
        end
    end
end