function ShimmerViconTestExampleCode(testNum, method, compHighPass, height)
%% Shimmer Vicon test data import using classes
% This script use the classes used to organize the data collected from the
% shimmer-vicon test. This script is used to test the data stored in the
% classes and plot them. 
%
% script by erick nunez

%% clean up
% clear; 
% clc; 
close all;
 
%% Variables
% 5 shimmer-vicon tests performed. 
% testNum = 212;  method = 10;    compHighPass = 0.97;  height = 1.778;
accRange = 2:4; % Low Noise
% accRange = 5:7; % Wide Range           
gyroRange = 9:11;       magRange = 12:14;
gyroCutOff = 5;        accCutOff = [0.25 5];  filterOrder = 4;
sigma = 3;

%% File Paths
addpath('Gait Data Classes');
addpath('Gait Data Functions');
path = imuTestSelection(testNum);
logFileName = 'shimmerViconTest2Times.xlsx';

%% Create measured data object and pull data
mData = shimmerGaitData(testNum,logFileName,true);
mData.addSensor('rightAnkle',path.RightAnkle);
mData.addSensor('leftAnkle',path.LeftAnkle);
mData.addSensor('waist',path.Waist);
mData.getAcc(accRange);
mData.getGyro(gyroRange);
mData.getMag(magRange);

%% Filter data 
mData.getSamplingFreq();
mData.createFilter(gyroCutOff, filterOrder);
mData.filterData('gyro');
mData.createFilter(accCutOff, filterOrder);
mData.filterData('acc');
mData.createFilter(accCutOff, filterOrder);
mData.filterData('mag');

%% Method N
if method == 'N'
    fig0 = figure;
    subplot(2,1,1)
    plot(mData.rightAnkle.times, mData.rightAnkle.fGyro); grid on; xlim([20,40]);
    legend('X','Y','Z'); title(['Test ',num2str(testNum),' Right Ankle Filtered Gyro Before Rotating']);
    xlabel('Time (sec)'); ylabel('Angular Velocities (deg/sec)'); 
    subplot(2,1,2)
    plot(mData.leftAnkle.times, mData.leftAnkle.fGyro); grid on; xlim([20,40]);
    legend('X','Y','Z'); title(['Test ',num2str(testNum),' Left Ankle Filtered Gyro Before Rotating']);
    xlabel('Time (sec)'); ylabel('Angular Velocities (deg/sec)'); 
    set(fig0, 'PaperOrientation', 'landscape', 'PaperUnits', 'normalized', 'PaperPosition',[0,0,1,1]);
    saveas(fig0, ['Shimmer Vicon Results/subj',num2str(testNum),'PreRotationGyro.pdf']);
end

%% Rotate sensors
% Rotate sensors to match the following orientation:
% +x: direction of walking
% +y: medial to Left Side
% +z: normal to ground (vertical)
mData.waist.rotateIMU('y', 90);
mData.waist.rotateIMU('x',-90);
mData.rightAnkle.rotateIMU('x',-90);
mData.leftAnkle.rotateIMU('y', 180);
mData.leftAnkle.rotateIMU('x',-90);

%% Method R
if method == 'R'
    fig0 = figure;
    subplot(2,1,1)
    plot(mData.rightAnkle.times, mData.rightAnkle.fGyro); grid on; xlim([20,40]);
    legend('X','Y','Z'); title(['Test ',num2str(testNum),' Right Ankle Filtered Gyro After Rotating']);
    xlabel('Time (sec)'); ylabel('Angular Velocities (deg/sec)'); 
    subplot(2,1,2)
    plot(mData.leftAnkle.times, mData.leftAnkle.fGyro); grid on; xlim([20,40]);
    legend('X','Y','Z'); title(['Test ',num2str(testNum),' Left Ankle Filtered Gyro After Rotating']);
    xlabel('Time (sec)'); ylabel('Angular Velocities (deg/sec)'); 
    set(fig0, 'PaperOrientation', 'landscape', 'PaperUnits', 'normalized', 'PaperPosition',[0,0,1,1]);
    saveas(fig0, ['Shimmer Vicon Results/subj',num2str(testNum),'RotationGyro.pdf']);
end

%% Finds Angles using complimentary filter 
waist.preAngles = compFusionAngles(compHighPass, mData.waist.fPreAcc, mData.waist.fPreGyro, mData.waist.fPreMag, mData.waist.sampleFreq);
waist.angles = compFusionAngles(compHighPass, mData.waist.fAcc, mData.waist.fGyro, mData.waist.fMag, mData.waist.sampleFreq);

rAnkle.preAngles = compFusionAngles(compHighPass, mData.rightAnkle.fPreAcc, mData.rightAnkle.fPreGyro, mData.rightAnkle.fPreMag, mData.rightAnkle.sampleFreq);
rAnkle.angles = compFusionAngles(compHighPass, mData.rightAnkle.fAcc, mData.rightAnkle.fGyro, mData.rightAnkle.fMag, mData.rightAnkle.sampleFreq);

lAnkle.preAngles = compFusionAngles(compHighPass, mData.leftAnkle.fPreAcc, mData.leftAnkle.fPreGyro, mData.leftAnkle.fPreMag, mData.leftAnkle.sampleFreq);
lAnkle.angles = compFusionAngles(compHighPass, mData.leftAnkle.fAcc, mData.leftAnkle.fGyro, mData.leftAnkle.fMag, mData.leftAnkle.sampleFreq);

%% Calculate Global Accelerations and removes baseline
% +x: direction of walking 
% +y: medial to Left Side
% +z: normal to ground (vertical)
if method ~= 5 || method ~= 6
    waist.preGlobalAcc  = imuGlobalAcc(waist.preAngles, mData.waist.fPreAcc);
    waist.globalAcc     = imuGlobalAcc(waist.angles, mData.waist.fAcc);
    waist.globalAcc     = waist.globalAcc - mean(waist.preGlobalAcc);
    rAnkle.preGlobalAcc = imuGlobalAcc(rAnkle.preAngles, mData.rightAnkle.fPreAcc);
    rAnkle.globalAcc    = imuGlobalAcc(rAnkle.angles, mData.rightAnkle.fAcc);
    rAnkle.globalAcc    = rAnkle.globalAcc - mean(rAnkle.preGlobalAcc);
    lAnkle.preGlobalAcc = imuGlobalAcc(lAnkle.preAngles, mData.leftAnkle.fPreAcc);
    lAnkle.globalAcc    = imuGlobalAcc(lAnkle.angles, mData.leftAnkle.fAcc);
    lAnkle.globalAcc    = lAnkle.globalAcc - mean(lAnkle.preGlobalAcc);
end

%% Method #1 
% Finds Heel Strikes and Resamples data
if method == 1
    % Segmentation of Data
    rAnkle.steps = findHSIndexes(mData.rightAnkle.fGyro, mData.rightAnkle.times, sigma);
    lAnkle.steps = findHSIndexes(mData.leftAnkle.fGyro, mData.leftAnkle.times, sigma);
    
    % Data Resampling
    rAnkle.Step = imuDataResampling(rAnkle.steps.indexes, mData.rightAnkle.times, rAnkle.globalAcc, mData.rightAnkle.fGyro, rAnkle.angles);
    lAnkle.Step = imuDataResampling(lAnkle.steps.indexes, mData.leftAnkle.times, lAnkle.globalAcc, mData.leftAnkle.fGyro, lAnkle.angles);
    
    % Walking Speed
    [waist, rAnkle, lAnkle] = gaitSpeed(mData.waist, waist, rAnkle, lAnkle);
    
    % Saves data
    writematrix([waist.avgSpeed, rAnkle.avgSpeed, lAnkle.avgSpeed],'Shimmer Vicon Results/shimmerViconVel.xlsx',...
                'Sheet',['Sheet',num2str(floor(testNum/100))],'Range',['E',num2str(testNum-100*floor(testNum/100)+2)])
    
    % Plot the calculated data
    plotStepIMUdata(testNum, 'RightAnkle', rAnkle.steps, 'HS', ...
        rAnkle.globalAcc, mData.rightAnkle.fGyro, rAnkle.angles, mData.rightAnkle.times, ...
        rAnkle.Step.meanAcc, rAnkle.Step.meanGyro, rAnkle.Step.meanAngle, ...
        rAnkle.Step.stdAcc, rAnkle.Step.stdGyro, rAnkle.Step.stdAngle, rAnkle.Step.avgTimes)
    plotStepIMUdata(testNum, 'LeftAnkle', lAnkle.steps, 'HS', ...
        lAnkle.globalAcc, mData.leftAnkle.fGyro, lAnkle.angles, mData.leftAnkle.times, ...
        lAnkle.Step.meanAcc, lAnkle.Step.meanGyro, lAnkle.Step.meanAngle, ...
        lAnkle.Step.stdAcc, lAnkle.Step.stdGyro, lAnkle.Step.stdAngle, lAnkle.Step.avgTimes)
end

%% Method #2
% Finds Heel strikes and DOES NOT resample data. 
if method == 2
    % Segmentation of Data
    rAnkle.steps = findHSIndexes(mData.rightAnkle.fGyro, mData.rightAnkle.times, sigma);
    lAnkle.steps = findHSIndexes(mData.leftAnkle.fGyro, mData.leftAnkle.times, sigma);
    
    % Speeds per step
    rAnkle.steps = cycleIMUdata(rAnkle.steps, rAnkle.globalAcc, mData.rightAnkle.times);
    lAnkle.steps = cycleIMUdata(lAnkle.steps,lAnkle.globalAcc, mData.leftAnkle.times);
    if rAnkle.steps.indexes{1}(1) < lAnkle.steps.indexes{1}(1)
        waist.steps = cycleIMUdata(rAnkle.steps, waist.globalAcc, mData.waist.times);
    else
        waist.steps = cycleIMUdata(lAnkle.steps, waist.globalAcc, mData.waist.times);
    end
   
    % Saves data
    writematrix([mean(waist.steps.Speed), mean(rAnkle.steps.Speed), mean(lAnkle.steps.Speed)],'Shimmer Vicon Results/shimmerViconVel.xlsx',...
                'Sheet',['Sheet',num2str(floor(testNum/100))],'Range',['H',num2str(testNum-100*floor(testNum/100)+2)])
end

%% Method #3
% Finds Mid Stance and Resamaples data. 
if method == 3
    % Segmentation of Data
    rAnkle.steps = findMSIndexes(mData.rightAnkle.fGyro, rAnkle.angles, mData.rightAnkle.times, sigma);
    lAnkle.steps = findMSIndexes(mData.leftAnkle.fGyro, lAnkle.angles, mData.leftAnkle.times, sigma);
    
    % Data Resampling
    rAnkle.Step = imuDataResampling(rAnkle.steps.indexes, mData.rightAnkle.times, rAnkle.globalAcc, mData.rightAnkle.fGyro, rAnkle.angles);
    lAnkle.Step = imuDataResampling(lAnkle.steps.indexes, mData.leftAnkle.times, lAnkle.globalAcc, mData.leftAnkle.fGyro, lAnkle.angles);
    
    % Walking Speed
    [waist, rAnkle, lAnkle] = gaitSpeed(mData.waist, waist, rAnkle, lAnkle);
    
    % Saves data
    writematrix([waist.avgSpeed, rAnkle.avgSpeed, lAnkle.avgSpeed],'Shimmer Vicon Results/shimmerViconVel.xlsx',...
                'Sheet',['Sheet',num2str(floor(testNum/100))],'Range',['K',num2str(testNum-100*floor(testNum/100)+2)])

    % Plot the calculated data
    plotStepIMUdata(testNum, 'RightAnkle', rAnkle.steps, 'MS', ...
        rAnkle.globalAcc, mData.rightAnkle.fGyro, rAnkle.angles, mData.rightAnkle.times, ...
        rAnkle.Step.meanAcc, rAnkle.Step.meanGyro, rAnkle.Step.meanAngle, ...
        rAnkle.Step.stdAcc, rAnkle.Step.stdGyro, rAnkle.Step.stdAngle, rAnkle.Step.avgTimes)
    plotStepIMUdata(testNum, 'LeftAnkle', lAnkle.steps, 'MS', ...
        lAnkle.globalAcc, mData.leftAnkle.fGyro, lAnkle.angles, mData.leftAnkle.times, ...
        lAnkle.Step.meanAcc, lAnkle.Step.meanGyro, lAnkle.Step.meanAngle, ...
        lAnkle.Step.stdAcc, lAnkle.Step.stdGyro, lAnkle.Step.stdAngle, lAnkle.Step.avgTimes)
end

%% Method #4
% Finds Mid Stance and DOES NOT resample data.
if method == 4
    % Segmentation of Data
    rAnkle.steps = findMSIndexes(mData.rightAnkle.fGyro, rAnkle.angles, mData.rightAnkle.times, sigma);
    lAnkle.steps = findMSIndexes(mData.leftAnkle.fGyro, lAnkle.angles, mData.leftAnkle.times, sigma);
    
    % Speeds per step
    rAnkle.steps = cycleIMUdata(rAnkle.steps, rAnkle.globalAcc, mData.rightAnkle.times);
    lAnkle.steps = cycleIMUdata(lAnkle.steps,lAnkle.globalAcc, mData.leftAnkle.times);
    if rAnkle.steps.indexes{1}(1) < lAnkle.steps.indexes{1}(1)
        waist.steps = cycleIMUdata(rAnkle.steps, waist.globalAcc, mData.waist.times);
    else
        waist.steps = cycleIMUdata(lAnkle.steps, waist.globalAcc, mData.waist.times);
    end
    
    % Saves data
    writematrix([mean(waist.steps.Speed), mean(rAnkle.steps.Speed), mean(lAnkle.steps.Speed)],'Shimmer Vicon Results/shimmerViconVel.xlsx',...
                'Sheet',['Sheet',num2str(floor(testNum/100))],'Range',['N',num2str(testNum-100*floor(testNum/100)+2)])
end

%% Method #5
% Finds Global accelerations using 2 angles, then finds mid stance 
if method == 5
    % Calculate Global Accelerations and removes baseline
    % +x: direction of walking
    % +y: medial to Left Side
    % +z: normal to ground (vertical)
    waist.preGlobalAcc  = imuGlobalAcc2Angles(waist.preAngles, mData.waist.fPreAcc);
    waist.globalAcc     = imuGlobalAcc2Angles(waist.angles, mData.waist.fAcc);
    waist.globalAcc     = waist.globalAcc - mean(waist.preGlobalAcc);
    rAnkle.preGlobalAcc = imuGlobalAcc2Angles(rAnkle.preAngles, mData.rightAnkle.fPreAcc);
    rAnkle.globalAcc    = imuGlobalAcc2Angles(rAnkle.angles, mData.rightAnkle.fAcc);
    rAnkle.globalAcc    = rAnkle.globalAcc - mean(rAnkle.preGlobalAcc);
    lAnkle.preGlobalAcc = imuGlobalAcc2Angles(lAnkle.preAngles, mData.leftAnkle.fPreAcc);
    lAnkle.globalAcc    = imuGlobalAcc2Angles(lAnkle.angles, mData.leftAnkle.fAcc);
    lAnkle.globalAcc    = lAnkle.globalAcc - mean(lAnkle.preGlobalAcc);
    
    % Segmentation of Data
    rAnkle.steps = findMSIndexes(mData.rightAnkle.fGyro, rAnkle.angles, mData.rightAnkle.times, sigma);
    lAnkle.steps = findMSIndexes(mData.leftAnkle.fGyro, lAnkle.angles, mData.leftAnkle.times, sigma);
    
    % Speeds per step
    rAnkle.steps = cycleIMUdata(rAnkle.steps, rAnkle.globalAcc, mData.rightAnkle.times);
    lAnkle.steps = cycleIMUdata(lAnkle.steps,lAnkle.globalAcc, mData.leftAnkle.times);
    if rAnkle.steps.indexes{1}(1) < lAnkle.steps.indexes{1}(1)
        waist.steps = cycleIMUdata(rAnkle.steps, waist.globalAcc, mData.waist.times);
    else
        waist.steps = cycleIMUdata(lAnkle.steps, waist.globalAcc, mData.waist.times);
    end
    
    % Saves data
    writematrix([mean(waist.steps.Speed), mean(rAnkle.steps.Speed), mean(lAnkle.steps.Speed)],'Shimmer Vicon Results/shimmerViconVel.xlsx',...
                'Sheet',['Sheet',num2str(floor(testNum/100))],'Range',['Q',num2str(testNum-100*floor(testNum/100)+2)])
end

%% Method #6
% Finds global accelerations using 3 angles, then finds mid stance
if method == 6
    % Calculate Global Accelerations and removes baseline
    % +x: direction of walking
    % +y: medial to Left Side
    % +z: normal to ground (vertical)
    waist.preGlobalAcc  = imuGlobalAcc3Angles(waist.preAngles, mData.waist.fPreAcc);
    waist.globalAcc     = imuGlobalAcc3Angles(waist.angles, mData.waist.fAcc);
    waist.globalAcc     = waist.globalAcc - mean(waist.preGlobalAcc);
    rAnkle.preGlobalAcc = imuGlobalAcc3Angles(rAnkle.preAngles, mData.rightAnkle.fPreAcc);
    rAnkle.globalAcc    = imuGlobalAcc3Angles(rAnkle.angles, mData.rightAnkle.fAcc);
    rAnkle.globalAcc    = rAnkle.globalAcc - mean(rAnkle.preGlobalAcc);
    lAnkle.preGlobalAcc = imuGlobalAcc3Angles(lAnkle.preAngles, mData.leftAnkle.fPreAcc);
    lAnkle.globalAcc    = imuGlobalAcc3Angles(lAnkle.angles, mData.leftAnkle.fAcc);
    lAnkle.globalAcc    = lAnkle.globalAcc - mean(lAnkle.preGlobalAcc);
    
    % Segmentation of Data
    rAnkle.steps = findMSIndexes(mData.rightAnkle.fGyro, rAnkle.angles, mData.rightAnkle.times, sigma);
    lAnkle.steps = findMSIndexes(mData.leftAnkle.fGyro, lAnkle.angles, mData.leftAnkle.times, sigma);
    
    % Speeds per step
    rAnkle.steps = cycleIMUdata(rAnkle.steps, rAnkle.globalAcc, mData.rightAnkle.times);
    lAnkle.steps = cycleIMUdata(lAnkle.steps,lAnkle.globalAcc, mData.leftAnkle.times);
    if rAnkle.steps.indexes{1}(1) < lAnkle.steps.indexes{1}(1)
        waist.steps = cycleIMUdata(rAnkle.steps, waist.globalAcc, mData.waist.times);
    else
        waist.steps = cycleIMUdata(lAnkle.steps, waist.globalAcc, mData.waist.times);
    end
    
    % Saves data
    writematrix([mean(waist.steps.Speed), mean(rAnkle.steps.Speed), mean(lAnkle.steps.Speed)],'Shimmer Vicon Results/shimmerViconVel.xlsx',...
                'Sheet',['Sheet',num2str(floor(testNum/100))],'Range',['T',num2str(testNum-100*floor(testNum/100)+2)])
end

%% Method #7
% Finds mid stance based on mid swings first
if method == 7
    rAnkle.gyro = mData.rightAnkle.fGyro;
    lAnkle.gyro = mData.leftAnkle.fGyro;
    waist.gyro = mData.waist.fGyro;
    rAnkle.times = mData.rightAnkle.times;
    lAnkle.times = mData.leftAnkle.times;
    waist.times = mData.waist.times;
    removeY = false;
    % Segmentation of Data
    rAnkle.steps = findIndexes(mData.rightAnkle.fGyro, mData.rightAnkle.times, sigma, removeY);
    lAnkle.steps = findIndexes(mData.leftAnkle.fGyro, mData.leftAnkle.times, sigma,removeY);
    
    % Speeds per step
    rAnkle.steps = cycleIMUdata(rAnkle.steps, rAnkle.globalAcc, mData.rightAnkle.times);
    lAnkle.steps = cycleIMUdata(lAnkle.steps,lAnkle.globalAcc, mData.leftAnkle.times);
    if rAnkle.steps.indexes{1}(1) < lAnkle.steps.indexes{1}(1)
        waist.steps = cycleIMUdata(rAnkle.steps, waist.globalAcc, mData.waist.times);
    else
        waist.steps = cycleIMUdata(lAnkle.steps, waist.globalAcc, mData.waist.times);
    end
    
    % Finds Cycle profile
    cycle = gaitCycle(rAnkle, lAnkle, waist);
    
    % Plot calulated data
    plotMoreIMUdata(testNum, 'MS', rAnkle, lAnkle, waist, cycle,removeY)
    
    % Saves data
    writematrix([mean([rAnkle.steps.Distance,lAnkle.steps.Distance]),...
                 mean([rAnkle.steps.Duration,lAnkle.steps.Duration]),100,...
                 mean(waist.steps.Speed),mean(rAnkle.steps.Speed),mean(lAnkle.steps.Speed)],...
                'Shimmer Vicon Results/shimmerViconVel.xlsx',...
                'Sheet',['Sheet',num2str(floor(testNum/100))],...
                'Range',['W',num2str(testNum-100*floor(testNum/100)+2)])
end

%% Method #8
% Finds mid stance based on mid swings first and then removes y outliers
if method == 8
    rAnkle.gyro = mData.rightAnkle.fGyro;
    lAnkle.gyro = mData.leftAnkle.fGyro;
    waist.gyro = mData.waist.fGyro;
    rAnkle.times = mData.rightAnkle.times;
    lAnkle.times = mData.leftAnkle.times;
    waist.times = mData.waist.times;
    removeY = true;
    % Segmentation of Data
    rAnkle.steps = findIndexes(mData.rightAnkle.fGyro, mData.rightAnkle.times, sigma, removeY);
    lAnkle.steps = findIndexes(mData.leftAnkle.fGyro, mData.leftAnkle.times, sigma,removeY);
    
    % Speeds per step
    rAnkle.steps = cycleIMUdata(rAnkle.steps, rAnkle.globalAcc, mData.rightAnkle.times);
    lAnkle.steps = cycleIMUdata(lAnkle.steps,lAnkle.globalAcc, mData.leftAnkle.times);
    if rAnkle.steps.indexes{1}(1) < lAnkle.steps.indexes{1}(1)
        waist.steps = cycleIMUdata(rAnkle.steps, waist.globalAcc, mData.waist.times);
    else
        waist.steps = cycleIMUdata(lAnkle.steps, waist.globalAcc, mData.waist.times);
    end
    
    % Finds Cycle profile
    cycle = gaitCycle(rAnkle, lAnkle, waist);
    
    % Plot calulated data
    plotMoreIMUdata(testNum, 'MS', rAnkle, lAnkle, waist, cycle,removeY)
    
    % Saves data
    writematrix([mean([rAnkle.steps.Distance,lAnkle.steps.Distance]),...
                 mean([rAnkle.steps.Duration,lAnkle.steps.Duration]),100,...
                 mean(waist.steps.Speed),mean(rAnkle.steps.Speed),mean(lAnkle.steps.Speed)],...
                'Shimmer Vicon Results/shimmerViconVel.xlsx',...
                'Sheet',['Sheet',num2str(floor(testNum/100))],...
                'Range',['AC',num2str(testNum-100*floor(testNum/100)+2)])
end

%% Methods #9
if method == 9
    rAnkle.gyro = mData.rightAnkle.fGyro;
    lAnkle.gyro = mData.leftAnkle.fGyro;
    waist.gyro = mData.waist.fGyro;
    rAnkle.times = mData.rightAnkle.times;
    lAnkle.times = mData.leftAnkle.times;
    waist.times = mData.waist.times;
    removeY = false;
    
    % Segmentation of Data
    rAnkle.steps = findIndexes(mData.rightAnkle.fGyro, mData.rightAnkle.times, sigma, removeY);
    lAnkle.steps = findIndexes(mData.leftAnkle.fGyro, mData.leftAnkle.times, sigma,removeY);
    
    % Finds Speeds
    rAnkle.steps.speeds = imuSpeed(rAnkle);
    lAnkle.steps.speeds = imuSpeed(lAnkle);
    if rAnkle.steps.indexes{1}(1) < lAnkle.steps.indexes{1}(1)
        waist.steps.speeds  = imuSpeed(waist, rAnkle.steps.indexes);
    else
        waist.steps.speeds  = imuSpeed(waist, lAnkle.steps.indexes);
    end
    % Saves data
    writematrix([mean(waist.steps.speeds.meanSpeeds),mean(rAnkle.steps.speeds.meanSpeeds),mean(lAnkle.steps.speeds.meanSpeeds)],...
                'Shimmer Vicon Results/shimmerViconVel.xlsx',...
                'Sheet',['Sheet',num2str(floor(testNum/100))],...
                'Range',['AI',num2str(testNum-100*floor(testNum/100)+2)])
    
end

%% Method #10
if method == 10
    rAnkle.gyro = mData.rightAnkle.fGyro;
    lAnkle.gyro = mData.leftAnkle.fGyro;
    waist.gyro = mData.waist.fGyro;
    rAnkle.times = mData.rightAnkle.times;
    lAnkle.times = mData.leftAnkle.times;
    waist.times = mData.waist.times;
    removeY = false;
    
    % Segmentation of Data
    rAnkle.steps = findIndexes(mData.rightAnkle.fGyro, mData.rightAnkle.times, sigma, removeY);
    lAnkle.steps = findIndexes(mData.leftAnkle.fGyro, mData.leftAnkle.times, sigma,removeY);
    
    % Inverted Pendulum
    rWaist = invertedPendulum(waist, rAnkle.steps, height);
    lWaist = invertedPendulum(waist, lAnkle.steps, height);
    
    % Save data
    writematrix([mean(rWaist.steps.speed), mean(lWaist.steps.speed), mean(rWaist.steps.stepLen), mean(lWaist.steps.stepLen)],...
                'Shimmer Vicon Results/shimmerViconVel.xlsx',...
                'Sheet',['Sheet',num2str(floor(testNum/100))],...
                'Range',['AL',num2str(testNum-100*floor(testNum/100)+2)])
    
end

%% Troubleshooting Method
% if method == 100
%     rAnkle.steps = findMSIndexes(mData.rightAnkle.fGyro, rAnkle.angles, mData.rightAnkle.times, sigma);
%     rAnkle.angAcc = gradient(mData.rightAnkle.fGyro, mean(diff(mData.rightAnkle.times)));
%     
%     figure
%     subplot(3,1,1)
%     plot(mData.rightAnkle.times,rAnkle.angles(:,2));
%     grid on; hold on;
%     plot(mData.rightAnkle.times(rAnkle.steps.MSindexes),rAnkle.angles(rAnkle.steps.MSindexes,2),'rs')
%     xlim([0,20]);
%     subplot(3,1,2)
%     plot(mData.rightAnkle.times,mData.rightAnkle.fGyro(:,2));
%     grid on; hold on;
%     plot(mData.rightAnkle.times(rAnkle.steps.MSindexes),mData.rightAnkle.fGyro(rAnkle.steps.MSindexes,2),'rs')
%     plot(mData.rightAnkle.times,abs(mData.rightAnkle.fGyro(:,2)),'r');
%     xlim([0,20]);
%     subplot(3,1,3)
%     plot(mData.rightAnkle.times,rAnkle.angAcc(:,2));
%     grid on; hold on;
%     plot(mData.rightAnkle.times(rAnkle.steps.MSindexes),rAnkle.angAcc(rAnkle.steps.MSindexes,2),'rs')
%     xlim([0,20]);
%     
% end

end