%% IMU Feature extraction
% Script by Erick Nunez

%% clean up work space
clear; clc; close all;
%% Subject and session variables
% subjects 1, 2, 4, 5, 6, 7, 12, 13 have good gait data. 
subj.Num = 13;
sensorLoc = {'leftAnkle','rightAnkle','waist'};
% sensorLoc = {'leftAnkle','leftWrist','rightAnkle','rightWrist','waist'};
session.freqC = 2; session.filtOrder = 2; gravity = 9.81;
session.compHigh = 0.98; session.compLow = 1-session.compHigh;

%% Load Mat-file
load(['../../../Google Drive/School/PhD - Biomedical/Rutgers - Rotation/Gait Data/subjectMatLabData-Rot/subject',num2str(subj.Num),'RawDataLogs-Rot.mat']);

%% Session times
for i=1:length(sensorLoc)
    session.(sensorLoc{i}).startTime = raw.(sensorLoc{i})(1,1);
    session.(sensorLoc{i}).unixTimes = raw.(sensorLoc{i})(:,1);
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
%fuseIMU = imufilter('SampleRate',session.freqS,'DecimationFactor',1);
[session.filtB, session.filtA] = butter(session.filtOrder,session.freqC/(session.freqS/2));

%% Subject times and indexes
subj.preUnixTime = logFile(subj.Num,21);
subj.unixStartTime = logFile(subj.Num,23);
subj.unixEndTime = logFile(subj.Num,25);
if subj.preUnixTime > 1552201200
    subj.preUnixTime = (subj.preUnixTime + 14400)*1000;
    subj.unixStartTime = (subj.unixStartTime + 14400)*1000;
    subj.unixEndTime = (subj.unixEndTime + 14400)*1000;
else
    subj.preUnixTime = (subj.preUnixTime + 18000)*1000;
    subj.unixStartTime = (subj.unixStartTime + 18000)*1000;
    subj.unixEndTime = (subj.unixEndTime + 18000)*1000;
end

for i=1:length(sensorLoc)
    subj.(sensorLoc{i}).preIndex = 1;
    while session.(sensorLoc{i}).unixTimes(subj.(sensorLoc{i}).preIndex) < subj.preUnixTime
        subj.(sensorLoc{i}).preIndex = subj.(sensorLoc{i}).preIndex + 1;
    end
    subj.(sensorLoc{i}).preIndex = subj.(sensorLoc{i}).preIndex - 1;
    subj.(sensorLoc{i}).startIndex = subj.(sensorLoc{i}).preIndex;
    while session.(sensorLoc{i}).unixTimes(subj.(sensorLoc{i}).startIndex) < subj.unixStartTime
        subj.(sensorLoc{i}).startIndex = subj.(sensorLoc{i}).startIndex + 1;
    end
    subj.(sensorLoc{i}).startIndex = subj.(sensorLoc{i}).startIndex - 1;
    subj.(sensorLoc{i}).endIndex = subj.(sensorLoc{i}).startIndex;
    while session.(sensorLoc{i}).unixTimes(subj.(sensorLoc{i}).endIndex) < subj.unixEndTime
        subj.(sensorLoc{i}).endIndex = subj.(sensorLoc{i}).endIndex + 1;
    end
    subj.(sensorLoc{i}).preTimes = (raw.(sensorLoc{i})(subj.(sensorLoc{i}).preIndex:subj.(sensorLoc{i}).startIndex,1) - raw.(sensorLoc{i})(subj.(sensorLoc{i}).preIndex,1))/1000;
    subj.(sensorLoc{i}).times = (raw.(sensorLoc{i})(subj.(sensorLoc{i}).startIndex:subj.(sensorLoc{i}).endIndex,1) - raw.(sensorLoc{i})(subj.(sensorLoc{i}).startIndex,1))/1000;
end

%% Pulls Subject IMU data
for i = 1:length(sensorLoc)
    % Baseline IMU Data for subject
    subj.(sensorLoc{i}).preAcc = raw.(sensorLoc{i})(subj.(sensorLoc{i}).preIndex:subj.(sensorLoc{i}).startIndex,2:4);
    subj.(sensorLoc{i}).preGyro = raw.(sensorLoc{i})(subj.(sensorLoc{i}).preIndex:subj.(sensorLoc{i}).startIndex,6:8);
    subj.(sensorLoc{i}).preMag = raw.(sensorLoc{i})(subj.(sensorLoc{i}).preIndex:subj.(sensorLoc{i}).startIndex,9:11);
    % IMU Data for subject trial (m/s^2)
    subj.(sensorLoc{i}).acc = raw.(sensorLoc{i})(subj.(sensorLoc{i}).startIndex:subj.(sensorLoc{i}).endIndex,2:4);
    subj.(sensorLoc{i}).gyro = raw.(sensorLoc{i})(subj.(sensorLoc{i}).startIndex:subj.(sensorLoc{i}).endIndex,6:8);
    subj.(sensorLoc{i}).mag = raw.(sensorLoc{i})(subj.(sensorLoc{i}).startIndex:subj.(sensorLoc{i}).endIndex,9:11);
end

%% Finds orientation of data
% for i = 1:length(sensorLoc)
%    subj.(sensorLoc{i}).preQuaternion = fuseIMU(subj.(sensorLoc{i}).preAcc,subj.(sensorLoc{i}).preGyro);
%    subj.(sensorLoc{i}).preAngle = eulerd(subj.(sensorLoc{i}).preQuaternion,'XYZ','frame');
%    subj.(sensorLoc{i}).quaternion = fuseIMU(subj.(sensorLoc{i}).acc,subj.(sensorLoc{i}).gyro);
%    subj.(sensorLoc{i}).angle = eulerd(subj.(sensorLoc{i}).quaternion,'XYZ','frame');
% end

%% Filters subject IMU data
for i = 1:length(sensorLoc)
    for j = 1:3
        % Filter Baseline IMU Data for subject
        subj.(sensorLoc{i}).fPreAcc(:,j) = filtfilt(session.filtB,session.filtA,subj.(sensorLoc{i}).preAcc(:,j));
        subj.(sensorLoc{i}).fPreGyro(:,j) = filtfilt(session.filtB,session.filtA,subj.(sensorLoc{i}).preGyro(:,j));
        subj.(sensorLoc{i}).fPreMag(:,j) = filtfilt(session.filtB,session.filtA,subj.(sensorLoc{i}).preMag(:,j));
        % Filter trial IMU Data for subject
        subj.(sensorLoc{i}).fAcc(:,j) = filtfilt(session.filtB,session.filtA,subj.(sensorLoc{i}).acc(:,j));
        subj.(sensorLoc{i}).fGyro(:,j) = filtfilt(session.filtB,session.filtA,subj.(sensorLoc{i}).gyro(:,j));
        subj.(sensorLoc{i}).fMag(:,j) = filtfilt(session.filtB,session.filtA,subj.(sensorLoc{i}).mag(:,j));
    end
end

%% Finding angles and fuses data from IMUs
% gryoAng are angles integrated from the gyrocope readings
% accAng are angles derived from the acceleration vectors
% compAng(i) = highPass*(compAng(i-1)+gyro*dt)+lowPass*accAng
for i = 1:length(sensorLoc)
    subj.(sensorLoc{i}).preGyroAng(1,:) = subj.(sensorLoc{i}).fPreGyro(1,:) .* subj.(sensorLoc{i}).preTimes(1);
    for j = 1:3
        for k = 2:length(subj.(sensorLoc{i}).preTimes)
            subj.(sensorLoc{i}).preGyroAng(k,j) = ...
                subj.(sensorLoc{i}).preGyroAng(k-1,j) + ...
                (subj.(sensorLoc{i}).fPreGyro(k,j)+subj.(sensorLoc{i}).fPreGyro(k-1,j))/2 * (subj.(sensorLoc{i}).preTimes(k)-subj.(sensorLoc{i}).preTimes(k-1));
        end
    end
end

for i = 1:length(sensorLoc)
    for j = 1:length(subj.(sensorLoc{i}).preTimes)
        if strcmp(sensorLoc{i},'waist')
            subj.(sensorLoc{i}).preAccAng(j,1) = atan2d(-subj.(sensorLoc{i}).fPreAcc(j,3),sqrt((subj.(sensorLoc{i}).fPreAcc(j,1)^2)+(subj.(sensorLoc{i}).fPreAcc(j,2)^2)));
            subj.(sensorLoc{i}).preAccAng(j,3) = atan2d(subj.(sensorLoc{i}).fPreAcc(j,1),subj.(sensorLoc{i}).fPreAcc(j,2));
            Mx = subj.(sensorLoc{i}).fPreMag(j,3) * sind(subj.(sensorLoc{i}).preAccAng(j,3)) * sind(subj.(sensorLoc{i}).preAccAng(j,1)) + ...
                 subj.(sensorLoc{i}).fPreMag(j,1) * cosd(subj.(sensorLoc{i}).preAccAng(j,3)) - ...
                 subj.(sensorLoc{i}).fPreMag(j,2) * sind(subj.(sensorLoc{i}).preAccAng(j,3)) * cosd(subj.(sensorLoc{i}).preAccAng(j,1));
            Mz = subj.(sensorLoc{i}).fPreMag(j,3) * cosd(subj.(sensorLoc{i}).preAccAng(j,1)) + ...
                 subj.(sensorLoc{i}).fPreMag(j,2) * sind(subj.(sensorLoc{i}).preAccAng(j,1));
            subj.(sensorLoc{i}).preAccAng(j,2) = atan2d(Mx,Mz);
        else
            subj.(sensorLoc{i}).preAccAng(j,1) = atan2d(subj.(sensorLoc{i}).fPreAcc(j,3),subj.(sensorLoc{i}).fPreAcc(j,2));
            subj.(sensorLoc{i}).preAccAng(j,3) = atan2d(-subj.(sensorLoc{i}).fPreAcc(j,1),sqrt((subj.(sensorLoc{i}).fPreAcc(j,2)^2)+(subj.(sensorLoc{i}).fPreAcc(j,3)^2)));
            Mx = subj.(sensorLoc{i}).fPreMag(j,1) * cosd(subj.(sensorLoc{i}).preAccAng(j,3)) + ...
                 subj.(sensorLoc{i}).fPreMag(j,2) * sind(subj.(sensorLoc{i}).preAccAng(j,3));
            Mz = subj.(sensorLoc{i}).fPreMag(j,1) * sind(subj.(sensorLoc{i}).preAccAng(j,1)) * sind(subj.(sensorLoc{i}).preAccAng(j,3)) + ...
                 subj.(sensorLoc{i}).fPreMag(j,3) * cosd(subj.(sensorLoc{i}).preAccAng(j,1)) - ...
                 subj.(sensorLoc{i}).fPreMag(j,2) * sind(subj.(sensorLoc{i}).preAccAng(j,1)) * cosd(subj.(sensorLoc{i}).preAccAng(j,3));
            subj.(sensorLoc{i}).preAccAng(j,2) = atan2d(Mz,Mx);
        end
    end
end

for i = 1:length(sensorLoc)
    subj.(sensorLoc{i}).preCompAng(1,:) = session.compLow .* subj.(sensorLoc{i}).preAccAng(1,:);
    for j = 1:3
        for k = 2:length(subj.(sensorLoc{i}).preTimes)
            subj.(sensorLoc{i}).preCompAng(k,j) = ...
                session.compHigh * (subj.(sensorLoc{i}).preCompAng(k-1,j) + subj.(sensorLoc{i}).preGyro(k,j) * (subj.(sensorLoc{i}).preTimes(k)-subj.(sensorLoc{i}).preTimes(k-1))) + ...
                session.compLow * subj.(sensorLoc{i}).preAccAng(k,j);
        end
    end
end

for i = 1:length(sensorLoc)
    subj.(sensorLoc{i}).gyroAng(1,:) = subj.(sensorLoc{i}).fGyro(1,:) .* subj.(sensorLoc{i}).times(1);
    for j = 1:3
        for k = 2:length(subj.(sensorLoc{i}).times)
            subj.(sensorLoc{i}).gyroAng(k,j) = ...
                subj.(sensorLoc{i}).gyroAng(k-1,j) + ...
                (subj.(sensorLoc{i}).fGyro(k,j)+subj.(sensorLoc{i}).fGyro(k-1,j))/2 * (subj.(sensorLoc{i}).times(k)-subj.(sensorLoc{i}).times(k-1));
        end
    end
end

for i = 1:length(sensorLoc)
    for j = 1:length(subj.(sensorLoc{i}).times)
        if strcmp(sensorLoc{i},'waist')
            subj.(sensorLoc{i}).accAng(j,1) = atan2d(-subj.(sensorLoc{i}).fAcc(j,3),sqrt((subj.(sensorLoc{i}).fAcc(j,1)^2)+(subj.(sensorLoc{i}).fAcc(j,2)^2)));
            subj.(sensorLoc{i}).accAng(j,3) = atan2d(subj.(sensorLoc{i}).fAcc(j,1),subj.(sensorLoc{i}).fAcc(j,2));
            Mx = subj.(sensorLoc{i}).fMag(j,3) * sind(subj.(sensorLoc{i}).accAng(j,3)) * sind(subj.(sensorLoc{i}).accAng(j,1)) + ...
                 subj.(sensorLoc{i}).fMag(j,1) * cosd(subj.(sensorLoc{i}).accAng(j,3)) - ...
                 subj.(sensorLoc{i}).fMag(j,2) * sind(subj.(sensorLoc{i}).accAng(j,3)) * cosd(subj.(sensorLoc{i}).accAng(j,1));
            Mz = subj.(sensorLoc{i}).fMag(j,3) * cosd(subj.(sensorLoc{i}).accAng(j,1)) + ...
                 subj.(sensorLoc{i}).fMag(j,2) * sind(subj.(sensorLoc{i}).accAng(j,1));
            subj.(sensorLoc{i}).accAng(j,2) = atan2d(Mx,Mz);
        else
            subj.(sensorLoc{i}).accAng(j,1) = atan2d(subj.(sensorLoc{i}).fAcc(j,3),subj.(sensorLoc{i}).fAcc(j,2));
            subj.(sensorLoc{i}).accAng(j,3) = atan2d(-subj.(sensorLoc{i}).fAcc(j,1),sqrt((subj.(sensorLoc{i}).fAcc(j,2)^2)+(subj.(sensorLoc{i}).fAcc(j,3)^2)));
            Mx = subj.(sensorLoc{i}).fMag(j,1) * cosd(subj.(sensorLoc{i}).accAng(j,3)) + ...
                 subj.(sensorLoc{i}).fMag(j,2) * sind(subj.(sensorLoc{i}).accAng(j,3));
            MZ = subj.(sensorLoc{i}).fMag(j,1) * sind(subj.(sensorLoc{i}).accAng(j,1)) * sind(subj.(sensorLoc{i}).accAng(j,3)) + ...
                 subj.(sensorLoc{i}).fMag(j,3) * cosd(subj.(sensorLoc{i}).accAng(j,1)) - ...
                 subj.(sensorLoc{i}).fMag(j,2) * sind(subj.(sensorLoc{i}).accAng(j,1)) * cosd(subj.(sensorLoc{i}).accAng(j,3));
            subj.(sensorLoc{i}).accAng(j,2) = atan2d(Mz,Mx);
        end
    end
end

for i = 1:length(sensorLoc)
    subj.(sensorLoc{i}).compAng(1,:) = session.compLow .* subj.(sensorLoc{i}).accAng(1,:);
    for j = 1:3
        for k = 2:length(subj.(sensorLoc{i}).times)
            subj.(sensorLoc{i}).compAng(k,j) = ...
                session.compHigh * (subj.(sensorLoc{i}).compAng(k-1,j) + subj.(sensorLoc{i}).gyro(k,j) * (subj.(sensorLoc{i}).times(k)-subj.(sensorLoc{i}).times(k-1))) + ...
                session.compLow * subj.(sensorLoc{i}).accAng(k,j);
        end
    end
end

%% Baseline phase
for i = 1:length(sensorLoc)
   for j = 1:length(subj.(sensorLoc{i}).preTimes)
       if strcmp(sensorLoc{i},'waist')
           subj.(sensorLoc{i}).fPreAcc(j,2) = subj.(sensorLoc{i}).fPreAcc(j,2) - gravity * cosd(subj.(sensorLoc{i}).preCompAng(j,1));
           subj.(sensorLoc{i}).fPreAcc(j,3) = subj.(sensorLoc{i}).fPreAcc(j,3) - gravity * sind(subj.(sensorLoc{i}).preCompAng(j,1));
       else
           subj.(sensorLoc{i}).fPreAcc(j,2) = subj.(sensorLoc{i}).fPreAcc(j,2) - gravity * cosd(subj.(sensorLoc{i}).preCompAng(j,3));
           subj.(sensorLoc{i}).fPreAcc(j,1) = subj.(sensorLoc{i}).fPreAcc(j,1) - gravity * sind(subj.(sensorLoc{i}).preCompAng(j,3));
       end
   end
end

for i = 1:length(sensorLoc)
    for j = 1:3
        subj.(sensorLoc{i}).avgPreAcc(1,j) = sum(subj.(sensorLoc{i}).fPreAcc(:,j))/length(subj.(sensorLoc{i}).preTimes);
    end
end

%% Subject trial phase
for i = 1:length(sensorLoc)
   for j = 1:length(subj.(sensorLoc{i}).times)
       if strcmp(sensorLoc{i},'waist')
           subj.(sensorLoc{i}).fAcc(j,2) = subj.(sensorLoc{i}).fAcc(j,2) - gravity * cos(subj.(sensorLoc{i}).compAng(j,1));
           subj.(sensorLoc{i}).fAcc(j,3) = subj.(sensorLoc{i}).fAcc(j,3) - gravity * sin(subj.(sensorLoc{i}).compAng(j,1));
       else
           subj.(sensorLoc{i}).fAcc(j,2) = subj.(sensorLoc{i}).fAcc(j,2) - gravity * cos(subj.(sensorLoc{i}).compAng(j,3));
           subj.(sensorLoc{i}).fAcc(j,1) = subj.(sensorLoc{i}).fAcc(j,1) - gravity * sin(subj.(sensorLoc{i}).compAng(j,3));
       end
   end
end

for i = 1:length(sensorLoc)
    for j = 1:3
        subj.(sensorLoc{i}).fAcc(:,j) = subj.(sensorLoc{i}).fAcc(:,j) - subj.(sensorLoc{i}).avgPreAcc(1,j);
    end
end

%% Integrate for velocities during trial
for i = 1:length(sensorLoc)
    subj.(sensorLoc{i}).vel(1,:) = subj.(sensorLoc{i}).fAcc(1,:) .* subj.(sensorLoc{i}).times(1);
    for j = 1:3
        for k = 2:length(subj.(sensorLoc{i}).times)
            subj.(sensorLoc{i}).vel(k,j) = ...
                subj.(sensorLoc{i}).vel(k-1,j) + ...
                (subj.(sensorLoc{i}).fAcc(k,j)+subj.(sensorLoc{i}).fAcc(k-1,j))/2 * (subj.(sensorLoc{i}).times(k)-subj.(sensorLoc{i}).times(k-1));
        end
    end
end

for i = 1:length(sensorLoc)
    subj.(sensorLoc{i}).vel(1,:) = subj.(sensorLoc{i}).fAcc(1,:) .* subj.(sensorLoc{i}).times(1);
    for j = 1:3
        subj.(sensorLoc{i}).avgVelNoC(j) = sum(subj.(sensorLoc{i}).vel(:,j))/length(subj.(sensorLoc{i}).vel(:,j));
    end
end

for i = 1:length(sensorLoc)
    subj.(sensorLoc{i}).vel(1,:) = subj.(sensorLoc{i}).fAcc(1,:) .* subj.(sensorLoc{i}).times(1);
    for j = 1:3
        for k = 2:length(subj.(sensorLoc{i}).times)
            subj.(sensorLoc{i}).vel(k,j) = ...
                (subj.(sensorLoc{i}).fAcc(k,j)+subj.(sensorLoc{i}).fAcc(k-1,j))/2 * (subj.(sensorLoc{i}).times(k)-subj.(sensorLoc{i}).times(k-1));
        end
    end
end

for i = 1:length(sensorLoc)
    subj.(sensorLoc{i}).vel(1,:) = subj.(sensorLoc{i}).fAcc(1,:) .* subj.(sensorLoc{i}).times(1);
    for j = 1:3
        subj.(sensorLoc{i}).avgVelC(j) = sum(subj.(sensorLoc{i}).vel(:,j))/length(subj.(sensorLoc{i}).vel(:,j));
    end
end

for i = 1:length(sensorLoc)
    subj.(sensorLoc{i}).vel(1,:) = subj.(sensorLoc{i}).fAcc(1,:) .* subj.(sensorLoc{i}).times(1);
    for j = 1:3
        subj.(sensorLoc{i}).velTrapz(j) = trapz(subj.(sensorLoc{i}).times,subj.(sensorLoc{i}).vel(:,j));
    end
end

writematrix([subj.Num, subj.waist.avgVelNoC, subj.waist.avgVelC, subj.waist.velTrapz],'../../../Google Drive/School/PhD - Biomedical/Rutgers - Rotation/Gait Data/avgVel.xlsx','Range',['A',num2str(subj.Num+2)])

%% Plotting of Velocities and Accelerations
for i = 1:length(sensorLoc)
    figure(i+5)
    subplot(3,2,1)
    plot(subj.(sensorLoc{i}).times, subj.(sensorLoc{i}).accAng); 
    grid on; xlabel('seconds'); ylabel('deg'); legend('X','Y','Z');
    title(['Subject ',num2str(subj.Num),' ',(sensorLoc{i}),' Angles(Accel)'])
    
    subplot(3,2,2)
    plot(subj.(sensorLoc{i}).times, subj.(sensorLoc{i}).vel); 
    grid on; xlabel('seconds'); ylabel('m/s'); legend('X','Y','Z');
    title(['Subject ',num2str(subj.Num),' ',(sensorLoc{i}),' Velocities'])
    
    subplot(3,2,3)
    plot(subj.(sensorLoc{i}).times, subj.(sensorLoc{i}).gyroAng); 
    grid on; xlabel('seconds'); ylabel('deg'); legend('X','Y','Z');
    title(['Subject ',num2str(subj.Num),' ',(sensorLoc{i}),' Angles(Gyro)'])
    
    subplot(3,2,4)
    plot(subj.(sensorLoc{i}).times, subj.(sensorLoc{i}).fAcc); 
    grid on; xlabel('seconds'); ylabel('m/s^2'); legend('X','Y','Z');
    title(['Subject ',num2str(subj.Num),' ',(sensorLoc{i}),' Accelerations'])
    
    subplot(3,2,5)
    plot(subj.(sensorLoc{i}).times, subj.(sensorLoc{i}).compAng); 
    grid on; xlabel('seconds'); ylabel('deg'); legend('X','Y','Z');
    title(['Subject ',num2str(subj.Num),' ',(sensorLoc{i}),' Angles(Comp)'])
    
    subplot(3,2,6)
    plot(subj.(sensorLoc{i}).times, subj.(sensorLoc{i}).fGyro);
    grid on; xlabel('seconds'); ylabel('deg/sec'); legend('X','Y','Z');
    title(['Subject ',num2str(subj.Num),' ',(sensorLoc{i}),' Angular Velocities'])
end
