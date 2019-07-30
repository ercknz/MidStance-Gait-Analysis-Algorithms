%% IMU Feature extraction
% Script by Erick Nunez

%% clean up work space
clear; clc; close all;
%% Subject and session variables
% subjects 1, 2, 4, 5, 6, and 7 have good data. 
subj.Num = 1;
sensorLoc = {'leftAnkle','rightAnkle','waist'};
session.freqC = 5; session.filtOrder = 2; gravity = 9.81;

%% Import Data
[logFile, raw] = imuDataImport(subj);

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
[session.filtB, session.filtA] = butter(session.filtOrder,session.freqC/(session.freqS/2));

%% Subject times and indexes
subj.baseUnixTime = (logFile(subj.Num,21)+18000)*1000;
subj.unixStartTime = (logFile(subj.Num,23)+18000)*1000;
subj.unixEndTime = (logFile(subj.Num,25)+18000)*1000;

for i=1:length(sensorLoc)
    subj.(sensorLoc{i}).baseIndex = 1;
    while session.(sensorLoc{i}).unixTimes(subj.(sensorLoc{i}).baseIndex) < subj.baseUnixTime
        subj.(sensorLoc{i}).baseIndex = subj.(sensorLoc{i}).baseIndex + 1;
    end
    subj.(sensorLoc{i}).baseIndex = subj.(sensorLoc{i}).baseIndex - 1;
    subj.(sensorLoc{i}).startIndex = subj.(sensorLoc{i}).baseIndex;
    while session.(sensorLoc{i}).unixTimes(subj.(sensorLoc{i}).startIndex) < subj.unixStartTime
        subj.(sensorLoc{i}).startIndex = subj.(sensorLoc{i}).startIndex + 1;
    end
    subj.(sensorLoc{i}).startIndex = subj.(sensorLoc{i}).startIndex - 1;
    subj.(sensorLoc{i}).endIndex = subj.(sensorLoc{i}).startIndex;
    while session.(sensorLoc{i}).unixTimes(subj.(sensorLoc{i}).endIndex) < subj.unixEndTime
        subj.(sensorLoc{i}).endIndex = subj.(sensorLoc{i}).endIndex + 1;
    end
    subj.(sensorLoc{i}).baseTimes = (raw.(sensorLoc{i})(subj.(sensorLoc{i}).baseIndex:subj.(sensorLoc{i}).startIndex,1) - raw.(sensorLoc{i})(subj.(sensorLoc{i}).baseIndex,1))/1000;
    subj.(sensorLoc{i}).times = (raw.(sensorLoc{i})(subj.(sensorLoc{i}).startIndex:subj.(sensorLoc{i}).endIndex,1) - raw.(sensorLoc{i})(subj.(sensorLoc{i}).startIndex,1))/1000;
end

%% Pulls Subject IMU data
for i = 1:length(sensorLoc)
    % Baseline IMU Data for subject
    subj.(sensorLoc{i}).baseAccel = raw.(sensorLoc{i})(subj.(sensorLoc{i}).baseIndex:subj.(sensorLoc{i}).startIndex,2:4);
    subj.(sensorLoc{i}).baseGyro = raw.(sensorLoc{i})(subj.(sensorLoc{i}).baseIndex:subj.(sensorLoc{i}).startIndex,6:8);
    subj.(sensorLoc{i}).baseMag = raw.(sensorLoc{i})(subj.(sensorLoc{i}).baseIndex:subj.(sensorLoc{i}).startIndex,9:11);
    % IMU Data for subject trial (m/s^2)
    subj.(sensorLoc{i}).accel = raw.(sensorLoc{i})(subj.(sensorLoc{i}).startIndex:subj.(sensorLoc{i}).endIndex,2:4);
    subj.(sensorLoc{i}).gyro = raw.(sensorLoc{i})(subj.(sensorLoc{i}).startIndex:subj.(sensorLoc{i}).endIndex,6:8);
    subj.(sensorLoc{i}).mag = raw.(sensorLoc{i})(subj.(sensorLoc{i}).startIndex:subj.(sensorLoc{i}).endIndex,9:11);
end

%% Filters subject IMU data
for i = 1:length(sensorLoc)
    for j = 1:3
        % Filter Baseline IMU Data for subject
        subj.(sensorLoc{i}).fBaseAccel(:,j) = filtfilt(session.filtB,session.filtA,subj.(sensorLoc{i}).baseAccel(:,j));
        subj.(sensorLoc{i}).fBaseGyro(:,j) = filtfilt(session.filtB,session.filtA,subj.(sensorLoc{i}).baseGyro(:,j));
        subj.(sensorLoc{i}).fBaseMag(:,j) = filtfilt(session.filtB,session.filtA,subj.(sensorLoc{i}).baseMag(:,j));
        % Filter trial IMU Data for subject
        subj.(sensorLoc{i}).fAccel(:,j) = filtfilt(session.filtB,session.filtA,subj.(sensorLoc{i}).accel(:,j));
        subj.(sensorLoc{i}).fGyro(:,j) = filtfilt(session.filtB,session.filtA,subj.(sensorLoc{i}).gyro(:,j));
        subj.(sensorLoc{i}).fMag(:,j) = filtfilt(session.filtB,session.filtA,subj.(sensorLoc{i}).mag(:,j));
    end
end

%% Finding angles of IMUs
for i = 1:length(sensorLoc)
    subj.(sensorLoc{i}).baseAngle(1,:) = subj.(sensorLoc{i}).fBaseGyro(1,:) .* subj.(sensorLoc{i}).baseTimes(1);
    for j = 1:3
        for k = 2:length(subj.(sensorLoc{i}).baseTimes)
            subj.(sensorLoc{i}).baseAngle(k,j) = ...
                subj.(sensorLoc{i}).baseAngle(k-1,j) + ...
                (subj.(sensorLoc{i}).fBaseGyro(k,j)+subj.(sensorLoc{i}).fBaseGyro(k-1,j))/2 * (subj.(sensorLoc{i}).baseTimes(k)-subj.(sensorLoc{i}).baseTimes(k-1));
        end
    end
end

% for i = 1:length(sensorLoc)
%     for j = 1:length(subj.(sensorLoc{i}).baseTimes)
%         if strcmp(sensorLoc{i},'waist')
%             subj.(sensorLoc{i}).baseAngle(j,1) = atand(-subj.(sensorLoc{i}).fBaseAccel(j,3)/sqrt((subj.(sensorLoc{i}).fBaseAccel(j,1)^2)+(subj.(sensorLoc{i}).fBaseAccel(j,2)^2)));
%             subj.(sensorLoc{i}).baseAngle(j,3) = atand(subj.(sensorLoc{i}).fBaseAccel(j,1)/subj.(sensorLoc{i}).fBaseAccel(j,2));
%             Mx = subj.(sensorLoc{i}).fBaseMag(j,3) * sind(subj.(sensorLoc{i}).baseAngle(j,3)) * sind(subj.(sensorLoc{i}).baseAngle(j,1)) + ...
%                  subj.(sensorLoc{i}).fBaseMag(j,1) * cosd(subj.(sensorLoc{i}).baseAngle(j,3)) - ...
%                  subj.(sensorLoc{i}).fBaseMag(j,2) * sind(subj.(sensorLoc{i}).baseAngle(j,3)) * cosd(subj.(sensorLoc{i}).baseAngle(j,1));
%             Mz = subj.(sensorLoc{i}).fBaseMag(j,3) * cosd(subj.(sensorLoc{i}).baseAngle(j,1)) + ...
%                  subj.(sensorLoc{i}).fBaseMag(j,2) * sind(subj.(sensorLoc{i}).baseAngle(j,1));
%             subj.(sensorLoc{i}).baseAngle(j,2) = atand(Mx, Mz);
%         else
%             subj.(sensorLoc{i}).angleBase(j,1) = atand(subj.(sensorLoc{i}).fBaseAccel(j,3)/subj.(sensorLoc{i}).fBaseAccel(j,2));
%             subj.(sensorLoc{i}).baseAngle(j,3) = atand(-subj.(sensorLoc{i}).fBaseAccel(j,1)/sqrt((subj.(sensorLoc{i}).fBaseAccel(j,2)^2)+(subj.(sensorLoc{i}).fBaseAccel(j,3)^2)));
%             Mx = subj.(sensorLoc{i}).fBaseMag(j,1) * cosd(subj.(sensorLoc{i}).baseAngle(j,3)) + ...
%                  subj.(sensorLoc{i}).fBaseMag(j,2) * sind(subj.(sensorLoc{i}).baseAngle(j,3));
%             MZ = subj.(sensorLoc{i}).fBaseMag(j,1) * sind(subj.(sensorLoc{i}).baseAngle(j,1)) * sind(subj.(sensorLoc{i}).baseAngle(j,3)) + ...
%                  subj.(sensorLoc{i}).fBaseMag(j,3) * cosd(subj.(sensorLoc{i}).baseAngle(j,1)) - ...
%                  subj.(sensorLoc{i}).fBaseMag(j,2) * sind(subj.(sensorLoc{i}).baseAngle(j,1)) * cosd(subj.(sensorLoc{i}).baseAngle(j,3));
%             subj.(sensorLoc{i}).baseAngle(j,2) = atand(Mz, Mx);
%         end
%     end
% end

for i = 1:length(sensorLoc)
    subj.(sensorLoc{i}).angle(1,:) = subj.(sensorLoc{i}).fGyro(1,:) .* subj.(sensorLoc{i}).times(1);
    for j = 1:3
        for k = 2:length(subj.(sensorLoc{i}).times)
            subj.(sensorLoc{i}).angle(k,j) = ...
                subj.(sensorLoc{i}).angle(k-1,j) + ...
                (subj.(sensorLoc{i}).fGyro(k,j)+subj.(sensorLoc{i}).fGyro(k-1,j))/2 * (subj.(sensorLoc{i}).times(k)-subj.(sensorLoc{i}).times(k-1));
        end
    end
end

% for i = 1:length(sensorLoc)
%     for j = 1:length(subj.(sensorLoc{i}).times)
%         if strcmp(sensorLoc{i},'waist')
%             subj.(sensorLoc{i}).angle(j,1) = atand(-subj.(sensorLoc{i}).fAccel(j,3)/sqrt((subj.(sensorLoc{i}).fAccel(j,1)^2)+(subj.(sensorLoc{i}).fAccel(j,2)^2)));
%             subj.(sensorLoc{i}).angle(j,3) = atand(subj.(sensorLoc{i}).fAccel(j,1)/subj.(sensorLoc{i}).fAccel(j,2));
%             Mx = subj.(sensorLoc{i}).fMag(j,3) * sind(subj.(sensorLoc{i}).angle(j,3)) * sind(subj.(sensorLoc{i}).angle(j,1)) + ...
%                  subj.(sensorLoc{i}).fMag(j,1) * cosd(subj.(sensorLoc{i}).angle(j,3)) - ...
%                  subj.(sensorLoc{i}).fMag(j,2) * sind(subj.(sensorLoc{i}).angle(j,3)) * cosd(subj.(sensorLoc{i}).angle(j,1));
%             Mz = subj.(sensorLoc{i}).fMag(j,3) * cosd(subj.(sensorLoc{i}).angle(j,1)) + ...
%                  subj.(sensorLoc{i}).fMag(j,2) * sind(subj.(sensorLoc{i}).angle(j,1));
%             subj.(sensorLoc{i}).angle(j,2) = atand(Mx, Mz);
%         else
%             subj.(sensorLoc{i}).angle(j,1) = atand(subj.(sensorLoc{i}).fAccel(j,3)/subj.(sensorLoc{i}).fAccel(j,2));
%             subj.(sensorLoc{i}).angle(j,3) = atand(-subj.(sensorLoc{i}).fAccel(j,1)/sqrt((subj.(sensorLoc{i}).fAccel(j,2)^2)+(subj.(sensorLoc{i}).fAccel(j,3)^2)));
%             Mx = subj.(sensorLoc{i}).fMag(j,1) * cosd(subj.(sensorLoc{i}).angle(j,3)) + ...
%                  subj.(sensorLoc{i}).fMag(j,2) * sind(subj.(sensorLoc{i}).angle(j,3));
%             MZ = subj.(sensorLoc{i}).fMag(j,1) * sind(subj.(sensorLoc{i}).angle(j,1)) * sind(subj.(sensorLoc{i}).angle(j,3)) + ...
%                  subj.(sensorLoc{i}).fMag(j,3) * cosd(subj.(sensorLoc{i}).angle(j,1)) - ...
%                  subj.(sensorLoc{i}).fMag(j,2) * sind(subj.(sensorLoc{i}).angle(j,1)) * cosd(subj.(sensorLoc{i}).angle(j,3));
%             subj.(sensorLoc{i}).angle(j,2) = atand(Mz, Mx);
%         end
%     end
% end

%% Baseline phase
for i = 1:length(sensorLoc)
   for j = 1:length(subj.(sensorLoc{i}).baseTimes)
       if strcmp(sensorLoc{i},'waist')
           subj.(sensorLoc{i}).fBaseAccel(j,2) = subj.(sensorLoc{i}).fBaseAccel(j,2) - gravity * cosd(subj.(sensorLoc{i}).baseAngle(j,1));
           subj.(sensorLoc{i}).fBaseAccel(j,3) = subj.(sensorLoc{i}).fBaseAccel(j,3) - gravity * sind(subj.(sensorLoc{i}).baseAngle(j,1));
       else
           subj.(sensorLoc{i}).fBaseAccel(j,2) = subj.(sensorLoc{i}).fBaseAccel(j,2) - gravity * cosd(subj.(sensorLoc{i}).baseAngle(j,3));
           subj.(sensorLoc{i}).fBaseAccel(j,1) = subj.(sensorLoc{i}).fBaseAccel(j,1) - gravity * sind(subj.(sensorLoc{i}).baseAngle(j,3));
       end
   end
end

for i = 1:length(sensorLoc)
    for j = 1:3
        subj.(sensorLoc{i}).avgBaseAccel(1,j) = sum(subj.(sensorLoc{i}).fBaseAccel(:,j))/length(subj.(sensorLoc{i}).baseTimes);
    end
end

%% Subject trial phase
for i = 1:length(sensorLoc)
   for j = 1:length(subj.(sensorLoc{i}).times)
       if strcmp(sensorLoc{i},'waist')
           subj.(sensorLoc{i}).fAccel(j,2) = subj.(sensorLoc{i}).fAccel(j,2) - gravity * cos(subj.(sensorLoc{i}).angle(j,1));
           subj.(sensorLoc{i}).fAccel(j,3) = subj.(sensorLoc{i}).fAccel(j,3) - gravity * sin(subj.(sensorLoc{i}).angle(j,1));
       else
           subj.(sensorLoc{i}).fAccel(j,2) = subj.(sensorLoc{i}).fAccel(j,2) - gravity * cos(subj.(sensorLoc{i}).angle(j,3));
           subj.(sensorLoc{i}).fAccel(j,1) = subj.(sensorLoc{i}).fAccel(j,1) - gravity * sin(subj.(sensorLoc{i}).angle(j,3));
       end
   end
end

for i = 1:length(sensorLoc)
    for j = 1:3
        subj.(sensorLoc{i}).fAccel(:,j) = subj.(sensorLoc{i}).fAccel(:,j) - subj.(sensorLoc{i}).avgBaseAccel(1,j);
    end
end

%% Integrate for velocities during trial
% for i = 1:length(sensorLoc)
%     for j = 1:3
%         subj.(sensorLoc{i}).vel(:,j) = cumtrapz(subj.waist.times, subj.waist.fAccel(:,j));
%     end
% end

for i = 1:length(sensorLoc)
    subj.(sensorLoc{i}).vel(1,:) = subj.(sensorLoc{i}).fAccel(1,:) .* subj.(sensorLoc{i}).times(1);
    for j = 1:3
        for k = 2:length(subj.(sensorLoc{i}).times)
            subj.(sensorLoc{i}).vel(k,j) = ...
                subj.(sensorLoc{i}).vel(k-1,j) + ...
                (subj.(sensorLoc{i}).fAccel(k,j)+subj.(sensorLoc{i}).fAccel(k-1,j))/2 * (subj.(sensorLoc{i}).times(k)-subj.(sensorLoc{i}).times(k-1));
        end
    end
end

%% Plotting of Baseline data
for i = 1:length(sensorLoc)
    figure(i)
    subplot(3,2,1)
    plot(subj.(sensorLoc{i}).baseTimes, subj.(sensorLoc{i}).baseAccel); 
    grid on; xlabel('seconds'); ylabel('m/s^2'); legend('X','Y','Z');
    title(['Subject ',num2str(subj.Num),' ',(sensorLoc{i}),' Baseline Accelerations'])
    
    subplot(3,2,3)
    plot(subj.(sensorLoc{i}).baseTimes, subj.(sensorLoc{i}).baseAngle); 
    grid on; xlabel('seconds'); ylabel('deg'); legend('X','Y','Z');
    title(['Subject ',num2str(subj.Num),' ',(sensorLoc{i}),' Baseline Angles'])
    
    subplot(3,2,2)
    plot(subj.(sensorLoc{i}).baseTimes, subj.(sensorLoc{i}).fBaseAccel); 
    grid on; xlabel('seconds'); ylabel('m/s^2'); legend('X','Y','Z');
    title(['Subject ',num2str(subj.Num),' ',(sensorLoc{i}),' filtered Baseline Accelerations'])
    
    subplot(3,2,4)
    plot(subj.(sensorLoc{i}).baseTimes, subj.(sensorLoc{i}).fBaseGyro);
    grid on; hold on; xlabel('seconds'); ylabel('deg/sec'); legend('X','Y','Z');
    title(['Subject ',num2str(subj.Num),' ',(sensorLoc{i}),' Baseline Angular Velocities'])
    
    subplot(3,2,5:6)
    plot(subj.(sensorLoc{i}).baseTimes, subj.(sensorLoc{i}).fBaseMag);
    grid on; xlabel('seconds'); ylabel('flux'); legend('X','Y','Z');
    title(['Subject ',num2str(subj.Num),' ',(sensorLoc{i}),' Baseline Local Flux'])
end

%% Plotting of Velocities and Accelerations
for i = 1:length(sensorLoc)
    figure(i+5)
    subplot(3,2,1)
    plot(subj.(sensorLoc{i}).times, subj.(sensorLoc{i}).vel); 
    grid on; xlabel('seconds'); ylabel('m/s'); legend('X','Y','Z');
    title(['Subject ',num2str(subj.Num),' ',(sensorLoc{i}),' Velocities'])
    
    subplot(3,2,3)
    plot(subj.(sensorLoc{i}).times, subj.(sensorLoc{i}).angle); 
    grid on; xlabel('seconds'); ylabel('deg'); legend('X','Y','Z');
    title(['Subject ',num2str(subj.Num),' ',(sensorLoc{i}),' Angles'])
    
    subplot(3,2,2)
    plot(subj.(sensorLoc{i}).times, subj.(sensorLoc{i}).fAccel); 
    grid on; xlabel('seconds'); ylabel('m/s^2'); legend('X','Y','Z');
    title(['Subject ',num2str(subj.Num),' ',(sensorLoc{i}),' Filtered Accelerations'])
    
    subplot(3,2,4)
    plot(subj.(sensorLoc{i}).times, subj.(sensorLoc{i}).accel); 
    grid on; xlabel('seconds'); ylabel('m/s^2'); legend('X','Y','Z');
    title(['Subject ',num2str(subj.Num),' ',(sensorLoc{i}),' Raw Accelerations'])
    
    subplot(3,2,5)
    plot(subj.(sensorLoc{i}).times, subj.(sensorLoc{i}).fGyro);
    grid on; xlabel('seconds'); ylabel('deg/sec'); legend('X','Y','Z');
    title(['Subject ',num2str(subj.Num),' ',(sensorLoc{i}),' Angular Velocities'])
    
    subplot(3,2,6)
    plot(subj.(sensorLoc{i}).times, subj.(sensorLoc{i}).fMag);
    grid on; xlabel('seconds'); ylabel('flux'); legend('X','Y','Z');
    title(['Subject ',num2str(subj.Num),' ',(sensorLoc{i}),' Local Flux'])
end
