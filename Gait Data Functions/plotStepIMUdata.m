function plotStepIMUdata(subjNum, location, HSindexes, accXYZ, gyroXYZ, anglesXYZ, times, indexes, meanAccXYZ, meanGyroXYZ, meanAngleXYZ, stdAccXYZ, stdGyroXYZ, stdAngleXYZ, meanDur)
%% Plot IMU Mean Step Data
% This function takes in resampled IMU data thats has been segmented heel 
% strike to heel strike and plots the profiles for global accelerations, 
% gyroscope, and IMU angle data. 
% 
% The function assumes the following orientation of the data:
% +x: direction of walking
% +y: medial-to-lateral
% +z: lower-to-upper body
% 
% Function by Erick Nunez

%% Plot the Data
offset = randi(1000);
% Figure 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig1 = figure(20 + subjNum + offset);
set(fig1, 'Units', 'Normalized', 'OuterPosition', [0,0, 1, 1]);
subplot(3,1,1)
plot(times, gyroXYZ(:,2),'-b'); grid on; hold on;
plot(times(HSindexes), gyroXYZ(HSindexes,2), 'rs') 
xlabel('seconds'); ylabel('deg/sec');
title([location,' Gyro w/ Heel-Strikes for Subject ',num2str(subjNum)])

subplot(3,1,2)
plot(times, accXYZ(:,3),'-b'); grid on; hold on;
plot(times(HSindexes), accXYZ(HSindexes,3), 'rs') 
xlabel('seconds'); ylabel('m/s^2');
title([location,' Global Acceleration w/ Heel-Strikes for Subject ',num2str(subjNum)])

subplot(3,1,3)
plot(times, anglesXYZ(:,2),'-b'); grid on; hold on;
plot(times(HSindexes), anglesXYZ(HSindexes,2), 'rs') 
xlabel('seconds'); ylabel('degrees');
title([location,' IMU Angle w/ Heel-Strikes for Subject ',num2str(subjNum)])

% Figure 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig2 = figure(40 + subjNum + offset);
set(fig2, 'Units', 'Normalized', 'OuterPosition', [0,0, 1, 1]);
subplot(3,2,1)
hold on; grid on;
for i = 1:length(indexes)
    plot(gyroXYZ(indexes{i},2))
end
xlabel('frames'); ylabel('deg/sec');
title([location,' Gyroscope Data of Steps for Subject ',num2str(subjNum)])

subplot(3,2,3)
hold on; grid on;
for i = 1:length(indexes)
    plot(accXYZ(indexes{i},1))
end
xlabel('frames'); ylabel('m/s^2');
title([location,' Global Acceleration of Steps for Subject',num2str(subjNum)])

subplot(3,2,5)
hold on; grid on;
for i = 1:length(indexes)
    plot(anglesXYZ(indexes{i},2))
end
xlabel('frames'); ylabel('degrees');
title([location,' IMU Angle of Steps for Subject ', num2str(subjNum)])

subplot(3,2,2)
errorbar(meanDur, meanGyroXYZ(:,2), stdGyroXYZ(:,2))
grid on;
xlabel('Sec'); ylabel('deg/sec');
title([location,' Mean Gyroscope Step Profile for Subject ',num2str(subjNum)])

subplot(3,2,4)
errorbar(meanDur, meanAccXYZ(:,1), stdAccXYZ(:,1))
grid on;
xlabel('Sec'); ylabel('m/s^2');
title([location,' Mean Global Accelerometer Step Profile for Subject ',num2str(subjNum)])

subplot(3,2,6)
errorbar(meanDur, meanAngleXYZ(:,2), stdAngleXYZ(:,2))
grid on;
xlabel('Sec'); ylabel('m/s^2');
title([location,' Mean IMU Angle Step Profile for Subject ',num2str(subjNum)])

%% Saves figures 
set(fig1, 'PaperOrientation', 'landscape', 'PaperUnits', 'normalized', 'PaperPosition',[0,0,1,1])
saveas(fig1, ['Shimmer Vicon Results/subj',num2str(subjNum),'HeelStrike',location,'.pdf']);
set(fig2, 'PaperOrientation', 'landscape', 'PaperUnits', 'normalized', 'PaperPosition',[0,0,1,1])
saveas(fig2, ['Shimmer Vicon Results/subj',num2str(subjNum),'StepProfiles',location,'.pdf']);

end