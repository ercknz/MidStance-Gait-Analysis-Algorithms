function plotStepIMUdata(subjNum, location, steps, HSorMS, accXYZ, gyroXYZ, anglesXYZ, times, meanAccXYZ, meanGyroXYZ, meanAngleXYZ, stdAccXYZ, stdGyroXYZ, stdAngleXYZ, meanDur)
%% Plot IMU Mean Step Data
% This function takes in resampled IMU data thats has been segmented, 
% either Heel Strike to Heel Strike or Mid Stance to Mid Stance and plots 
% the profiles for global accelerations, gyroscope, and IMU angle data. 
% 
% The function assumes the following orientation of the data:
% +x: direction of walking
% +y: medial to Left Side
% +z: normal to ground (vertical)
% 
% Function by Erick Nunez

%% Heel Strike or Mid Stance
switch HSorMS
    case 'HS'
        stancePts = steps.HSindexes;
        stance = 'Heel-Strike';     stc = 'HS';
    case 'MS'
        stancePts = steps.MSindexes;
        stance = 'Mid-Stance';      stc = 'MS';
    otherwise
        error('Invalid step segmentation');
end

%% Plot the Data
offset = randi(1000);
% Figure 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig1 = figure(20 + subjNum + offset);
set(fig1, 'Units', 'Normalized', 'OuterPosition', [0,0, 1, 1]);
subplot(3,1,1)
plot(times, gyroXYZ(:,2)); grid on; hold on;
plot(times(stancePts), gyroXYZ(stancePts,2), 'rs') 
xlabel('seconds'); ylabel('deg/sec');
title([location,' Gyro w/ ', stance,' for Subject ',num2str(subjNum)])

subplot(3,1,2)
plot(times, accXYZ(:,3)); grid on; hold on;
plot(times(stancePts), accXYZ(stancePts,3), 'rs') 
xlabel('seconds'); ylabel('m/s^2');
title([location,' Global Vertical Acceleration w/ ', stance,' for Subject ',num2str(subjNum)])

subplot(3,1,3)
plot(times, anglesXYZ(:,2)); grid on; hold on;
plot(times(stancePts), anglesXYZ(stancePts,2), 'rs') 
xlabel('seconds'); ylabel('degrees');
title([location,' IMU Angle w/ ', stance,' for Subject ',num2str(subjNum)])

% Figure 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig2 = figure(40 + subjNum + offset);
set(fig2, 'Units', 'Normalized', 'OuterPosition', [0,0, 1, 1]);
subplot(3,2,1)
hold on; grid on;
for i = 1:length(steps.indexes)
    plot(gyroXYZ(steps.indexes{i},2))
end
xlabel('frames'); ylabel('deg/sec');
title([location,' Gyroscope Data of Steps for Subject ',num2str(subjNum)])

subplot(3,2,3)
hold on; grid on;
for i = 1:length(steps.indexes)
    plot(accXYZ(steps.indexes{i},3))
end
xlabel('frames'); ylabel('m/s^2');
title([location,' Global Vertical Acceleration of Steps for Subject',num2str(subjNum)])

subplot(3,2,5)
hold on; grid on;
for i = 1:length(steps.indexes)
    plot(anglesXYZ(steps.indexes{i},2))
end
xlabel('frames'); ylabel('degrees');
title([location,' IMU Angle of Steps for Subject ', num2str(subjNum)])

subplot(3,2,2)
errorbar(meanDur, meanGyroXYZ(:,2), stdGyroXYZ(:,2))
grid on;
xlabel('Sec'); ylabel('deg/sec');
title([location,' Mean Gyroscope Step Profile for Subject ',num2str(subjNum)])

subplot(3,2,4)
errorbar(meanDur, meanAccXYZ(:,3), stdAccXYZ(:,3))
grid on;
xlabel('Sec'); ylabel('m/s^2');
title([location,' Mean Global Vertical Accelerometer Step Profile for Subject ',num2str(subjNum)])

subplot(3,2,6)
errorbar(meanDur, meanAngleXYZ(:,2), stdAngleXYZ(:,2))
grid on;
xlabel('Sec'); ylabel('m/s^2');
title([location,' Mean IMU Angle Step Profile for Subject ',num2str(subjNum)])

%% Saves figures 
set(fig1, 'PaperOrientation', 'landscape', 'PaperUnits', 'normalized', 'PaperPosition',[0,0,1,1])
saveas(fig1, ['Shimmer Vicon Results/subj',num2str(subjNum),stance,location,'.pdf']);
set(fig2, 'PaperOrientation', 'landscape', 'PaperUnits', 'normalized', 'PaperPosition',[0,0,1,1])
saveas(fig2, ['Shimmer Vicon Results/subj',num2str(subjNum),'StepProfiles',location,stc,'.pdf']);

end