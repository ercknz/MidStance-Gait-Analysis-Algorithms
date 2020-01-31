function plotStepIMUdata(testNumber, location, steps, HSorMS, accXYZ, gyroXYZ, anglesXYZ, times, meanAccXYZ, meanGyroXYZ, meanAngleXYZ, stdAccXYZ, stdGyroXYZ, stdAngleXYZ, meanDur)
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

%% Variables
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
subjNum = floor(testNumber/100);
testNum = testNumber-subjNum*100;

%% Plot the Data
% Figure 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig1 = figure(20 + subjNum + randi(100));
set(fig1, 'Units', 'Normalized', 'OuterPosition', [0,0, 1, 1]);
subplot(4,1,1)
plot(times, gyroXYZ(:,2)); grid on; hold on;
plot(times(stancePts), gyroXYZ(stancePts,2), 'rs') 
xlabel('seconds'); ylabel('deg/sec'); xlim([20,40]);
title([location,' Angular Vel w/ ', stance,'. Subject:',num2str(subjNum),' Trial:',num2str(testNum)])

subplot(4,1,2)
plot(times, accXYZ(:,3)); grid on; hold on;
plot(times(stancePts), accXYZ(stancePts,3), 'rs') 
xlabel('seconds'); ylabel('m/s^2'); xlim([20,40]);
title([location,' Global Vertical Acc w/ ', stance,'. Subject:',num2str(subjNum),' Trial:',num2str(testNum)])

subplot(4,1,3)
plot(times, accXYZ(:,1)); grid on; hold on;
plot(times(stancePts), accXYZ(stancePts,1), 'rs') 
xlabel('seconds'); ylabel('m/s^2'); xlim([20,40]);
title([location,' Global Horizontal Acc w/ ', stance,'. Subject:',num2str(subjNum),' Trial:',num2str(testNum)])

subplot(4,1,4)
plot(times, anglesXYZ(:,2)); grid on; hold on;
plot(times(stancePts), anglesXYZ(stancePts,2), 'rs') 
xlabel('seconds'); ylabel('degrees'); xlim([20,40]);
title([location,' Angle w/ ', stance,'. Subject:',num2str(subjNum),' Trial:',num2str(testNum)])

% Figure 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig2 = figure(40 + subjNum + randi(100));
set(fig2, 'Units', 'Normalized', 'OuterPosition', [0,0, 1, 1]);
subplot(4,2,1)
hold on; grid on;
for i = 1:length(steps.indexes)
    plot(gyroXYZ(steps.indexes{i},2))
end
xlabel('frames'); ylabel('deg/sec');
title([location,' Angular Vel of Steps. Subject:',num2str(subjNum),' Trial:',num2str(testNum)])

subplot(4,2,3)
hold on; grid on;
for i = 1:length(steps.indexes)
    plot(accXYZ(steps.indexes{i},3))
end
xlabel('frames'); ylabel('m/s^2');
title([location,' Global Vertical Acc of Steps. Subject:',num2str(subjNum),' Trial:',num2str(testNum)])

subplot(4,2,5)
hold on; grid on;
for i = 1:length(steps.indexes)
    plot(accXYZ(steps.indexes{i},1))
end
xlabel('frames'); ylabel('m/s^2');
title([location,' Global Horizontal Acc of Steps. Subject:',num2str(subjNum),' Trial:',num2str(testNum)])

subplot(4,2,7)
hold on; grid on;
for i = 1:length(steps.indexes)
    plot(anglesXYZ(steps.indexes{i},2))
end
xlabel('frames'); ylabel('degrees');
title([location,' Angle of Steps. Subject:',num2str(subjNum),' Trial:',num2str(testNum)])

subplot(4,2,2)
errorbar(meanDur, meanGyroXYZ(:,2), stdGyroXYZ(:,2))
grid on;
xlabel('Sec'); ylabel('deg/sec');
title([location,' Angular Vel Step Profile. Subject:',num2str(subjNum),' Trial:',num2str(testNum)])

subplot(4,2,4)
errorbar(meanDur, meanAccXYZ(:,3), stdAccXYZ(:,3))
grid on;
xlabel('Sec'); ylabel('m/s^2');
title([location,' Global Vertical Acc Step Profile. Subject:',num2str(subjNum),' Trial:',num2str(testNum)])

subplot(4,2,6)
errorbar(meanDur, meanAccXYZ(:,1), stdAccXYZ(:,1))
grid on;
xlabel('Sec'); ylabel('m/s^2');
title([location,' Global Horizontal Acc Step Profile. Subject:',num2str(subjNum),' Trial:',num2str(testNum)])

subplot(4,2,8)
errorbar(meanDur, meanAngleXYZ(:,2), stdAngleXYZ(:,2))
grid on;
xlabel('Sec'); ylabel('degrees');
title([location,' Angle Step Profile. Subject:',num2str(subjNum),' Trial:',num2str(testNum)])

%% Saves figures 
set(fig1, 'PaperOrientation', 'portrait', 'PaperUnits', 'normalized', 'PaperPosition',[0,0,1,1]);
saveas(fig1, ['Shimmer Vicon Results/subj',num2str(subjNum),'Trial',num2str(testNum),stance,location,'.pdf']);
set(fig2, 'PaperOrientation', 'portrait', 'PaperUnits', 'normalized', 'PaperPosition',[0,0,1,1]);
saveas(fig2, ['Shimmer Vicon Results/subj',num2str(subjNum),'Trial',num2str(testNum),'Steps',location,stc,'.pdf']);

end