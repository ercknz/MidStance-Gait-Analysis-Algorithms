function plotMoreIMUdata(testNumber, HSorMS, cRA, cLA, cW, cycle, yRemoved)
%% Plots More IMU Mean Step Data
% This function takes in IMU data thats has been segmented, either by Heel 
% Strike to Heel Strike or Mid Stance to Mid Stance and plots 
% the velocity, stride length, stride time, and double/single stance 
% information. 
% 
% The function assumes the following orientation of the data:
% +x: direction of walking
% +y: medial to Left Side
% +z: normal to ground (vertical)
% 
% Function by Erick Nunez

%% Variables to be used
switch HSorMS
    case 'HS'
%         stancePts = steps.HSindexes;
        stance = 'Heel-Strike';     stc = 'HS';
    case 'MS'
%         stancePts = steps.MSindexes;
        stance = 'Mid-Stance';      stc = 'MS';
    otherwise
        error('Invalid step segmentation');
end
subjNum = floor(testNumber/100);
testNum = testNumber-subjNum*100;

%% Figure 1 
% Figure 1 plots bar graphs showing the stride velocities, distances, and
% durations.
X = categorical({'Right Ankle','Left Ankle','Waist'});
X = reordercats(X,{'Right Ankle','Left Ankle','Waist'});
fig1 = figure(20 + subjNum);
set(fig1, 'Units', 'Normalized', 'OuterPosition', [0,0, 1, 1]);
subplot(3,12,1:4)
bar(X,[mean(cRA.steps.Distance),mean(cLA.steps.Distance),mean(cW.steps.Distance)]); 
grid on; hold on;
errorbar(X,[mean(cRA.steps.Distance),mean(cLA.steps.Distance),mean(cW.steps.Distance)],...
    [std(cRA.steps.Distance),std(cLA.steps.Distance),std(cW.steps.Distance)],'.',...
    'Capsize',25,'Linewidth',2) 
ylabel('Meters');
title([' Calculated Distance between ', stance,'s. Subject:',num2str(subjNum),' Trial:',num2str(testNum)])

subplot(3,12,5:8)
bar(X,[mean(cRA.steps.Duration),mean(cLA.steps.Duration),mean(cW.steps.Duration)]); 
grid on; hold on;
errorbar(X,[mean(cRA.steps.Duration),mean(cLA.steps.Duration),mean(cW.steps.Duration)],...
    [std(cRA.steps.Duration),std(cLA.steps.Duration),std(cW.steps.Duration)],'.',...
    'Capsize',25,'Linewidth',2) 
ylabel('Seconds');
title([' Measured Duration between ', stance,'s. Subject:',num2str(subjNum),' Trial:',num2str(testNum)])

subplot(3,12,9:12)
bar(X,[mean(cRA.steps.Speed),mean(cLA.steps.Speed),mean(cW.steps.Speed)]); 
grid on; hold on;
errorbar(X,[mean(cRA.steps.Speed),mean(cLA.steps.Speed),mean(cW.steps.Speed)],...
    [std(cRA.steps.Speed),std(cLA.steps.Speed),std(cW.steps.Speed)],'.',...
    'Capsize',25,'Linewidth',2) 
ylabel('Meters/Seconds');
title([' Calculated Walking Speed between ', stance,'s. Subject:',num2str(subjNum),' Trial:',num2str(testNum)])

%% Gait Cycle graph
subplot(3,12,13:18)
hold on; grid on;
errorbar(cycle.percent,cycle.rA.meanAcc(:,3),cycle.rA.stdAcc(:,3),'Color','#0072BD')
errorbar(cycle.percent,cycle.lA.meanAcc(:,3),cycle.lA.stdAcc(:,3),'Color','#D95319')
xline(cycle.percent(cycle.rA.HSIdx),'-.','Color','#0072BD');xline(cycle.percent(cycle.rA.TOIdx),'--','Color','#0072BD');
xline(cycle.percent(cycle.lA.HSIdx),'-.','Color','#D95319');xline(cycle.percent(cycle.lA.TOIdx),'--','Color','#D95319');
xlabel('Cycle Percentage'); ylabel('m/s');
title(['Subject:',num2str(subjNum),' Trial:',num2str(testNum),' Vertical Acceleration Profile'])
legend('rAnkle','lAnkle','Location','southeast')

subplot(3,12,19:24)
hold on; grid on;
errorbar(cycle.percent,cycle.rA.meanAcc(:,1),cycle.rA.stdAcc(:,1),'Color','#0072BD')
errorbar(cycle.percent,cycle.lA.meanAcc(:,1),cycle.lA.stdAcc(:,1),'Color','#D95319')
xline(cycle.percent(cycle.rA.HSIdx),'-.','Color','#0072BD');xline(cycle.percent(cycle.rA.TOIdx),'--','Color','#0072BD');
xline(cycle.percent(cycle.lA.HSIdx),'-.','Color','#D95319');xline(cycle.percent(cycle.lA.TOIdx),'--','Color','#D95319');
xlabel('Cycle Percentage'); ylabel('m/s');
title(['Subject:',num2str(subjNum),' Trial:',num2str(testNum),' Horizontal Acceleration Profile'])
legend('rAnkle','lAnkle','Location','southeast')

subplot(3,12,25:36)
hold on; grid on;
errorbar(cycle.percent,cycle.rA.meanGyro(:,2),cycle.rA.stdGyro(:,2),'Color','#0072BD')
errorbar(cycle.percent,cycle.lA.meanGyro(:,2),cycle.lA.stdGyro(:,2),'Color','#D95319')
xline(cycle.percent(cycle.rA.HSIdx),'-.','Color','#0072BD');xline(cycle.percent(cycle.rA.TOIdx),'--','Color','#0072BD');
xline(cycle.percent(cycle.lA.HSIdx),'-.','Color','#D95319');xline(cycle.percent(cycle.lA.TOIdx),'--','Color','#D95319');
xlabel('Cycle Percentage'); ylabel('Deg/Sec');
title(['Subject:',num2str(subjNum),' Trial:',num2str(testNum),' Sagital Angular Velocity Profile'])
legend('rAnkle','lAnkle','Location','southeast')

% subplot(3,12,31:36)
% hold on; grid on;
% errorbar(cycle.percent,cycle.rA.meanAng(:,2),cycle.rA.stdAng(:,2),'Color','#0072BD')
% errorbar(cycle.percent,cycle.lA.meanAng(:,2),cycle.lA.stdAng(:,2),'Color','#D95319')
% xline(cycle.percent(cycle.rA.HSIdx),'-.','Color','#0072BD');xline(cycle.percent(cycle.rA.TOIdx),'--','Color','#0072BD');
% xline(cycle.percent(cycle.lA.HSIdx),'-.','Color','#D95319');xline(cycle.percent(cycle.lA.TOIdx),'--','Color','#D95319');
% xlabel('Cycle Percentage'); ylabel('Degrees');
% title(['Subject:',num2str(subjNum),' Trial:',num2str(testNum),' Sagital Angles Profile'])
% legend('rAnkle','lAnkle','Location','southeast')

%% Figure 2
fig2 = figure(40 + subjNum);
set(fig2, 'Units', 'Normalized', 'OuterPosition', [0,0, 1, 1]);
subplot(3,2,1:2)
plot(cRA.times,cRA.gyro(:,2)); hold on; grid on;
plot(cRA.times(cRA.steps.MSwIdx),cRA.gyro(cRA.steps.MSwIdx,2),'rs');
plot(cRA.times(cRA.steps.TOIdx),cRA.gyro(cRA.steps.TOIdx,2),'bo');
plot(cRA.times(cRA.steps.HSIdx),cRA.gyro(cRA.steps.HSIdx,2),'m+');
plot(cRA.times(cRA.steps.MStIdx),cRA.gyro(cRA.steps.MStIdx,2),'go','Markersize',5);
xlabel('Time(Sec)'); ylabel('Deg/Sec');xlim([20,40]);
title(['Subject:',num2str(subjNum),' Trial:',num2str(testNum),' Right Ankle Ang Vel with MSt, TO, HS, and MSw'])

subplot(3,2,3:4)
plot(cLA.times,cLA.gyro(:,2)); hold on; grid on;
plot(cLA.times(cLA.steps.MSwIdx),cLA.gyro(cLA.steps.MSwIdx,2),'rs');
plot(cLA.times(cLA.steps.TOIdx),cLA.gyro(cLA.steps.TOIdx,2),'bo');
plot(cLA.times(cLA.steps.HSIdx),cLA.gyro(cLA.steps.HSIdx,2),'m+');
plot(cLA.times(cLA.steps.MStIdx),cLA.gyro(cLA.steps.MStIdx,2),'go','Markersize',5);
xlabel('Time(Sec)'); ylabel('Deg/Sec');xlim([20,40]);
title(['Subject:',num2str(subjNum),' Trial:',num2str(testNum),' Left Ankle Ang Vel with MSt, TO, HS, and MSw'])

subplot(3,2,5)
hold on; grid on;
for i = 1:length(cRA.steps.indexes)
    plot(cRA.gyro(cRA.steps.indexes{i},2))
end
xlabel('Frames'); ylabel('Deg/Sec');
title(['Subject:',num2str(subjNum),' Trial:',num2str(testNum),' Right Ankle Steps'])

subplot(3,2,6)
hold on; grid on;
for i = 1:length(cLA.steps.indexes)
    plot(cLA.gyro(cLA.steps.indexes{i},2))
end
xlabel('Frames'); ylabel('Deg/Sec');
title(['Subject:',num2str(subjNum),' Trial:',num2str(testNum),' Right Ankle Steps'])

%% Saves figures 
set(fig1, 'PaperOrientation', 'portrait', 'PaperUnits', 'normalized', 'PaperPosition',[0,0,1,1]);
if yRemoved
    saveas(fig1, ['Shimmer Vicon Results/subj',num2str(subjNum),'Trial',num2str(testNum),'Cycle',stc,'Yaxis.pdf']);
else
    saveas(fig1, ['Shimmer Vicon Results/subj',num2str(subjNum),'Trial',num2str(testNum),'Cycle',stc,'.pdf']);
end


end