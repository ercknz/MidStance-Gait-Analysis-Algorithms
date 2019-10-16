%% Clean workspace
clear; clc; close all;

%% Creates c3dData object
subjNumber = 1;
fileName = ['C3D Files/Gait0',num2str(subjNumber),'.c3d'];
addpath('Gait Data Classes');
gaitData = viconC3Ddata(fileName);

%% Findings data points
x = (gaitData.data(:,1,3)+gaitData.data(:,1,4))/2;
y = (gaitData.data(:,2,3)+gaitData.data(:,2,4))/2;
z = (gaitData.data(:,3,3)+gaitData.data(:,3,4))/2;
waist = [x,y,z];
rightAnkle = gaitData.data(:,:,14);
leftAnkle = gaitData.data(:,:,8);

[waistVel, waistSpeed] = viconPointVel(waist, gaitData.times);
[rightAnkleVel, rightAnkleSpeed] = viconPointVel(rightAnkle, gaitData.times);
[leftAnkleVel, leftAnkleSpeed] = viconPointVel(leftAnkle, gaitData.times);

%% plot data
if ~gaitData.dataFound
    for i=1:gaitData.numberOfFrames
        tic; cla;
        axis([-1,2,-3.5,3.5,0,2.4])
        xlabel('X'); ylabel('Y'); zlabel('Z');
        hold on; grid on;
        for j =1:16
            plot3(gaitData.data(i,1,j), gaitData.data(i,2,j), gaitData.data(i,3,j), '.', 'MarkerSize',10);   
        end
        quiver3(waist(i,1), waist(i,2), waist(i,3), waistVel(i,1), waistVel(i,2), waistVel(i,3), 'o', 'MarkerSize',5)
        quiver3(rightAnkle(i,1), rightAnkle(i,2), rightAnkle(i,3), rightAnkleVel(i,1), rightAnkleVel(i,2), rightAnkleVel(i,3), 'o', 'MarkerSize',5)
        quiver3(leftAnkle(i,1), leftAnkle(i,2), leftAnkle(i,3), leftAnkleVel(i,1), leftAnkleVel(i,2), leftAnkleVel(i,3), 'o', 'MarkerSize',5)
        while toc < 1/gaitData.frameRate
        end
        drawnow
    end
else
    fig1 = figure(300+subjNumber);    
    set(fig1, 'Units', 'Normalized', 'OuterPosition', [0,0, 1, 1]);
    plot(gaitData.times,waistSpeed)
    grid on; legend('Waist')
    title(['Speed of Points from Subject ',num2str(subjNumber)])
    xlabel('Times (sec)'); ylabel('Speed (m/s)');
end

%% Save Workspace and figure
set(fig1, 'PaperOrientation', 'landscape', 'PaperUnits', 'normalized', 'PaperPosition',[0,0,1,1])
saveas(fig1, ['Shimmer Vicon Results/waistSpeedSubject',num2str(subjNumber),'.pdf']);
save(['Shimmer Workspaces/shimmerViconC3D',num2str(subjNumber),'.mat'])
writematrix([mean(waistSpeed), mean(rightAnkleSpeed), mean(leftAnkleSpeed)], 'Shimmer Vicon Results/shimmerViconVel.xlsx', 'Range', ['E',num2str(subjNumber+2)])
