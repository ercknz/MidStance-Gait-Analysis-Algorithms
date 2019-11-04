function [waist, rightAnkle, leftAnkle] = gaitVelocity(subjNum, measuredWaist, waist, rightAnkle, leftAnkle)
%% GaitVelocity
% This functions takes in the measured waist data, the separated right
% steps, and the separated  left steps. Then, based on the leading foot, it
% pulls the waist accelerations per cycle (HS-HS). Next it calculates the
% velocities for the waist and ankles. Data is plotted and saved locally.
% Function assumed the following orientation of the data:
% +x: direction of walking
% +y: medial to Left Side
% +z: normal to ground (vertical)
%
% Function by erick nunez

%% Variables 

%% Uses leading Foot to find waist data.
if rightAnkle.steps.indexes{1}(1) < leftAnkle.steps.indexes{1}(1)
    waist.Step = imuDataResampling(rightAnkle.steps.indexes, ...
        measuredWaist.times, waist.globalAcc);
else
    waist.Step = imuDataResampling(leftAnkle.steps.indexes, ...
        measuredWaist.times, waist.globalAcc);
end

%% Finds Velocities
waist.Step.velocity =       calculateVel(waist.Step.meanAcc, mean(diff(waist.Step.avgTimes)));
rightAnkle.Step.velocity =  calculateVel(rightAnkle.Step.meanAcc, mean(diff(rightAnkle.Step.avgTimes)));
leftAnkle.Step.velocity =   calculateVel(leftAnkle.Step.meanAcc, mean(diff(leftAnkle.Step.avgTimes)));

%% Plot data 
fig1 = figure(subjNum);
set(fig1, 'Units', 'Normalized', 'OuterPosition', [0,0, 1, 1]);
subplot(3,1,1)
plot(waist.Step.avgTimes, waist.Step.velocity); grid on; 
xlabel('seconds'); ylabel('m/s');
title(['Velocity of Waist for Subject ',num2str(subjNum),' for Average Step'])

subplot(3,1,2)
plot(rightAnkle.Step.avgTimes, rightAnkle.Step.velocity); grid on; 
xlabel('seconds'); ylabel('m/s');
title(['Velocity of Right Ankle for Subject ',num2str(subjNum),' for Average Step'])

subplot(3,1,3)
plot(leftAnkle.Step.avgTimes, leftAnkle.Step.velocity); grid on; 
xlabel('seconds'); ylabel('m/s');
title(['Velocity of Left Ankle for Subject ',num2str(subjNum),' for Average Step'])

%% Saves figure and data
set(fig1, 'PaperOrientation', 'landscape', 'PaperUnits', 'normalized', 'PaperPosition',[0,0,1,1])
saveas(fig1, ['Shimmer Vicon Results/subj',num2str(subjNum),'Velocities.pdf']);

writematrix([mean(waist.Step.velocity(:,1)), mean(rightAnkle.Step.velocity(:,1)), mean(leftAnkle.Step.velocity(:,1))], 'Shimmer Vicon Results/shimmerViconVel.xlsx', 'Range', ['B',num2str(subjNum+2)])

end

function velocity = calculateVel(acc, dt)
velocity = zeros(size(acc));
for i = 2:length(velocity)
    velocity(i,:) = velocity(i-1,:) + (acc(i,:) + acc(i-1,:)).*dt/2;
end
end