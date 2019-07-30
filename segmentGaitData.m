function subj = segmentGaitData(subj, timesSigma, yAxis)
%% SegmentGaitData
% This function is used to seperate the data angular velocity data into
% steps so we can normalize the data afterwards with respect to the gait
% stance.
%
% timesSigma should be a value between 1-3 to remove outliers. 
% yAxis is a boolean value to enable to removal of yAxis outliers. 
% 
% Function by Erick Nunez

%% Variables set up.
sensorLoc = {'leftAnkle','rightAnkle','waist'};
side = {'left','right'};
subj.right.TOidx = [];      subj.left.TOidx = [];
subj.right.steps = {};      subj.left.steps = {};
subj.right.accSteps = {};   subj.left.accSteps = {};
subj.right.stepTimes = {};  subj.left.stepTimes = {};
subj.right.stepDur = [];    subj.left.stepDur = [];
subj.right.stepVecLen = []; subj.left.stepVecLen = [];
subj.right.RSsteps = {};    subj.left.RSsteps = {};
subj.right.avgStep = [];    subj.left.avgStep = [];
subj.right.stdStep = [];    subj.left.stdStep = [];
subj.right.RSAccSteps = {}; subj.left.RSAccSteps = {};
subj.right.avgAccStep = []; subj.left.avgAccStep = [];
subj.right.stdAccStep = []; subj.left.stdAccStep = [];
subj.right.senseTimes = subj.rightAnkle.times;
subj.left.senseTimes = subj.leftAnkle.times;
subj.waist.cycleAccX = {};   subj.waist.cycleAccY = {};
subj.waist.cycleTimes = {};  subj.waist.cycleDur = [];   
subj.waist.cycleVecLen = [];


%% Finds Mins
subj.left.minsIdx = find(islocalmin(subj.leftAnkle.fGyro(:,3),'MinSeparation',0.4,'SamplePoints',subj.leftAnkle.times));
subj.right.minsIdx = find(islocalmin(subj.rightAnkle.fGyro(:,3),'MinSeparation',0.4,'SamplePoints',subj.rightAnkle.times));
subj.left.minsVal = subj.leftAnkle.fGyro(subj.left.minsIdx,3);
subj.right.minsVal = subj.rightAnkle.fGyro(subj.right.minsIdx,3);

%% Find toeOffs
for i = 1:length(side)
   for j = 2:length(subj.(side{i}).minsVal)-1
       if (subj.(side{i}).minsVal(j) < subj.(side{i}).minsVal(j-1)) && (subj.(side{i}).minsVal(j) < subj.(side{i}).minsVal(j+1))
           subj.(side{i}).TOidx(end+1) = subj.(side{i}).minsIdx(j);
       end
   end
end

% finds leading side
if subj.rightAnkle.times(subj.right.TOidx(1)) < subj.leftAnkle.times(subj.left.TOidx(1))
    lagSide = 'left'; leadSide = 'right';
else
    lagSide = 'right'; leadSide = 'left';
end


%% Seperate steps and times
for i = 1:length(side)
    for j = 1:length(subj.(side{i}).TOidx)-1
        subj.(side{i}).steps{end+1} = subj.(sensorLoc{i}).fGyro(subj.(side{i}).TOidx(j):subj.(side{i}).TOidx(j+1),3);
        subj.(side{i}).accSteps{end+1} = subj.(sensorLoc{i}).globalAcc(subj.(side{i}).TOidx(j):subj.(side{i}).TOidx(j+1),1);
        subj.(side{i}).stepTimes{end+1} = subj.(sensorLoc{i}).times(subj.(side{i}).TOidx(j):subj.(side{i}).TOidx(j+1));
        subj.(side{i}).stepDur(end+1) = subj.(side{i}).stepTimes{j}(end) - subj.(side{i}).stepTimes{j}(1);
        subj.(side{i}).stepVecLen(end+1) = length(subj.(side{i}).steps{j});
    end
    subj.(side{i}).avgStepVecLen = mean(subj.(side{i}).stepVecLen);
    subj.(side{i}).stdStepVecLen = std(subj.(side{i}).stepVecLen);
    subj.(side{i}).avgStepTime = mean(subj.(side{i}).stepDur);
    subj.(side{i}).stdStepTime = std(subj.(side{i}).stepDur);
    subj.(side{i}).minVecLen = floor(subj.(side{i}).avgStepVecLen - subj.(side{i}).stdStepVecLen);
end

%% Finds data of lagging foot
subj.(lagSide).lagIdx = zeros(1,length(subj.(leadSide).TOidx));
for i = 1:length(subj.(leadSide).TOidx)
    [value, index] = min(abs(subj.(lagSide).senseTimes - subj.(leadSide).senseTimes(subj.(leadSide).TOidx(i))));
    subj.(lagSide).lagIdx(i) = index;
end

%% Finds waist data w/r to leading foot
subj.waist.cycIdx = zeros(1,length(subj.(leadSide).TOidx));
for i = 1:length(subj.(leadSide).TOidx)
    [value, index] = min(abs(subj.waist.times - subj.(leadSide).senseTimes(subj.(leadSide).TOidx(i))));
    subj.waist.cycIdx(i) = index;
end

%% Seperate waist data and times
for j = 1:length(subj.waist.cycIdx)-1
    subj.waist.cycleAccX{end+1} = subj.waist.globalAcc(subj.waist.cycIdx(j):subj.waist.cycIdx(j+1),1);
    subj.waist.cycleAccY{end+1} = subj.waist.globalAcc(subj.waist.cycIdx(j):subj.waist.cycIdx(j+1),2);
    subj.waist.cycleTimes{end+1} = subj.waist.times(subj.waist.cycIdx(j):subj.waist.cycIdx(j+1));
    subj.waist.cycleDur(end+1) = subj.waist.cycleTimes{j}(end) - subj.waist.cycleTimes{j}(1);
    subj.waist.cycleVecLen(end+1) = length(subj.waist.cycleAccX{j});
end
subj.waist.avgCycleVecLen = mean(subj.waist.cycleVecLen);
subj.waist.stdCycleVecLen = std(subj.waist.cycleVecLen);
subj.waist.avgCycleTime = mean(subj.waist.cycleDur);
subj.waist.stdCycleTime = std(subj.waist.cycleDur);
subj.waist.minCycleVecLen = floor(subj.waist.avgCycleVecLen - subj.waist.stdCycleVecLen);


%% Removes outlier steps
for i = 1:length(side)
    removed = 0;
    for j = 1:length(subj.(side{i}).stepDur)
        if (subj.(side{i}).stepDur(j-removed) > (subj.(side{i}).avgStepTime + timesSigma*subj.(side{i}).stdStepTime)) || (subj.(side{i}).stepDur(j-removed) < (subj.(side{i}).avgStepTime - timesSigma*subj.(side{i}).stdStepTime))
            subj.(side{i}).steps(j-removed) = [];
            subj.(side{i}).accSteps(j-removed) = [];
            subj.(side{i}).stepTimes(j-removed) = [];
            subj.(side{i}).stepDur(j-removed) = [];
            subj.(side{i}).stepVecLen(j-removed) = [];
            removed = removed + 1;
        end
    end
end
removed = 0;
for j = 1:length(subj.waist.cycleDur)
    if (subj.waist.cycleDur(j-removed) > (subj.waist.avgCycleTime + timesSigma*subj.waist.stdCycleTime)) || (subj.waist.cycleDur(j-removed) < (subj.waist.avgCycleTime - timesSigma*subj.waist.stdCycleTime))
        subj.waist.cycleAccX(j-removed) = [];
        subj.waist.cycleAccY(j-removed) = [];
        subj.waist.cycleTimes(j-removed) = [];
        subj.waist.cycleDur(j-removed) = [];
        subj.waist.cycleVecLen(j-removed) = [];
        removed = removed + 1;
    end
end

%% Resamples the steps
for i = 1:length(side)
   for j = 1:length(subj.(side{i}).steps)
      subj.(side{i}).RSsteps{end+1} = resample(subj.(side{i}).steps{j},subj.(side{i}).minVecLen,length(subj.(side{i}).steps{j}));
      subj.(side{i}).RSAccSteps{end+1} = resample(subj.(side{i}).accSteps{j},subj.(side{i}).minVecLen,length(subj.(side{i}).accSteps{j}));
   end
end
subj.waist.RScycleAccX = {};    subj.waist.RScycleAccY = {};
for j = 1:length(subj.waist.cycleAccX)
    subj.waist.RScycleAccX{end+1} = resample(subj.waist.cycleAccX{j},subj.waist.minCycleVecLen,length(subj.waist.cycleAccX{j}));
    subj.waist.RScycleAccY{end+1} = resample(subj.waist.cycleAccY{j},subj.waist.minCycleVecLen,length(subj.waist.cycleAccY{j}));
end

%% finds average and std per step.
for i = 1:length(side)
   for j = 1:subj.(side{i}).minVecLen
       samples = []; accSamples = [];
       for k = 1:length(subj.(side{i}).RSsteps)
           samples(end+1) = subj.(side{i}).RSsteps{k}(j);
           accSamples(end+1) = subj.(side{i}).RSAccSteps{k}(j);
       end
       subj.(side{i}).avgStep(end+1) = mean(samples);
       subj.(side{i}).stdStep(end+1) = std(samples);
       subj.(side{i}).avgAccStep(end+1) = mean(accSamples);
       subj.(side{i}).stdAccStep(end+1) = std(accSamples);
   end
end
subj.waist.avgCycleAccX = [];   subj.waist.avgCycleAccY = [];
subj.waist.stdCycleAccX = [];   subj.waist.stdCycleAccY = [];
for j = 1:subj.waist.minCycleVecLen
    samplesX = []; samplesY = [];
    for k = 1:length(subj.waist.RScycleAccX)
        samplesX(end+1) = subj.waist.RScycleAccX{k}(j);
        samplesY(end+1) = subj.waist.RScycleAccY{k}(j);
    end
    subj.waist.avgCycleAccX(end+1) = mean(samplesX);
    subj.waist.stdCycleAccX(end+1) = std(samplesX);
    subj.waist.avgCycleAccY(end+1) = mean(samplesY);
    subj.waist.stdCycleAccY(end+1) = std(samplesY);
end

%% Removes outliers steps based on y-axis
if yAxis
    for i = 1:length(side)
        removed = 0;
        for j = 1:length(subj.(side{i}).RSAccSteps)
            for k = 1:subj.(side{i}).minVecLen
                if (subj.(side{i}).RSAccSteps{j-removed}(k) < (subj.(side{i}).avgAccStep(k) - timesSigma * subj.(side{i}).stdAccStep(k))) || (subj.(side{i}).RSAccSteps{j-removed}(k) > (subj.(side{i}).avgAccStep(k) + timesSigma * subj.(side{i}).stdAccStep(k)))
                    subj.(side{i}).RSAccSteps(j-removed) = [];
                    subj.(side{i}).RSsteps(j-removed) = [];
                    removed = removed + 1;
                    break
                end
            end
        end
    end
    removed = 0;
    for j = 1:length(subj.waist.RScycleAccX)
        for k = 1:subj.waist.minCycleVecLen
            if (subj.waist.RScycleAccX{j-removed}(k) < (subj.waist.avgCycleAccX(k) - timesSigma * subj.waist.stdCycleAccX(k))) || (subj.waist.RScycleAccX{j-removed}(k) > (subj.waist.avgCycleAccX(k) + timesSigma * subj.waist.stdCycleAccX(k)))
                subj.waist.RScycleAccX(j-removed) = [];
                subj.waist.RScycleAccY(j-removed) = [];
                removed = removed + 1;
                break
            end
        end
    end
end

%% If yAxis outlaiers removed, then Recalculate average and std per step.
if yAxis
    for i = 1:length(side)
        subj.(side{i}).avgStep = [];        subj.(side{i}).stdStep = [];
        subj.(side{i}).avgAccStep = [];     subj.(side{i}).stdAccStep = [];
        for j = 1:subj.(side{i}).minVecLen
            samples = []; accSamples = [];
            for k = 1:length(subj.(side{i}).RSsteps)
                samples(end+1) = subj.(side{i}).RSsteps{k}(j);
                accSamples(end+1) = subj.(side{i}).RSAccSteps{k}(j);
            end
            subj.(side{i}).avgStep(end+1) = mean(samples);
            subj.(side{i}).stdStep(end+1) = std(samples);
            subj.(side{i}).avgAccStep(end+1) = mean(accSamples);
            subj.(side{i}).stdAccStep(end+1) = std(accSamples);
        end
    end
    subj.waist.avgCycleAccX = [];   subj.waist.avgCycleAccY = [];
    subj.waist.stdCycleAccX = [];   subj.waist.stdCycleAccY = [];
    for j = 1:subj.waist.minCycleVecLen
        samplesX = []; samplesY = [];
        for k = 1:length(subj.waist.RScycleAccX)
            samplesX(end+1) = subj.waist.RScycleAccX{k}(j);
            samplesY(end+1) = subj.waist.RScycleAccY{k}(j);
        end
        subj.waist.avgCycleAccX(end+1) = mean(samplesX);
        subj.waist.stdCycleAccX(end+1) = std(samplesX);
        subj.waist.avgCycleAccY(end+1) = mean(samplesY);
        subj.waist.stdCycleAccY(end+1) = std(samplesY);
    end
end

%% Velocity of average waist acceleration
subj.waist.cycleTimes = linspace(0,subj.waist.avgCycleTime,subj.waist.minCycleVecLen);
subj.waist.cycleVel(1) = subj.waist.avgCycleAccX(1) .* subj.waist.cycleTimes(1);
for i = 2:length(subj.waist.cycleTimes)
    subj.waist.cycleVel(i) = ...
        subj.waist.cycleVel(i-1) + ...
        (subj.waist.avgCycleAccX(i)+subj.waist.avgCycleAccX(i-1))/2 * (subj.waist.cycleTimes(i)-subj.waist.cycleTimes(i-1));
end
subj.waist.cycleAvgVel = mean(subj.waist.cycleVel);
mean(subj.waist.cycleVel)

%% Velocity of average Ankle Accelerations
for i = 1:length(side)
    subj.(side{i}).avgStepTimes = linspace(0,subj.(side{i}).avgStepTime,subj.(side{i}).minVecLen);
    subj.(side{i}).stepVel(1) = subj.(side{i}).avgAccStep(1) .* subj.(side{i}).avgStepTimes(1);
    for j = 2:length(subj.(side{i}).avgStepTimes)
        subj.(side{i}).stepVel(j) = ...
            subj.(side{i}).stepVel(j-1) + ...
            (subj.(side{i}).avgAccStep(j)+subj.(side{i}).avgAccStep(j-1))/2 * (subj.(side{i}).avgStepTimes(j)-subj.(side{i}).avgStepTimes(j-1));
    end
    subj.(side{i}).avgVel = mean(subj.(side{i}).stepVel);
    mean(subj.(side{i}).stepVel)
end

%% Plots gait phases
close all
% Figure 1 shows gyro data used to seperate cycles %%%%%%%%%%%%%%%%%%%%%%%%
fig1 = figure(100+subj.Num);    
set(fig1, 'Units', 'Normalized', 'OuterPosition', [0,0, 0.5, 1]);
subplot(3,1,1)
plot(subj.leftAnkle.times, subj.leftAnkle.fGyro(:,3),'-b')
grid on; hold on;
plot(subj.rightAnkle.times, subj.rightAnkle.fGyro(:,3),'-r')
plot(subj.leftAnkle.times(subj.left.minsIdx), subj.leftAnkle.fGyro(subj.left.minsIdx,3),'*g')
plot(subj.rightAnkle.times(subj.right.minsIdx), subj.rightAnkle.fGyro(subj.right.minsIdx,3),'*m')
plot(subj.leftAnkle.times(subj.left.TOidx), subj.leftAnkle.fGyro(subj.left.TOidx,3),'bs')
plot(subj.rightAnkle.times(subj.right.TOidx), subj.rightAnkle.fGyro(subj.right.TOidx,3),'bs')
xlabel('seconds'); ylabel('deg/sec'); legend('left','right');
title(['Subject ',num2str(subj.Num),' Angular Vel'])

subplot(3,1,2)
plot(subj.right.RSsteps{1});
hold on; grid on;
for i=1:length(subj.right.RSsteps)
    plot(subj.right.RSsteps{i});
end
xlabel('frame'); ylabel('deg/sec');
title(['Subject ',num2str(subj.Num),' Resampled Leading TO-TO cycle'])

subplot(3,1,3)
errorbar(linspace(1,100,subj.right.minVecLen), subj.right.avgStep, subj.right.stdStep);
hold on; grid on;
xlabel('Percent of Gait Cycle'); ylabel('deg/sec');
title(['Subject ',num2str(subj.Num),' Leading foot Angular Velocity Gait Cycle'])

% Figure 2 shows left and right ankle Accel data.%%%%%%%%%%%%%%%%%%%%%%%%%%
fig2 = figure(200+subj.Num);    
set(fig2, 'Units', 'Normalized', 'OuterPosition', [0,0, 1, 1]);
subplot(3,2,1)
plot(subj.right.RSAccSteps{1});
hold on; grid on;
for i=1:length(subj.right.RSAccSteps)
    plot(subj.right.RSAccSteps{i});
end
xlabel('frame'); ylabel('m/s^2');
title(['Subject ',num2str(subj.Num),' Resampled Right Ankle Accels in TO - TO cycle'])

subplot(3,2,2)
plot(subj.left.RSAccSteps{1});
hold on; grid on;
for i=1:length(subj.left.RSAccSteps)
    plot(subj.left.RSAccSteps{i});
end
xlabel('frame'); ylabel('m/s^2');
title(['Subject ',num2str(subj.Num),' Resampled Left Ankle Accels in TO - TO cycle'])

subplot(3,2,3)
errorbar(linspace(1,100,subj.right.minVecLen), subj.right.avgAccStep, subj.right.stdAccStep);
hold on; grid on;
xlabel('Percent of Gait Cycle'); ylabel('m/s^2');
title(['Subject ',num2str(subj.Num),' Right Ankle Accels in TO-TO Cycle'])

subplot(3,2,4)
errorbar(linspace(1,100,subj.left.minVecLen), subj.left.avgAccStep, subj.left.stdAccStep);
hold on; grid on;
xlabel('Percent of Gait Cycle'); ylabel('m/s^2');
title(['Subject ',num2str(subj.Num),' Left Ankle Accels in TO-TO Cycle'])

subplot(3,2,5)
plot(subj.right.avgStepTimes, subj.right.stepVel); grid on;
xlabel('Times(sec)'); ylabel('m/s');
title(['Subject ',num2str(subj.Num),' Mean Right Ankle Velocity in TO-TO Cycle'])

subplot(3,2,6)
plot(subj.left.avgStepTimes, subj.left.stepVel); grid on;
xlabel('Times(sec)'); ylabel('m/s');
title(['Subject ',num2str(subj.Num),' Mean Left Ankle Velocity in TO-TO Cycle'])

% Figure 3 shows waist accels in cycles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig3 = figure(300+subj.Num);    
set(fig3, 'Units', 'Normalized', 'OuterPosition', [0,0, 1, 1]);
subplot(2,2,1)
plot(subj.waist.cycleAccX{1});
hold on; grid on;
for i=1:length(subj.waist.cycleAccX)
    plot(subj.waist.cycleAccX{i});
end
xlabel('frame'); ylabel('m/s^2');
title(['Subject ',num2str(subj.Num),' Waist X Accel W/R to Leading TO-TO cycle'])

subplot(2,2,2)
plot(subj.waist.RScycleAccX{1});
hold on; grid on;
for i=1:length(subj.waist.RScycleAccX)
    plot(subj.waist.RScycleAccX{i});
end
xlabel('frame'); ylabel('m/s^2');
title(['Subject ',num2str(subj.Num),' Resampled Waist X Accel W/R to Leading TO-TO cycle'])

subplot(2,2,3)
errorbar(linspace(1,100,subj.waist.minCycleVecLen), subj.waist.avgCycleAccX, subj.waist.stdCycleAccX);
hold on; grid on;
errorbar(linspace(1,100,subj.waist.minCycleVecLen), subj.waist.avgCycleAccY, subj.waist.stdCycleAccY, 'Color' , 'red');
xlabel('Percent of Gait Cycle'); ylabel('m/s^2'); legend('Horizontal','Vertical')
title(['Subject ',num2str(subj.Num),' Global Waist Accel Gait Cycle'])

subplot(2,2,4)
plot(subj.waist.cycleTimes,subj.waist.cycleVel); grid on;
xlabel('Times(sec)'); ylabel('m/s');
title(['Subject ',num2str(subj.Num),' Mean Waist Horizontal Velocity Cycle'])

%% outputs the data found to excel file and saves figures
% determines folder/worksheet to save data/pictures
if timesSigma == 3
    if yAxis
        section = 'Sigma3XY';
    else
        section = 'Sigma3';
    end
else
    if yAxis
        section = 'Sigma2XY';
    else
        section = 'Sigma2';
    end
end
% saves data and pictures in folder/worksheet
writematrix([subj.Num, mean(subj.waist.avgCycleAccX) mean(subj.waist.avgCycleAccY) subj.waist.cycleAvgVel, subj.right.avgVel, subj.right.avgStepTime, subj.right.stdStepTime, length(subj.right.steps), length(subj.right.RSsteps),subj.left.avgVel, subj.left.avgStepTime, subj.left.stdStepTime, length(subj.left.steps), length(subj.left.RSsteps)],'cycleSeparationData.xlsx','Sheet',section,'Range',['A',num2str(subj.Num+2)])
set(fig1, 'PaperOrientation', 'landscape', 'PaperUnits', 'normalized', 'PaperPosition',[0,0,1,1])
saveas(fig1, ['IMU Data Images/',section,'/subj',num2str(subj.Num),'LeadingAngVel',section,'.pdf']);
set(fig2, 'PaperOrientation', 'landscape', 'PaperUnits', 'normalized', 'PaperPosition',[0,0,1,1])
saveas(fig2, ['IMU Data Images/',section,'/subj',num2str(subj.Num),'AnkleAccels',section,'.pdf']);
set(fig3, 'PaperOrientation', 'landscape', 'PaperUnits', 'normalized', 'PaperPosition',[0,0,1,1])
saveas(fig3, ['IMU Data Images/',section,'/subj',num2str(subj.Num),'WaistAccels',section,'.pdf']);

end
