function cycle = gaitCycle(cRA, cLA, cW)
%% IMU Data Resampling
% This function takes in IMU data and calculated angles. Then using the
% steps indexes to resample the data. This is used to normalized the data
% based on heel strike to heel strike. StepsCell should contain the step
% indexes.  Gyro data and angles are optional arguments.
%
% The function assumes the following orientation of the data:
% +x: direction of walking
% +y: medial to Left Side
% +z: normal to ground (vertical)
%
% Function by Erick Nunez

%% Finds leading side
if cRA.steps.indexes{1}(1) < cLA.steps.indexes{1}(1)
    idxs = cRA.steps.indexes; 
    meanFrames = cRA.steps.meanFrames;
else
    idxs = cLA.steps.indexes;
    meanFrames = cLA.steps.meanFrames;
end

%% Pulls Cycle data
for i = 1:length(idxs)
    rA.acc{i} = cRA.globalAcc(idxs{i},:);
    rA.gyro{i} = cRA.gyro(idxs{i},:);
    rA.ang{i} = cRA.angles(idxs{i},:);
    lA.acc{i} = cLA.globalAcc(idxs{i},:);
    lA.gyro{i} = cLA.gyro(idxs{i},:);
    lA.ang{i} = cLA.angles(idxs{i},:);
    w.acc{i} = cW.globalAcc(idxs{i},:);
    w.gyro{i} = cW.gyro(idxs{i},:);
    w.ang{i} = cW.angles(idxs{i},:);
end

%% Resamples 
rsRAgyro = cell(1,length(idxs));
rsRAacc = cell(1,length(idxs));
rsRAang = cell(1,length(idxs));
rsLAgyro = cell(1,length(idxs));
rsLAacc = cell(1,length(idxs));
rsLAang = cell(1,length(idxs));
rsWgyro = cell(1,length(idxs));
rsWacc = cell(1,length(idxs));
rsWang = cell(1,length(idxs));
for i = 1:length(idxs)
    rsRAgyro{i} = resample(rA.gyro{i}, meanFrames, length(idxs{i}));
    rsRAacc{i} = resample(rA.acc{i}, meanFrames, length(idxs{i}));
    rsRAang{i} = resample(rA.ang{i}, meanFrames, length(idxs{i}));
    rsLAgyro{i} = resample(lA.gyro{i}, meanFrames, length(idxs{i}));
    rsLAacc{i} = resample(lA.acc{i}, meanFrames, length(idxs{i}));
    rsLAang{i} = resample(lA.ang{i}, meanFrames, length(idxs{i}));
    rsWgyro{i} = resample(w.gyro{i}, meanFrames, length(idxs{i}));
    rsWacc{i} = resample(w.acc{i}, meanFrames, length(idxs{i}));
    rsWang{i} = resample(w.ang{i}, meanFrames, length(idxs{i}));
end

%% Finds profiles
[rA.meanAcc,    rA.stdAcc]  = findProfile(meanFrames, rsRAacc);
[rA.meanGyro,   rA.stdGyro] = findProfile(meanFrames, rsRAgyro);
[rA.meanAng,    rA.stdAng]  = findProfile(meanFrames, rsRAang);
[lA.meanAcc,    lA.stdAcc]  = findProfile(meanFrames, rsLAacc);
[lA.meanGyro,   lA.stdGyro] = findProfile(meanFrames, rsLAgyro);
[lA.meanAng,    lA.stdAng]  = findProfile(meanFrames, rsLAang);
[w.meanAcc,     w.stdAcc]   = findProfile(meanFrames, rsWacc);
[w.meanGyro,    w.stdGyro]  = findProfile(meanFrames, rsWgyro);
[w.meanAng,     w.stdAng]   = findProfile(meanFrames, rsWang);

%% Finds indexes for gait phase transitions
% Mid Swing
[~,rA.MSwIdx] = min(rA.meanGyro(:,2));
[~,lA.MSwIdx] = min(lA.meanGyro(:,2));
% Toe Off
rA.TOIdx = findToeOff(rA.meanGyro,rA.MSwIdx);
lA.TOIdx = findToeOff(lA.meanGyro,lA.MSwIdx);
% Heel Strike
rA.HSIdx = findHeelStrike(rA.meanGyro, rA.MSwIdx);
lA.HSIdx = findHeelStrike(lA.meanGyro, lA.MSwIdx);

%% Save Data for output
cycle.rA.meanAcc = rA.meanAcc;      cycle.rA.stdAcc = rA.stdAcc;  
cycle.rA.meanGyro = rA.meanGyro;    cycle.rA.stdGyro = rA.stdGyro;
cycle.rA.meanAng = rA.meanAng;      cycle.rA.stdAng = rA.stdAng;  
cycle.lA.meanAcc = lA.meanAcc;      cycle.lA.stdAcc = lA.stdAcc; 
cycle.lA.meanGyro = lA.meanGyro;    cycle.lA.stdGyro = lA.stdGyro;
cycle.lA.meanAng = lA.meanAng;      cycle.lA.stdAng = lA.stdAng;  
cycle.w.meanAcc = w.meanAcc;        cycle.w.stdAcc = w.stdAcc;   
cycle.w.meanGyro = w.meanGyro;      cycle.w.stdGyro = w.stdGyro; 
cycle.w.meanAng = w.meanAng;        cycle.w.stdAng = w.stdAng;
cycle.percent = linspace(0,100,meanFrames)';
cycle.rA.MSwIdx = rA.MSwIdx;        cycle.lA.MSwIdx = lA.MSwIdx;        
cycle.rA.TOIdx = rA.TOIdx;          cycle.lA.TOIdx = lA.TOIdx;
cycle.rA.HSIdx = rA.HSIdx;          cycle.lA.HSIdx = lA.HSIdx;

end

%% Additional Functions
function [meanProfile, stdProfile] = findProfile(frames, dataXYZcell)
for i = 1:frames
    samples = [];
    for j = 1:length(dataXYZcell)
        samples(end+1,:) = dataXYZcell{j}(i,:);
    end
    meanProfile(i,:) = mean(samples);   stdProfile(i,:) = std(samples);
end
end

function TOIdx = findToeOff(meanGyro,MSwIdx)
gyro = [meanGyro; meanGyro];
Idx = MSwIdx + length(meanGyro);
for i = 2:30
    if gyro(Idx-1-i,2) < gyro(Idx-i,2) && gyro(Idx+1-i,2) < gyro(Idx-i,2)
        TOIdx = Idx-i;
        break
    end
end
if TOIdx > length(meanGyro)
    TOIdx = TOIdx - length(meanGyro);
end
end

function HSIdx = findHeelStrike(meanGyro, MSwIdx)
gyro = [meanGyro; meanGyro];
for i = 2:30
    if gyro(MSwIdx-1+i,2) < gyro(MSwIdx+i,2) && gyro(MSwIdx+1+i,2) < gyro(MSwIdx+i,2)
        HSIdx = MSwIdx+i;
        break
    end
end
if HSIdx > length(meanGyro)
    HSIdx = HSIdx - length(meanGyro);
end
end