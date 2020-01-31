function waist = invertedPendulum(waist, steps, height)
%% inverted Pendulum
% Function assumed the following orientation of the data:
% +x: direction of walking
% +y: medial to Left Side
% +z: normal to ground (vertical)
%
% Function by erick nunez

%% Variables
legLen = 0.53*height;
HSidx = steps.HSIdx;              TOidx = steps.TOIdx;
MStIdx = steps.MStIdx;
indexes = nan(length(HSidx),3);

%% Creates steps
for i=1:length(HSidx)
    HS = HSidx(i);
    [~, MSti] = min(abs(MStIdx-HS));
    if MStIdx(MSti) < HS && MSti ~= length(MStIdx)
        MSti = MSti + 1;
    end
    [~, TOi] = min(abs(TOidx-MStIdx(MSti)));
    if TOidx(TOi) < MStIdx(MSti) && TOi ~= length(TOidx)
        TOi = TOi + 1;
    end
    if (MSti < length(MStIdx) && TOi < length(TOidx)) && ((HSidx(i)<MStIdx(MSti)) && (MStIdx(MSti)<TOidx(TOi)))
        indexes(i,:) = [HSidx(i), MStIdx(MSti), TOidx(TOi)];
    end
end
indexes = reshape(indexes(~isnan(indexes)),[],3);

%% remove outliers
durStance = indexes(:,3) - indexes(:,1);
toRemove = durStance > steps.meanFrames;
indexes(toRemove,:) = [];
acc = cell(length(indexes),1);        times = cell(length(indexes),1);
vel = cell(length(indexes),1);        pos = cell(length(indexes),1);
vertDisp = nan(length(indexes),1);    stepLen = nan(length(indexes),1);
speed = nan(length(indexes),1); 

%% Pull Accelerations
for i = 1:length(indexes)
    times{i} = waist.times(indexes(i,1):indexes(i,3));
    acc{i} = waist.globalAcc(indexes(i,1):indexes(i,3),:);
end

%% Find Step velocotiy and positions 
for i = 1:length(indexes)
    vel{i} = calculateInt(acc{i}, mean(diff(times{i})));
    pos{i} = calculateInt(vel{i}, mean(diff(times{i})));
end

%% Uses inverted pendulum model to find speed
for i = 1:length(indexes)
    vertDisp(i) = max(pos{i}(:,3))- min(pos{i}(:,3));
    stepLen(i) = 2 * sqrt(2 * legLen * vertDisp(i) - vertDisp(i)^2);
    speed(i) = stepLen(i)/(times{i}(end)-times{i}(1));
end

%% Save outputs
waist.steps.vel = vel;
waist.steps.pos = pos;
waist.steps.vertDisp = vertDisp;
waist.steps.stepLen = stepLen;
waist.steps.speed = speed;

end

%% Additional Functions
function dataOut = calculateInt(rate, dt)
dataOut = zeros(size(rate));
for i = 2:length(dataOut)
    dataOut(i,:) = dataOut(i-1,:) + (rate(i,:) + rate(i-1,:)).*dt/2;
end
end