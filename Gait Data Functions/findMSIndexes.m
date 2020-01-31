function steps = findMSIndexes(gryoXYZ, anglesXYZ, times, sigma)
%% Finds MS Indexes
% This function is used to find the Mid Stance point for each step based on 
% the gyroscope data. It find the local maxs first which are the HS and TO
% indexes. Then it find the midpoint between them which ends up being
% midstance.
% +x: direction of walking
% +y: medial to Left Side
% +z: normal to ground (vertical)
% 
% Function by Erick Nunez

%% Variables to be used
indexes = {};

%% Finds Maxs (HS to TO)
HSTOindex = find(islocalmax(gryoXYZ(:,2), 'MinSeparation', 0.4, 'SamplePoints',times));

%% Find MS between HS and TO
if anglesXYZ(HSTOindex(1),2) > anglesXYZ(HSTOindex(2),2)
    HSTOindex(1) = [];
end
if anglesXYZ(HSTOindex(end),2) < anglesXYZ(HSTOindex(end-1),2)
    HSTOindex(end) = [];
end
MSindex = [];
for i = 1:length(HSTOindex)/2
    MSindex = [MSindex round((HSTOindex(i*2)+HSTOindex(i*2-1))/2)]; 
end

%% separates steps
for i =2:length(MSindex)
    indexes{end+1} = MSindex(i-1):MSindex(i)-1;
end

%% Finds mean and std for steps
for i = 1:length(indexes)
    frameLengths(i) = length(indexes{i});
end
meanFrames = mean(frameLengths);    stdFrames = std(frameLengths);

%% Removes outliers in steps
removed = 0;
for i = 1:length(indexes)
    if (length(indexes{i-removed}) < meanFrames - stdFrames*sigma) || (length(indexes{i-removed}) > meanFrames + stdFrames*sigma)
        indexes(i-removed) = [];
        removed = removed + 1;
    end
end

%% Save data for output
steps.removed = removed;            steps.indexes = indexes;
steps.meanFrames = meanFrames;      steps.stdFrames = stdFrames;
steps.MSindexes = MSindex;

end

%% Additional Functions
% function dataOut = calculateInt(rate, dt)
% dataOut = zeros(size(rate));
% for i = 2:length(dataOut)
%     dataOut(i,:) = dataOut(i-1,:) + (rate(i,:) + rate(i-1,:)).*dt/2;
% end
% end

