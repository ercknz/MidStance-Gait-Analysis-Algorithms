classdef sensorLocation < handle
   % sensorLocation is a subclass for the shimmerGaitData class. It stores
   % the local sensor data, filtered data and can rotate the data collected.
   %
   % class by erick nunez

   % ----------------------------------------------------------------------
   % Properties
   % ----------------------------------------------------------------------
   properties
       sampleFreq
       period
       % [x,y,z] for all data
       preAcc
       preGyro
       preMag
       fPreAcc
       fPreGyro
       fPreMag

       acc
       gyro
       mag
       fAcc
       fGyro
       fMag

       preTimes     %secs
       times    %secs
   end

   properties (Access = private)
      data
      unixTimes

      startIndex
      preEndIndex

      accRange = 2:4;
      gyroRange = 5:7;
      magRange = 9:11;

      filterOrder
      cutOffFreq
      filtA
      filtB
   end

   % ----------------------------------------------------------------------
   % Methods
   % ----------------------------------------------------------------------
   methods
       function obj = sensorLocation(dataFileName, preTime, startTime, endTime, preEndTime)
           rawData = readmatrix(dataFileName);
           [~,col] = size(rawData);
           [~, preIndex] = min(abs(rawData(:,1)- preTime));
           [~, endIndex] = min(abs(rawData(:,1)- endTime));
           for i = 1:col
               dataCol = rawData(~isnan(rawData(:,i)),i);
               if ~isempty(dataCol)
                   obj.data(:,i) = dataCol(preIndex:endIndex);
               end
           end
           [~, obj.startIndex] = min(abs(obj.data(:,1)- startTime));
           switch nargin
               case 4
                   obj.preEndIndex = obj.startIndex-1;
                   preEndTime = rawData(obj.preEndIndex,1);
               case 5
                   [~, obj.preEndIndex] = min(abs(obj.data(:,1)-preEndTime));
               otherwise
                   error('invalid number of arguments');
           end
           obj.unixTimes = obj.data(:,1);
           obj.preTimes = (obj.unixTimes(1:obj.preEndIndex)-preEndTime)/1000;
           obj.times = (obj.unixTimes(obj.startIndex:end)-startTime)/1000;
       end

       function out = getRaw(obj)
           out = obj.data;
       end

       function obj = getAcc(obj,range)
           if nargin == 2
               obj.accRange = range;
           elseif nargin ~= 1 && nargin ~= 2
               error('invalid number of arguemnts');
           end
           obj.preAcc = obj.data(1:obj.preEndIndex, obj.accRange);
           obj.acc = obj.data(obj.startIndex:end, obj.accRange);
       end

       function obj = getGyro(obj,range)
           if nargin == 2
               obj.gyroRange = range;
           elseif nargin ~= 1 && nargin ~= 2
               error('invalid number of arguemnts');
           end
           obj.preGyro = obj.data(1:obj.preEndIndex, obj.gyroRange);
           obj.gyro = obj.data(obj.startIndex:end, obj.gyroRange);
       end

       function obj = getMag(obj,range)
           if nargin == 2
               obj.magRange = range;
           elseif nargin ~= 1 && nargin ~= 2
               error('invalid number of arguemnts');
           end
           obj.preMag = obj.data(1:obj.preEndIndex, obj.magRange);
           obj.mag = obj.data(obj.startIndex:end, obj.magRange);
       end

       function obj = getSampFreq(obj)
           obj.period = mean(diff(obj.times));
           obj.sampleFreq = 1/obj.period;
       end

       function obj = createFilter(obj, cutOffFreq, filterOrder)
           obj.filterOrder = filterOrder;
           obj.cutOffFreq = cutOffFreq;
           if ~isempty(obj.sampleFreq)
               [obj.filtB, obj.filtA] = butter(obj.filterOrder, obj.cutOffFreq/(obj.sampleFreq/2));
           else
              error('find sampling freq first')
           end
       end

       function obj = filterData(obj)
           if ~isempty(obj.filtA) && ~isempty(obj.filtB)
               for i = 1:3
                   if ~isempty(obj.acc)
                       obj.fPreAcc(:,i) = filtfilt(obj.filtB, obj.filtA, obj.preAcc(:,i));
                       obj.fAcc(:,i) = filtfilt(obj.filtB, obj.filtA, obj.acc(:,i));
                   end
                   if ~isempty(obj.gyro)
                       obj.fPreGyro(:,i) = filtfilt(obj.filtB, obj.filtA, obj.preGyro(:,i));
                       obj.fGyro(:,i) = filtfilt(obj.filtB, obj.filtA, obj.gyro(:,i));
                   end
                   if ~isempty(obj.mag)
                       obj.fPreMag(:,i) = filtfilt(obj.filtB, obj.filtA, obj.preMag(:,i));
                       obj.fMag(:,i) = filtfilt(obj.filtB, obj.filtA, obj.mag(:,i));
                   end
               end
           else
               error('create filter first')
           end
       end

       function obj = rotateIMU(obj, axis, degrees)
           if ~isempty(obj.acc)
               obj.preAcc = rotate(obj.preAcc, axis, degrees);
               obj.acc = rotate(obj.acc, axis, degrees);
           end
           if ~isempty(obj.gyro)
               obj.preGyro = rotate(obj.preGyro, axis, degrees);
               obj.gyro = rotate(obj.gyro, axis, degrees);
           end
           if ~isempty(obj.mag)
               obj.preMag = rotate(obj.preMag, axis, degrees);
               obj.mag = rotate(obj.mag, axis, degrees);
           end
           if ~isempty(obj.fAcc)
               obj.fPreAcc = rotate(obj.fPreAcc, axis, degrees);
               obj.fAcc = rotate(obj.fAcc, axis, degrees);
           end
           if ~isempty(obj.fGyro)
               obj.fPreGyro = rotate(obj.fPreGyro, axis, degrees);
               obj.fGyro = rotate(obj.fGyro, axis, degrees);
           end
           if ~isempty(obj.fMag)
               obj.fPreMag = rotate(obj.fPreMag, axis, degrees);
               obj.fMag = rotate(obj.fMag, axis, degrees);
           end
       end
   end
end

% -------------------------------------------------------------------------
% Functions
% -------------------------------------------------------------------------
function xyzOUT = rotate(xyzIN, axis, theta)
[length, ~] = size(xyzIN);
xyzOUT = nan(length,3);
switch axis
    case 'x'
        rotMat = [1 0 0; 0 cosd(theta) sind(theta); 0 -sind(theta) cosd(theta)];
    case 'y'
        rotMat = [cosd(theta) 0 -sind(theta); 0 1 0; sind(theta) 0 cosd(theta)];
    case 'z'
        rotMat = [cosd(theta) sind(theta) 0; -sind(theta) cosd(theta) 0; 0 0 1];
    otherwise
        error('invalid axis, use "x", "y", or "z"')
end
for i = 1:length
    xyzOUT(i,:) = (rotMat * xyzIN(i,:)')';
end
end
