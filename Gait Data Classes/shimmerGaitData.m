classdef shimmerGaitData < handle
   % This object is used to keep track of the subject data collected using 
   % the Shimmer IMU's during gait. The methods help encapsulate the 
   % methods for all the subclasses used at each sensor location during a
   % trial. 
   % 
   % class by erick nunez
   
   % ----------------------------------------------------------------------
   % Properties 
   % ----------------------------------------------------------------------
   
   properties
       subjectNumber
       
       rightAnkle
       leftAnkle
       waist
       leftWrist
       rightWrist
   end
   
   properties (Access = private)
       syncJumpUsed = false;
       
       preUnixTime
       preEndUnixTime
       startUnixTime
       endUnixTime
   end
   
   % ----------------------------------------------------------------------
   % Methods
   % ----------------------------------------------------------------------
   
   methods
       function obj = shimmerGaitData(subjNum, logFilePath, syncJumpUsed)
           logFile = readmatrix(logFilePath);
           obj.subjectNumber = subjNum;
           obj.preUnixTime = logFile(subjNum,2);
           if syncJumpUsed
               obj.syncJumpUsed = syncJumpUsed;
               obj.preEndUnixTime = logFile(subjNum,3);
               obj.startUnixTime = logFile(subjNum,4);
               obj.endUnixTime = logFile(subjNum,5);
           else
               obj.startUnixTime = logFile(subjNum,3);
               obj.endUnixTime = logFile(subjNum,4);
           end
       end
       
       function obj = addSensor(obj, sensorName, fileName)
           switch sensorName
               case 'rightAnkle'
                   if obj.syncJumpUsed
                       obj.rightAnkle = sensorLocation(fileName, obj.preUnixTime, obj.startUnixTime, obj.endUnixTime, obj.preEndUnixTime);
                   else
                       obj.rightAnkle = sensorLocation(fileName, obj.preUnixTime, obj.startUnixTime, obj.endUnixTime);
                   end
               case 'leftAnkle'
                   if obj.syncJumpUsed
                       obj.leftAnkle = sensorLocation(fileName, obj.preUnixTime, obj.startUnixTime, obj.endUnixTime, obj.preEndUnixTime);
                   else
                       obj.leftAnkle = sensorLocation(fileName, obj.preUnixTime, obj.startUnixTime, obj.endUnixTime);
                   end
               case 'waist'
                   if obj.syncJumpUsed
                       obj.waist = sensorLocation(fileName, obj.preUnixTime, obj.startUnixTime, obj.endUnixTime, obj.preEndUnixTime);
                   else
                       obj.waist = sensorLocation(fileName, obj.preUnixTime, obj.startUnixTime, obj.endUnixTime);
                   end
               case 'leftWrist'
                   if obj.syncJumpUsed
                       obj.leftWrist = sensorLocation(fileName, obj.preUnixTime, obj.startUnixTime, obj.endUnixTime, obj.preEndUnixTime);
                   else
                       obj.leftWrist = sensorLocation(fileName, obj.preUnixTime, obj.startUnixTime, obj.endUnixTime);
                   end
               case 'rightWrist'
                   if obj.syncJumpUsed
                       obj.rightWrist = sensorLocation(fileName, obj.preUnixTime, obj.startUnixTime, obj.endUnixTime, obj.preEndUnixTime);
                   else
                       obj.rightWrist = sensorLocation(fileName, obj.preUnixTime, obj.startUnixTime, obj.endUnixTime);
                   end
               otherwise
                   error('invalid sensor location');
           end
       end
       
       function obj = getAcc(obj,range)
           switch nargin
               case 1
                   if ~isempty(obj.rightAnkle)
                       obj.rightAnkle.getAcc();
                   end
                   if ~isempty(obj.leftAnkle)
                       obj.leftAnkle.getAcc();
                   end
                   if ~isempty(obj.rightWrist)
                       obj.rightWrist.getAcc();
                   end
                   if ~isempty(obj.leftWrist)
                       obj.leftWrist.getAcc();
                   end
                   if ~isempty(obj.waist)
                       obj.waist.getAcc();
                   end
               case 2
                   if ~isempty(obj.rightAnkle)
                       obj.rightAnkle.getAcc(range);
                   end
                   if ~isempty(obj.leftAnkle)
                       obj.leftAnkle.getAcc(range);
                   end
                   if ~isempty(obj.rightWrist)
                       obj.rightWrist.getAcc(range);
                   end
                   if ~isempty(obj.leftWrist)
                       obj.leftWrist.getAcc(range);
                   end
                   if ~isempty(obj.waist)
                       obj.waist.getAcc(range);
                   end
               otherwise
                   error('invalid number of arguemnts');
           end
       end
       
       function obj = getGyro(obj,range)
           switch nargin
               case 1
                   if ~isempty(obj.rightAnkle)
                       obj.rightAnkle.getGyro();
                   end
                   if ~isempty(obj.leftAnkle)
                       obj.leftAnkle.getGyro();
                   end
                   if ~isempty(obj.rightWrist)
                       obj.rightWrist.getGyro();
                   end
                   if ~isempty(obj.leftWrist)
                       obj.leftWrist.getGyro();
                   end
                   if ~isempty(obj.waist)
                       obj.waist.getGyro();
                   end
               case 2
                   if ~isempty(obj.rightAnkle)
                       obj.rightAnkle.getGyro(range);
                   end
                   if ~isempty(obj.leftAnkle)
                       obj.leftAnkle.getGyro(range);
                   end
                   if ~isempty(obj.rightWrist)
                       obj.rightWrist.getGyro(range);
                   end
                   if ~isempty(obj.leftWrist)
                       obj.leftWrist.getGyro(range);
                   end
                   if ~isempty(obj.waist)
                       obj.waist.getGyro(range);
                   end
               otherwise
                   error('invalid number of arguemnts');
           end
       end
       
       function obj = getMag(obj,range)
           switch nargin
               case 1
                   if ~isempty(obj.rightAnkle)
                       obj.rightAnkle.getMag();
                   end
                   if ~isempty(obj.leftAnkle)
                       obj.leftAnkle.getMag();
                   end
                   if ~isempty(obj.rightWrist)
                       obj.rightWrist.getMag();
                   end
                   if ~isempty(obj.leftWrist)
                       obj.leftWrist.getMag();
                   end
                   if ~isempty(obj.waist)
                       obj.waist.getMag();
                   end
               case 2
                   if ~isempty(obj.rightAnkle)
                       obj.rightAnkle.getMag(range);
                   end
                   if ~isempty(obj.leftAnkle)
                       obj.leftAnkle.getMag(range);
                   end
                   if ~isempty(obj.rightWrist)
                       obj.rightWrist.getMag(range);
                   end
                   if ~isempty(obj.leftWrist)
                       obj.leftWrist.getMag(range);
                   end
                   if ~isempty(obj.waist)
                       obj.waist.getMag(range);
                   end
               otherwise
                   error('invalid number of arguemnts');
           end
       end
       
       function obj = getSamplingFreq(obj)
          if ~isempty(obj.rightAnkle)
              obj.rightAnkle.getSampFreq();
          end
          if ~isempty(obj.leftAnkle)
              obj.leftAnkle.getSampFreq();
          end
          if ~isempty(obj.rightWrist)
              obj.rightWrist.getSampFreq();
          end
          if ~isempty(obj.leftWrist)
              obj.leftWrist.getSampFreq();
          end
          if ~isempty(obj.waist)
              obj.waist.getSampFreq();
          end
       end
       
       function obj = createFilter(obj,cutOffFreq, filterOrder)
          if ~isempty(obj.rightAnkle)
              obj.rightAnkle.createFilter(cutOffFreq, filterOrder);
          end 
          if ~isempty(obj.leftAnkle)
              obj.leftAnkle.createFilter(cutOffFreq, filterOrder);
          end 
          if ~isempty(obj.rightWrist)
              obj.rightWrist.createFilter(cutOffFreq, filterOrder);
          end 
          if ~isempty(obj.leftWrist)
              obj.leftWrist.createFilter(cutOffFreq, filterOrder);
          end 
          if ~isempty(obj.waist)
              obj.waist.createFilter(cutOffFreq, filterOrder);
          end 
       end
       
       function obj = filterData(obj)
          if ~isempty(obj.rightAnkle)
              obj.rightAnkle.filterData();
          end
          if ~isempty(obj.leftAnkle)
              obj.leftAnkle.filterData();
          end 
          if ~isempty(obj.rightWrist)
              obj.rightWrist.filterData();
          end 
          if ~isempty(obj.leftWrist)
              obj.leftWrist.filterData();
          end 
          if ~isempty(obj.waist)
              obj.waist.filterData();
          end 
       end
   end 
end