classdef viconC3Ddata < handle
    % viconC3Ddata is a class that imports the data from a C3D file created
    % in Vicon. The class depends on Windows running the C3D server.
    % Functions are pulled from the C3Dserver PDF.
    % Tested in Windows 7 with Matlab R2019a.
    % Values are converted to doubles and strings to make it easier to work
    % with in Matlab.
    %
    % Scripted by erick nunez
    
    % ---------------------------------------------------------------------
    % Properties
    % ---------------------------------------------------------------------
    properties (Access = public)
        labels
        data
        times
        startFrame
        units
        frameRate
        numberOfPoints
        numberOfFrames
    end
    
    properties (Access = private)
        fileName
        sdk
        dataPulled = false;
        labelsIndex
        startFrameIndex
        unitIndex
        rateIndex
        usedIndex
        framesIndex
        startHeader
        endHeader
    end
    
    % ---------------------------------------------------------------------
    % Methods
    % ---------------------------------------------------------------------
    methods
        function obj = viconC3Ddata(fileName)
            if strcmp('PCWIN64',computer)
                addpath('C:\Program Files\Common Files\Motion Lab Systems\C3Dserver')
            end
            
            obj.sdk = actxserver('C3DServer.C3D');
            obj.fileName = fileName;
            
            if ~obj.sdk.Open(fileName,3)
                obj = getIndexFromC3D(obj);
                obj = getDataFromC3D(obj);
            end
            
            obj.times = (0:(obj.numberOfFrames-1))'/obj.frameRate;
        end
        
        function state = dataFound(obj)
            state = obj.dataPulled;
        end
        
        function viewSDKmethods(obj)
            methodsview(obj.sdk)
        end
        
        function out = getIndex(obj,groupName,parameterName)
            out = double(obj.sdk.GetParameterIndex(groupName,parameterName));
        end
        
        function out = getValue(obj, index, itemNumber)
           value = obj.sdk.GetParameterValue(index, itemNumber); 
           if ischar(value)
               out = convertCharsToStrings(value);
           else
               out = double(value);
           end
        end
        
        function out = getPointData(obj, point, dir, startFrame, endFrame)
            out = cell2mat(obj.sdk.GetPointDataEx(point, dir, startFrame, endFrame, '0'));
        end
                
%         function plot3D(obj)
%             
%         end
        
        function delete(obj)
            obj.sdk.Close()
        end
    end
    
end

% -------------------------------------------------------------------------
% Functions
% -------------------------------------------------------------------------
function obj = getDataFromC3D(obj)
obj.startHeader     = obj.sdk.GetVideoFrameHeader(0);
obj.endHeader       = obj.sdk.GetVideoFrameHeader(1);

obj.startFrame      = obj.getValue(obj.startFrameIndex,0);
obj.units           = obj.getValue(obj.unitIndex,0);
obj.frameRate       = obj.getValue(obj.rateIndex,0);
obj.numberOfPoints  = obj.getValue(obj.usedIndex,0);
obj.numberOfFrames  = obj.getValue(obj.framesIndex,0);

obj.labels  = cell(obj.numberOfPoints, 1);
obj.data    = nan(obj.numberOfFrames, 3, obj.numberOfPoints);

for i = 1:obj.numberOfPoints
    x = obj.getPointData(i-1, 0, obj.startHeader, obj.endHeader);
    y = obj.getPointData(i-1, 1, obj.startHeader, obj.endHeader);
    z = obj.getPointData(i-1, 2, obj.startHeader, obj.endHeader);
    obj.data(:,:,i) = [x,y,z];
    obj.labels{i} = obj.getValue(obj.labelsIndex, i-1);
end

if strcmp(obj.units,'mm')
    obj.data = obj.data/1000;
    obj.units = 'm';
end

obj.dataPulled = true;
end

function obj = getIndexFromC3D(obj)
obj.labelsIndex     = obj.getIndex('POINT','LABELS');
obj.startFrameIndex = obj.getIndex('POINT','DATA_START');
obj.unitIndex       = obj.getIndex('POINT','UNITS');
obj.rateIndex       = obj.getIndex('POINT','RATE');
obj.usedIndex       = obj.getIndex('POINT','USED');
obj.framesIndex     = obj.getIndex('POINT','FRAMES');
end
