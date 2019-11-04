# Rotation-gaitAnalysis

This code takes in the data collected from subjects during medical evaluations. 
Subjects wore IMU's on the waist, wrists and ankles.
Subjects were asked to perfrom different motor taskes inclduing 2 min walking. 
Subjects walked in open indoor enviroments. 
The code pulls the gait data from a log file using unix timestamps. 
Then it segments the data (heel-strike to heel-strike) or (Mid-stance to mid-stance).
Based on averaged step acceleration profiles, we determine speeds. 

Code was update into functions and classes for ease of use. 
Refer to ShimmerViconTestExample to see how functions can be used. 
Classes and functions were also added for processing vicon data within a C3D file.

All code was tested in Matlab 2019a on Windows 7. 

Code is a working progress. 

No subject data or results is pushed to github and only remains locally. 

Scripts and functions by erick nunez
