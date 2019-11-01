%% Shimmer Vicon All Trials
% Scripts and functions by erick nunez

%% Variables
trials = [1 2 3 4 5];
methods = [1 2 3 4];

%% Run analysis
for i = 1:length(trials)
   for j = 1:length(methods)
      ShimmerViconTestExampleCode(i,j); 
   end
end