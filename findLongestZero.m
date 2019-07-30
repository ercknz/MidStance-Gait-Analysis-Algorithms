function out = findLongestZero(A, B)
if length(A) < length(B)
    short = A; long = B; disp('b is longer')
else
    short = B; long = A; disp('a is longer')
end
lenDiff = length(long)-length(short);
short = [short zeros(1,lenDiff)];

for i = 1:lenDiff + 1
    cla
    disp(i)
    short = [zeros(1) short];
    short(end) = [];
    resultDiff = floor(long-short);
    %     for j = 2:length(resultDiff)
    %         if resultDiff(j) == 0
    %             if resultDiff(j) == resultDiff(j-1)
    %
    %     end
    plot(long,'b'); grid on; hold on; plot(short,'r')
    plot(resultDiff,'g');
    drawnow
    
end

out = resultDiff;
end