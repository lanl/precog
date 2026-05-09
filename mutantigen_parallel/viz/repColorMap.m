function [repMap] = repColorMap(total)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

cycleLength = 20;
baseMap = hsv(cycleLength);
repMap = zeros(total, 3);

startLoc = 1;
endLoc = cycleLength;
while (endLoc <= total)
    repMap(startLoc:endLoc,1:3) = baseMap;
    startLoc = startLoc + cycleLength;
    endLoc = endLoc + cycleLength;
end
lastLoc = startLoc;
fillTo = total - lastLoc + 1;
repMap(lastLoc:total,:) = baseMap(1:fillTo,:);




end

