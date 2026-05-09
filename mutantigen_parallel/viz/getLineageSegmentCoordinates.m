function [lineageSegs] = getLineageSegmentCoordinates(lineageSegs, X, Y)
%   Set the coordinates for the lineage segments
%   Detailed explanation goes here

lines = lineageSegs.numLineages;

lineageSegs.x = cell(lines,1);

%Set x coordinates for each segment in each lineage
for ln = 1:lines
    numSegs = size(lineageSegs.segLengths{ln},2);
    if (numSegs == 1)
        %There is only one segment in this lineage
        lineageSegs.x{ln} = X((1:2),ln);
    else
        %Multiple segments in lineage
        currXVal = X(1,ln);
        for seg = 1: numSegs
            nextXVal = currXVal - lineageSegs.segLengths{ln}(seg);
            lineageSegs.x{ln}(1,seg) = currXVal;
            lineageSegs.x{ln}(2,seg) = nextXVal;
            currXVal = nextXVal;
        end
    end
end

end

