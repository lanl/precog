function [void] = plot_cladePrevalence()

global clades;
global cmapTypes;
global times;
global filePath;

fileName = strcat(filePath, 'timeSeries');

load(fileName);
%times = out(:,1);
I = out(:,9);

fileName = strcat(filePath, 'out.cladeSamples');
in = load(fileName);
times = in(:,1);
Y = in(:,2:end);
clades = 1:length(Y(1,:));
%cmapTypes = jet(length(clades));
cmapTypes = repColorMap(length(clades));

counts = Y(:,clades);
totals = sum(counts,2);

nonZeroTotals = find(totals > 0);
startLoc = min(nonZeroTotals);

times = times(startLoc:end);
counts = counts(startLoc:end,:);
totals = totals(startLoc:end);
I = I(startLoc:end);

noSampleLocs = find(totals == 0);
for i = 1:length(noSampleLocs)
    prevLoc = noSampleLocs(i) - 1;
    prevLocSamples = totals(prevLoc);
    while(prevLocSamples == 0 && prevLoc > 1)
        prevLoc = prevLoc - 1;
        prevLocSamples = totals(prevLoc);
    end
    if(prevLocSamples == 0)
        display('WARNING: Could not interpolate between times with no samples!')
    end
    counts(noSampleLocs(i),:) = counts(prevLoc,:);
end

prev = zeros(size(counts));
for n = 1:length(times)
    prev(n,:) = I(n) * (counts(n,:) / sum(counts(n,:)));   
end

%prev = fliplr(prev);
h = area(times,prev); ylabel('Prevalence'); xlabel('Time (years)');
for i = 1:length(clades)
    set(h(i),'FaceColor', cmapTypes(i,:))
end
%set(gca,'Ydir','reverse');
%set(gca,'yticklabel', 1:-0.1:0)

end

