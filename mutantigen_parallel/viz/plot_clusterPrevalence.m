function [void] = plot_clusterPrevalence()

global clustMap;
global cmapTypes;
global times;
global filePath;

fileName = strcat(filePath, 'timeSeries');
load(fileName);
%times = out(:,1);
I = out(:,9);

fileName = strcat(filePath, 'out.antigenicSamples');

in = load(fileName);
times = in(:,1);

types = length(in(1,:)) - 2;
Y = in(:,2:(1+types));

nClusters = max(clustMap);
counts = zeros(length(times), nClusters);
for i = 1:types
   clust = clustMap(i);
   counts(:,clust) = counts(:,clust) + Y(:,i);
end

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

h = area(times,prev); ylabel('Prevalence'); xlabel('Time (years)');
for i = 1:nClusters
    set(h(i),'FaceColor', cmapTypes(i,:))
end

end

