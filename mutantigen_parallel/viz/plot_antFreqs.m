function [] = plot_antFreqs()

global majorTypes;
global cmapTypes;
global times;

in = load('out.antigenicSamples');
times = in(:,1);
Y = in(:,2:end);
totals = sum(Y,1);
majorTypes = 1:length(Y(1,:)); %find(totals > 100);
cmapTypes = jet(length(majorTypes));

counts = Y(:,majorTypes);
totals = sum(counts,2);
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

freqs = zeros(size(counts));
for n = 1:length(times)
    freqs(n,:) = counts(n,:) / sum(counts(n,:));   
end

h = area(times,freqs); ylim([0 1]); ylabel('Frequency'); xlabel('Time (years)');
for i = 1:length(majorTypes)
    set(h(i),'FaceColor', cmapTypes(i,:))
end
set(gca,'Ydir','reverse');
set(gca,'yticklabel', 1:-0.1:0)
end

