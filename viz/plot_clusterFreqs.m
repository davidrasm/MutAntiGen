function [void] = plot_clusterFreqs()

global clustMap;
global cmapTypes;
global times;
global filePath;

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
for i = 1:nClusters
    set(h(i),'FaceColor', cmapTypes(i,:))
end
set(gca,'Ydir','reverse');
set(gca,'yticklabel', 1:-0.1:0)

end
