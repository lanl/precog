function [ output_args ] = plot_antMapClusters()

%global majorTypes;
global cmapTypes;
global filePath;
global clusterCutOff;
global clustMap;

fileName = strcat(filePath, 'out.antigenicDistances');

in = load(fileName);
D = in; %in(majorTypes, majorTypes);

%tipClades = load('out.tipClades');
types = length(D(1,:));

% Cluster by distance
links = linkage(D);
%dendrogram(links);
clustMap = cluster(links, 'cutoff', clusterCutOff);
%clustMap = cluster(links,'maxclust',3);

nClusters = max(clustMap);
if (nClusters > 20)
    cmapTypes = repColorMap(nClusters);
else
    cmapTypes = hsv(nClusters);
end

% Project down to 2D
[Y,eigvals] = cmdscale(D);


for i = 1:types
    clust = clustMap(i); 
    plot(Y(i,1),Y(i,2),'o', 'MarkerSize', 10.0, 'MarkerEdgeColor', cmapTypes(clust,:), 'MarkerFaceColor', cmapTypes(clust,:));
    hold on;
end

% Check goodness of 2d and 3d reconstructions
%maxerr3 = max(abs(D - pdist(Y(:,1:3))))  % Good reconstruction in 3D
%maxerr2 = max(abs(D - pdist(Y(:,1:2))))  % Poor reconstruction in 2D

end


