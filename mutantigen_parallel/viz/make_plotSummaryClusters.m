function [] = make_plotSummaryClusters()

global filePath;
filePath = './example_data/';
global clusterCutOff;
clusterCutOff = 0.75;

% Plot time series
figure(1)
plot_timeSeries();

% Plot antigenic map colored by clusters
figure(2)
plot_antMapClusters();

% Plot strain freq dynamics
figure(3)
plot_clusterFreqs();

% Plot strain prev dynamics
figure(4)
plot_clusterPrevalence();

% Plot load time series
%plot_mutLoadSeries();

% Plot painted mut tree
mutTreeFile = strcat(filePath, 'out.simmapLoad');
simmapLoad_plot(mutTreeFile);

% Plot painted antigenic tree
cladeTreeFile = strcat(filePath, 'out.simmapAntigenic');
simmapAntigenic_plot(cladeTreeFile);

end

