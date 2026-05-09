function [ output_args ] = plot_timeSeriesMeanLoad( input_args )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

filePath = './FluLoadEvo_noAntigenic_N400_eps012/';


load(strcat(filePath,'timeseries'));
times = out(:,1);
I = out(:,9);
C = out(:,11);

%subplot(2,1,1);
%plot(times, I, 'k', 'LineWidth', 2.0);
%ylabel('Prevalence');
%xlabel('Time (years)')

Y = load(strcat(filePath,'out.mutationSeries'));
Y = Y(1:length(times),:);
[rows, cols] = size(Y);
totals = sum(Y,2);
D = zeros(size(Y));
means = zeros(rows,1);
for i = 1:rows
    D(i,:) = Y(i,:) / totals(i);
    means(i) = sum(D(i,:) * (1:1:50)');
end
%subplot(2,1,2);
%plot(times, means, 'r', 'LineWidth', 2.0);

textSize = 14;
[AX, H1, H2] = plotyy(times, I, times, means);
%load('timeseriesNoLoad');
%times = out(:,1);
%I = out(:,9);
%hold on; plot(times, I, 'LineWidth', 2.0, 'Color', [0.4 0.4 0.4], 'LineStyle', ':');
set(AX,{'ycolor'},{'k';'r'})
set(H1, 'LineWidth', 2.0, 'Color', 'k');
set(H2, 'LineWidth', 2.0, 'Color', 'r');
axes(AX(1));
ylabel('Prevalence', 'FontSize', textSize)
axes(AX(2));
ylabel('Mean load', 'FontSize', textSize)
box off;


end