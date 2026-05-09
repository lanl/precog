function [void] = plot_timeSeries()
%UNTITLED2 Summary of this function goes here

global filePath;

fileName = strcat(filePath, 'timeSeries');

load(fileName);
times = out(:,1);
I = out(:,9);
C = out(:,11);

plot(times, I, 'k', 'LineWidth', 2.0);
ylabel('Prevalence');
xlabel('Time (years)')
%hold on;
%plot(times, C, 'b', 'LineWidth', 2.0);


end

