function [ output_args ] = plot_attackRates()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%load('./FluLoadEvo_noAntigenic_N400_eps012/timeseries');
load('timeseries');
times = out(:,1);
I = out(:,9);
C = out(:,11);

T = length(times);
N = 40e06;
cycle = 36;
startLoc = 1;
endLoc = cycle;
attackRates = [];
cntr = 1;
while (endLoc < T)
   perc = sum(C(startLoc:endLoc)) / N;
   attackRates(cntr) = perc;
   startLoc = startLoc + cycle;
   endLoc = endLoc + cycle;
   cntr = cntr+1;
end

plot(attackRates, 'k', 'LineWidth', 2.0);
ylabel('Rate');
xlabel('Time (years)')
%hold on;
%plot(times, C, 'b', 'LineWidth', 2.0);


end
