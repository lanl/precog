function [] = plot_mutLoadSeries()

global times;
global filePath;

fileName = strcat(filePath, 'out.mutationSeries');

Y = load(fileName);
%Y = Y(1:length(times),:);
[rows, cols] = size(Y);
totals = sum(Y,2);
D = zeros(size(Y));
for i = 1:rows
    D(i,:) = Y(i,:) / totals(i);
end
hmo = HeatMap(D', 'ColorMap', cool, 'Symmetric', false);
addXLabel(hmo, 'Time (years)')
addYLabel(hmo, 'Deleterious mutations')

end

