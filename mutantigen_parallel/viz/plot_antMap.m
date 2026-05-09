function [ output_args ] = plot_antMap()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
global majorTypes;
global cmapTypes;

in = load('out.antigenicDistances');
D = in(majorTypes, majorTypes);
types = length(D(1,:));
[Y,eigvals] = cmdscale(D);

for i = 1:types
    %plot(Y(i,1),Y(i,2),'o', 'MarkerSize', 10.0, 'MarkerEdgeColor', cmap(find(i == order),:), 'MarkerFaceColor', cmap(find(i == order),:));
    plot(Y(i,1),Y(i,2),'o', 'MarkerSize', 10.0, 'MarkerEdgeColor', cmapTypes(i,:), 'MarkerFaceColor', cmapTypes(i,:));
    hold on;
end

% Check goodness of 2d and 3d reconstructions
%maxerr3 = max(abs(D - pdist(Y(:,1:3))))  % Good reconstruction in 3D
%maxerr2 = max(abs(D - pdist(Y(:,1:2))))  % Poor reconstruction in 2D

end

