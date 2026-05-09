function [ output_args ] = plot_antMapClades()

%global majorTypes;
global cmapTypes;
global filePath;

fileName = strcat(filePath, 'out.antigenicDistances');

in = load(fileName);
D = in; %in(majorTypes, majorTypes);

tipClades = load('out.tipClades');
types = length(D(1,:));

[Y,eigvals] = cmdscale(D);

for i = 1:types
    clade = tipClades(i); 
    plot(Y(i,1),Y(i,2),'o', 'MarkerSize', 10.0, 'MarkerEdgeColor', cmapTypes(clade,:), 'MarkerFaceColor', cmapTypes(clade,:));
    hold on;
end

% Check goodness of 2d and 3d reconstructions
%maxerr3 = max(abs(D - pdist(Y(:,1:3))))  % Good reconstruction in 3D
%maxerr2 = max(abs(D - pdist(Y(:,1:2))))  % Poor reconstruction in 2D

end

