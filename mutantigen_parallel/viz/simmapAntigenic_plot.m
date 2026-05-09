function [void] = simmapAntigenic_plot(filename)
%Plot SimMap tree in "filename" with lineages colored.
%   Detailed explanation goes here

tr = phytree_read(filename); %phytree_read has been modified to parse lineage annotations in SimMap format
phytree_plot_antigenic(tr); %phytree_plot has been modified to plot lineage line segments


end

