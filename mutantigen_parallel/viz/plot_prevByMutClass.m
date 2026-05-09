function [ output_args ] = plot_prevByMutClass(Z)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

h = bar3(Z);
colorbar

for k = 1:length(h)
    zdata = get(h(k),'ZData');
    set(h(k),'CData',zdata,'FaceColor','interp')
end

end

