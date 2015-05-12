%-------------------------------%
% function: PlotArrayLayout
%           Plots the desired values on the layout of the MEA array, e.g.
%           mean firing rate of channel, mean correlation etc.
%
% dependancy: ---
%
% input:    - channel key: [active channel no; x coordinate on MEA; y
%          coordinate on MEA];
%          - a vector with the value to be plotted;
%          - [optional] array size;
% 
%  DAP Feb 2014
% !!! no error control !!!
%-------------------------------%

function PlotArrayLayout(chankey,value,array_size)

if ~exist('array_size')
    MeanAct0 = NaN(65,65);
else 
    MeanAct0 = NaN(array_size+1,array_size+1);
end

for i=1:size(chankey,1)
    MeanAct0(chankey(i,2),chankey(i,3)) = value(i);
end
cscale = [min(min(MeanAct0)) max(max(MeanAct0))];
pcolor(MeanAct0')
shading flat
axis square
colormap jet
caxis(cscale)
    
end