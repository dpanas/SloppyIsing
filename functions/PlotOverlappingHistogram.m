%-------------------------------%
% function: PlotOverlappingHistogram
%           Plots a histogram of provided values in provided bins, using
%           only a line to outline the bars, with provided color (if no
%           color then default is black), and no borders between bars.
%
% dependancy: ---
%
% input:    - values of histogram bins;
%           - bin edges (so needs to be one element longer than values!);
%           - color for the line outlining histrogram shape; 
% 
%  DAP Nov 2014
% !!! no error control !!!
%-------------------------------%

function h = PlotOverlappingHistogram(values,bins,col)

bsize = (bins(2)-bins(1));

if (isempty(col))
    col = [0 0 0];
end

hold on
% first bin separately, because has to start form 0:
h = plot([bins(1) bins(1)],[0 values(1)],'lineWidth',2,'color',col);     % creates a handle to use for legend labelling
plot([bins(1) bins(2)],[values(1) values(1)],'lineWidth',2,'color',col)
% consecutive bins:
for i=2:length(bins)-1
    plot([bins(i) bins(i)],[values(i-1) values(i)],'lineWidth',2,'color',col)
    plot([bins(i) bins(i+1)],[values(i) values(i)],'lineWidth',2,'color',col)
end
% and finally the last bin going back to zero:
plot([bins(end) bins(end)],[values(end) 0],'lineWidth',2,'color',col)

end