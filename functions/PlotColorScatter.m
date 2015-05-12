%-------------------------------%
% function: PlotColorScatter
%           Plots two vectors as a scatterplot of second vs first, and uses
%           a third vector as a colormap for the points. Default it to use
%           the third vector whole - and then color just denotes the rank
%           of each point within the third vector. If an extra entry is
%           provided, it is assumed to be the edges for histogram binning
%           the third vector - and then the color denotes the bin of each
%           point from the third vector.
% 
% dependancy: ---
%
% input:   - a vector of arguments;
%          - a vector of values;
%          - a vector of associated values for the color scale;
%          - a vector of edges for a histogram for the third vector; [opt]
%          - a vector of colors for the histogram bins. [opt]
% 
%  DAP Aug 2014
% !!! no error control !!!
%-------------------------------%

function PlotColorScatter(vec1,vec2,vec3,edges,col)

hold on

if nargin==3
    col = jet(length(vec3));
    [~,ord] = sort(vec3);
    vec1 = vec1(ord);
    vec2 = vec2(ord);
    for i=1:length(vec1)
        plot(vec1(i),vec2(i),'.','color',col(i,:),'MarkerSize',5)
    end
elseif nargin==4
    col = jet(length(edges));
    [~,bin] = histc(vec3,edges);
    for i=1:length(vec1)
        plot(vec1(i),vec2(i),'.','color',col(bin(i),:),'MarkerSize',5)
    end
elseif nargin==5
    [~,bin] = histc(vec3,edges);
    for i=1:length(vec1)
        plot(vec1(i),vec2(i),'.','color',col(bin(i),:),'MarkerSize',5)
    end
end

end