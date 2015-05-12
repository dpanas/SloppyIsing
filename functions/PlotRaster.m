%-----------------------------------------------------------------------------------%
% function: PlotRaster
%           plots a raster of spiking activity
%
% dependancy: ---
%
% input:   - channel key: [active channel no; x coordinate on MEA; y
%          coordinate on MEA; number of spikes in channel];
%          - a cell array; each cell is a vector with spike times;
%          - [optional] size of time bin if spike count in bins is to be 
%            plotted instead (in sec);
%
%           DAP Feb 2015
%-----------------------------------------------------------------------------------%

function PlotRaster(chankey,spikes,bin_for_counts)

% sort the channels according to mean rate:
chky = chankey;
chky = sortrows(chky,[4]);

spike_count = [];              % spike count for the histogram
isis = [];                     % interspike intervals from each channel pooled together
H = [];

% go through sorted channels and plot raster if wanted and collect the
% spike count:
for i=1:length(chky)
    if ~exist('bin_for_counts')
        hold on
        plot(spikes{chky(i,1)},i*ones(length(spikes{chky(i,1)}),1),'.','color',[.2 .2 .2],'MarkerSize',3)
    end
    spike_count = [spike_count spikes{chky(i,1)}'];
    isis = [isis diff(spikes{chky(i,1)})'];
end
% if not raster, than spike count in time bins:
if exist('bin_for_counts')
    edges=0:bin_for_counts:ceil(max(spike_count));
    H = histc(spike_count,edges);
    bar(edges(1:end-1)+(bin_for_counts./2),H(1:end-1),'k')
end

end