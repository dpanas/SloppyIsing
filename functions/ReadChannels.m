%-------------------------------%
% function: ReadChannels
%           A function to read in data from specified channels in the form 
%           of an array of -1 and 1 representing spikes; columns correspond
%           to channels.
%
% dependancy: ---
%
% input:   - a cell array; each cell is a vector with indices of active
%          bins of a given channel;
%          - number of bins total;
%          - matrix with coordinates of chosen channels (no. x coord. y coord.);
%
% output:  - matrix with binned data [-1 for silent, 1 for spiking]; column
%            = channel
%
% DAP April 2013
% !!! almost no error control !!!
%-------------------------------%

function [data] = ReadChannels(spikes_binned,Nbin,cord)

data = [];

for i=1:length(cord)
   spikes = -1*ones(1,Nbin);
   spikes(spikes_binned{cord(i,1)}) = 1;
   data = [data; spikes];
end

data = data';

end

