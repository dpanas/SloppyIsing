%-------------------------------%
% function: GetSpikesText
%           Get spikes from a text file; if needed, already bin them into
%           desired time bins (then binsize must be supplied); if needed,
%           already filter channels (then desired criteria must be 
%           supplied). 
%           It assumes the text file to contain two columns, first with
%           channel numbers, second with time frames of recorded spikes.
%           Channel numbering starts from 1 and channels are assumed to be
%           arranged in a 63x63 lattice. 
%
% dependancy: - format of the .txt files
%
% input:  - name of the file (in current path only!!!);
%         - size of the bin to bin the data (if needed, otherwise []) [sec];
%         - lower frequency for filtering out channels  (if needed, otherwise [])[Hz];
%         - higher frequency for filtering out channels (if needed, otherwise [])[Hz];
%         - sampling frequency with which spikes were recorded [Hz];
%         - recording time [sec];
%
% output:  - array holding: channel number;
%                           x coordinate on MEA;
%                           y coordinate on MEA; 
%                           number of spikes;
%                           mean interspike interval [sec];
%                           standard deviation of the mean interspike interval [sec];
%          - a cell array; each cell is a vector with times of spikes [sec]
%          [and if binsize supplied:
%          - array holding: channel number;
%                           x coordinate on MEA;
%                           y coordinate on MEA; 
%                           number of spikes;
%                           number of active bins;
%          - a cell array; each cell is a vector with time bins of spikes]
%
% DAP Apr 2015
% !!! no error control !!!
%-------------------------------%

function [chankey,spikes,chankey_binned,spikes_binned] = GetSpikesText(filename,binsize,lo_thresh,up_thresh,freq,rtime)

disp(' ')
disp(['GetSpikes from: ',' '' ',filename,' '' '])

% reading in the data from a .txt file:
a = load(filename);
chans = a(:,1);
spiketimes = a(:,2);

% arrays / cells for the results:
chankey = zeros(64*64,6);
spikes = cell(64*64,1);
chankey_binned = zeros(64*64,5);      % key binding channel coordinates and statistics
spikes_binned = cell(64*64,1);        % cell array for active channels

i = 0;
j = 0;
% reading in the data:
for y = 1:64
    for x = 1:64
        j = j + 1;     
        prop = find(chans==j);
        if ~isempty(prop)
            spiketrain = double(spiketimes(prop))./double(freq);     % convert to seconds
            if (length(spiketrain)>=1) 
                i = i+1;
                chankey(i,1) = j;               % channel number
                chankey(i,2:3) = [x y];         % channel coordinates on MEA
                spikes{j} = spiketrain;         % spike times [sec]
                chankey(i,4) = length(spikes{j});              % number of spikes
                chankey(i,5) = mean(diff(spikes{j}));          % mean interspike interval
                chankey(i,6) = std(diff(spikes{j}));          % standard deviation of the interspike interval
                if (~isempty(binsize))
                    chankey_binned(i,1) = j;               % channel number
                    chankey_binned(i,2:3) = [x y];         % channel coordinates on MEA
                    spikes_binned{j} = ceil(spiketrain./binsize);         % spike bins
                    chankey_binned(i,4) = length(spikes_binned{j});  
                    prop3 = find(diff(spikes_binned{j})==0);
                    spikes_binned{j}(prop3+1) = [];              % getting rid of the extra spikes in the same bin
                    chankey_binned(i,5) = length(spikes_binned{j});       % number of active bins
                end
            end
        end                   
    end 
end
chankey((chankey(:,1)<1),:) = [];
chankey_binned((chankey_binned(:,1)<1),:) = [];

if (~isempty(binsize))
    chankey_binned((chankey_binned(:,1)<1),:) = [];
    disp(['Spikes binned: ',num2str(binsize),' sec.'])
    % if needed, getting rid of channels with too few active bins:
    if (~isempty(lo_thresh))
        disp(['Channels below ',num2str(lo_thresh),' Hz filtered out.'])
        chankey_binned = chankey_binned(find(chankey_binned(:,5)>=lo_thresh*rtime),:);
    end
    % if needed, getting rid of channels with too many active bins:
    if (~isempty(up_thresh))
        disp(['Channels above ',num2str(up_thresh),' Hz filtered out.'])
        chankey_binned = chankey_binned(find(chankey_binned(:,5)<=up_thresh*rtime),:);
    end  
else
    chankey_binned = [];
    spikes_binned = [];
end

end