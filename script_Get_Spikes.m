% %-------------------------------%
% script:   Get_Spikes
%           This is a script to format spike data for the use with 
%           SloppyIsing package and to compute correlations between
%           channnels. 
%           It assumes the spikes to be in a text format, with first column
%           corresponding to channel number, and second column
%           corresponding to the frame in which a spike occurred. Channel
%           layout is assumed to be a 64x64 matrix, with 4096 channels,
%           numbering starting from 1. Spike occurrences assumed to be integers,
%           i.e. a recording frame rather than time.
%           Time of the recording and sampling frequency are needed to be
%           manually input along with the file name.
%           Plots array laout of activity and raster plot.
%
% !!! This is a script, not a function - parameters need to be changed
% manually by the user upon each call !!!
%
% dependency: GetSpikesText, GetCorrelations, PlotArrayLayout, PlotRaster
%
% DAP Apr 2015
% %-------------------------------%

% first clear the workspace:
clear

% variables are here:
%-------------------%
% bin size and filtering if wanted:
binsize = 0.005;    % desired bin size for binning spikes and computing correlations
freq_lobound = 0.1;    % lower bound frequency for filtering out channels
freq_upbound = 10;    % upper bound frequency for filtering out channels
% sampling frequency of the recording [in Hertz]:
freq = 7022;
% duration of the recording [in seconds]:
rtime = 900;
% names of files of spike data, assumed to be simple .txt files:
% assumed to be in the current path
filenames{1} = './test_text_spikes.txt';
% and the names for files to be written out in a spike format compatible
% with the rest of the package:
fileouts{1} = './results_spikes/test_spikes_from_text.mat';
%-------------------%

% (below no more variables, just code)
%-------------------%
close all
addpath('./functions/')
NOF = length(filenames);
% going through data files on the list:
for i = 1:NOF
    % getting spikes from the data file (produced by Olivers algorithm):
    filename_data = filenames{i};
    [chankey,spikes,chankey_binned,spikes_binned] = GetSpikesText(filename_data,binsize,freq_lobound,freq_upbound,freq,rtime);
    % module for computing correlations:
    if (~isempty(chankey_binned))
        tstart = tic;
        [correlations,distances,correlations_probtest]=GetCorrelations(chankey_binned,spikes_binned,ceil(rtime./binsize));
        time_compute_corr=toc(tstart);
        disp(['GetCorrelations took: ',num2str(time_compute_corr/60),' minutes.'])
        % and average correlation of a given channel with all other avail:
        chankey_binned(1,6) = mean([correlations(1,2:end)]);
        for j = 2:length(chankey_binned)-1
            chankey_binned(j,6) = mean([correlations(1:j-1,j)' correlations(j,j+1:end)]);
        end
        chankey_binned(end,6) = mean([correlations([end-1:-1:1],end)']);
    end
    % file info and date of creation:
    Info_spikes = ['This file was created by script Get_Spikes on ',date];
    % saving results:
    filename_spikes = fileouts{i};
    if (~isempty(chankey_binned))
        save(filename_spikes,'freq','rtime','binsize','freq_lobound','freq_upbound','filename_data','filename_spikes','chankey','spikes','chankey_binned','spikes_binned','correlations','distances','correlations_probtest','time_compute_corr','Info_spikes')
    else
        save(filename_spikes,'freq','rtime','filename_data','filename_spikes','chankey','spikes','Info_spikes')
    end
    disp(' ')
    disp('Script finished, results written succesfully:')
    disp(filename_spikes)
    disp(' ')
    % and plot the results:
    figure
    set(gcf,'position',[100 400 1300 300])
    subplot('position',[.06 .1 .15 .8])
    PlotArrayLayout(chankey,chankey(:,4),64)
    subplot('position',[0.3 .15 .65 .7])
    PlotRaster(chankey,spikes)
end
%---------END----------%