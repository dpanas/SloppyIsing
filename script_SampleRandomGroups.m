% %-------------------------------%
% script:   SampleRandomGroups
%           A script to choose sample subsets from available channels, for
%           the Ising model fits; it automatically generates the file with
%           'key' of groups necessary to run script_Ising_Main.
%           It randomly chooses 'non' neurons for 'nos' groups and displays
%           on a MEA layout how many times each channel is chosen, for the
%           user to decide if the population sampling is appropriate. The 
%           script uses randperm and removes overrepresented channels after 
%           each iteration. An additional 'max_times' variable can be used 
%           to control how many times a channel can be chosen.
%           (different combinations of these variables might need to be
%           tested depending on how even the sampling is supposed to be)
%
% !!! This is a script, not a function - parameters need to be changed
% manually by the user upon each call !!!
%           
% dependency: GetSameChannels, format of the spike files
%
% DAP Apr 2015
% %-------------------------------%

% first clear the workspace:
clear

% variables are here:
%-------------------%
% number of neurons in a set (to model Ising for that group size):
non = 8;                  
% number of sets:
nos = 160;       
% maximum number of times a channel can be chosen [optional]:
max_times = [];
% names of files with recorded, filtered, and binned spikes:
filename{1} = './results_spikes/chip136_0_spikes_new_prob10_bin5_filt.mat';
filename{2} = './results_spikes/chip136_1_spikes_new_prob10_bin5_filt.mat';
filename{3} = './results_spikes/chip136_2_spikes_new_prob10_bin5_filt.mat';
filename{4} = './results_spikes/chip136_4_spikes_new_prob10_bin5_filt.mat';
filename{5} = './results_spikes/chip136_5_spikes_new_prob10_bin5_filt.mat';
filename{6} = './results_spikes/chip136_6_spikes_new_prob10_bin5_filt.mat';
% name for the output 'key':
fileout = './keys/chip136_all';
%-------------------%

% below no more adjustable parameters:
%-------------------%
addpath('./functions/')
NOF = length(filename);
disp(' ')
disp(['Randomly sampling  ',num2str(nos),'  sets of  ',num2str(non),'  neurons, using common channels of following files:'])
% loading the data:
chankeyAll = cell(NOF,1);
for i = 1:NOF
    filename_spikes = filename{i};
    eval(['load ',filename_spikes])
    disp(['Loaded file: ',filename_spikes])
    chankeyAll{i} = chankey_binned;
end
% if there is more than 1 file - compare channels to get common ones across
% all recordings:
if NOF>1
    chankey = GetSameChannels(chankeyAll);    
else
    chankey = chankeyAll{1};
end
% cell for the coordinates of each set of neurons:
key = cell(nos,1);               
% and array to monitor how many times channels were chosen already:
chancount=zeros(length(chankey),4);
chancount(:,1:3)=chankey(:,1:3);
% boolean variable for stopping the file write-out if the key is incomplete:
dont_write = 0; 
for i = 1:nos
    % first check there is enough channels:
    if size(chancount,1)<non
        if dont_write == 0
            disp(' ')
            disp('!! No file was written out, not enough channels left to sample from !!') 
            disp('Choose a different number of groups or relax the maximum number of times a channel is chosen.')
            disp(' ')
        end
        dont_write = 1;
        continue
    end
    % a random permutation:
    prop = randperm(size(chancount,1));
    pix = prop(1:non);
    % pick 'non' random channels and record which were chosen:
    for j = 1:non
        key{i}(j,:)= [chancount(pix(j),1:3)];
        prop2 = find(chancount(:,1)==key{i}(j,1));
        chancount(prop2,4) = chancount(prop2,4) + 1;
    end
    % find channels chosen too often already: 
    if ~(isempty(max_times))
        ind=find(chancount(:,4)>=max_times);
    else
        ind=find(chancount(:,4)>ceil(non*nos/length(chankey))-1);
    end
    % and get rid of them:
    if (length(ind)>=1)
        chancount(ind,:)=[];
    end   
end
% if 'key' is complete, display the layout and write out the file:
if dont_write == 0
    check = zeros(65,65);
    for i =1:nos; 
        for j=1:non;
            check(key{i}(j,2),key{i}(j,3)) = check(key{i}(j,2),key{i}(j,3))+1;
        end
    end
    figure
    pcolor(check')
    colors = jet(max(check(:))+1);
    colors(1,:) = [1 1 1];
    colormap(colors)
    colorbar
    shading flat
    disp(' ')
    disp(['Used for samples were ',num2str(length(find(check~=0))),'  out of  ',num2str(length(chankey)),'  channels.'])
    save(fileout,'key')
    disp('Resulting ''key'' with coordinates of channels in each group was saved:')
    disp(fileout)
    disp(' ')
end
%---------END----------%

