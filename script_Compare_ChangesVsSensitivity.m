% %-------------------------------%
% script:   Compare_ChangesVsSensitivity
%           Assuming that the first eigenvector of Fisher Information
%           Matrix is a good proxy for parameter sensitivity, this script
%           plots the mean absolute changes in parameters, spike counts and
%           eigenparameters between recordings as a function of sensitivity
%           bin (average relative change across all data points from a
%           given sensitivity range).
%           This can be done for a list of files in reference to either a
%           single basliene file, or to a list of paired baseline files.
%           It can plot for comparison randomly resampled and refit
%           data from same distributions (if such file is provided; this
%           gives an idea of limited sample size effect); and it can also
%           plot an average over 50 reshuffles of sets in the data (to give
%           an idea of how different from random are the changes).
%           Other options include specification of own sensitivity range to
%           use and plotting also separate figures with scatterplots of
%           relative changes vs sensitivity.
%
% !!! This is a script, not a function - parameters need to be changed
% manually by the user upon each call !!!
%
% dependancy: SetPanelAppearance, PlotColorScatter, GetBinnedStats,
%             colormaps.mat
%
%  DAP Apr 2015
% %-------------------------------%

% first clear the workspace:
clear

% variables are here:
%----------------------------%
% the file names of fisher results to be compared:
% either one baseline file, or a list of files as baselines for each
% comparison (that way consecutive recordings can be compared):
filenames{1} =  './results_fisher/chip136_0_fisher_8n_new_prob10_bin5_filt.mat';
filenames{2} =  './results_fisher/chip136_0_fisher_8n_new_prob10_bin5_filt.mat';
filenames{3} =  './results_fisher/chip136_0_fisher_8n_new_prob10_bin5_filt.mat';
filenames{4} =  './results_fisher/chip136_0_fisher_8n_new_prob10_bin5_filt.mat';
filenames{5} =  './results_fisher/chip136_0_fisher_8n_new_prob10_bin5_filt.mat';
% and a list of files to be compared to baseline(s):
filenames2{1} =  './results_fisher/chip136_1_fisher_8n_new_prob10_bin5_filt.mat';
filenames2{2} =  './results_fisher/chip136_2_fisher_8n_new_prob10_bin5_filt.mat';
filenames2{3} =  './results_fisher/chip136_4_fisher_8n_new_prob10_bin5_filt.mat';
filenames2{4} =  './results_fisher/chip136_5_fisher_8n_new_prob10_bin5_filt.mat';
filenames2{5} =  './results_fisher/chip136_6_fisher_8n_new_prob10_bin5_filt.mat';

% optional: the name of a refit-resample fisher results file:
% (keep it commented out if not needed)
 filename_rr = './results_fisher/chip136_0_fisher_resample-refit_8n_new_prob10_bin5_filt.mat';

% optional: whether to plot shuffled sets for comparison or not:
% (default setting, 1, is yes)
do_shuffle = 1;

% optional: bin edges for computing mean absolute relative change and
% variance of neuronal parameters within separate sensitivity bins:
% (otherwise a default will be used of percentiles)
bin_edges = 0:0.08:0.32; 

% optional: whether to plot separate figures with scatterplots:
% (default, 0, is not to plot, not to take up memory)
scatterplot = 0; 
%----------------------------%


% (below no more variables, just code)
%----------------------------%
close all
addpath('./functions')
load ./functions/colormaps.mat
NOF = length(filenames2);

cell_eigenparam = cell(1,NOF);    % cell to hold all pooled 'sensitivity' (measured in first rec) - i.e. eigenparameters - i.e. principal eigenvector entries (scaled by eigenvalue)
cell_rates = cell(1,NOF);         % cell to hold all pooled mean firing rate (for a couple of neurons, for one neuron it's just rate) (measured in first rec)
cell_param_diff = cell(1,NOF);    % cell to hold all pooled differences in parameters between two recordings (relative differences)
cell_spike_diff = cell(1,NOF);    % cell to hold all pooled differences in spike count between two recordings (relative differences)
cell_eig_diff = cell(1,NOF);      % cell to hold all pooled differences in eigenparameters between two recordings (relative differences)

% loop through files to collect the data:
for i_file=1:NOF
    
    % loading the second recording, to be compared to the first / previous one:
    eval(['load ',filenames2{i_file}])
    disp(' ')
    disp(['Loaded Fisher results: ',' '' ',filenames2{i_file},' '' '])
    eval(['load ',filename_ising])
    eval(['load ',filename_spikes])
    chankey_binned1 = chankey_binned;
    correlations1 = correlations;
    fisher_matrix1=fisher_matrix;
   % rates_model1 = rates_model;
   % interactions1 = interactions;
   % fields1 = fields;
    corrs1 = corrs;
    means1 = means;
    parameter_matrix1 = parameter_matrix;
    eigenparam_matrix1 = eigenparam_matrix;
    spike_matrix1 = spike_matrix;
    if ( (i_file==1) && (exist('filename_rr')) )
        eval(['load ',filename_rr])
        disp(' ')
        disp(['Loaded Fisher results: ',' '' ',filename_rr,' '' '])
        parameter_matrix2 = parameter_matrix;
        eigenparam_matrix2 = eigenparam_matrix;
        spike_matrix2 = spike_matrix;
        N_res = size(eigenparam_matrix,2);
        vec_eig3 = [];
        vec_spike3 = []; 
        vec_param3 = []; 
        vec_rate3 = [];  
        vec_sens3 = [];
    end
    
    % and loading the first / previous recording, treated as baseline for comparison:
    if length(filenames)==1
        eval(['load ',filenames{1}])
        disp(' ')
        disp(['Loaded Fisher results: ',' '' ',filenames{1},' '' '])
        eval(['load ',filename_ising])
        eval(['load ',filename_spikes])        
    elseif length(filenames)==length(filenames2)
        eval(['load ',filenames{i_file}])
        disp(' ')
        disp(['Loaded Fisher results: ',' '' ',filenames{i_file},' '' '])
        eval(['load ',filename_ising])
        eval(['load ',filename_spikes])                
    else 
        disp('Incorrect number of files for comparison')
    end
    
    
    N = rtime./binsize;         % number of bins
    nos = size(key,1);          % number of sets
    non = length(key{1});     % number of neurons in each set
    
    % auxilliary variables to be cleared later:
    vec_eig = [];
    vec_spike = [];
    vec_param = [];
    vec_rate = [];   
    vec_sens = [];
    % if shuffled data is to be created, generate 50 random permutations of
    % set order:
    N_shuff = 50;
    if do_shuffle==1
        for i_shuf=1:N_shuff
            shuffle_prop(i_shuf,:) = randperm(nos);   % random permutation for a shuffle
        end
        vec_eig2 = [];
        vec_spike2 = [];
        vec_param2 = [];
        vec_rate2 = []; 
        vec_sens2 = [];
        
    end
    
    % go through all sets to collect stats to plot:
    for i=1:nos
        prop_r = (means{i}+1)./2/binsize;
        [x,y]=meshgrid(prop_r);
        rates_matrix = (x+y)./2;        
        for k=1:non
            vec_rate = [vec_rate rates_matrix(k,k:end)];
            vec_sens = [vec_sens abs(eigenparam_matrix{i}(k,k:end))];
            vec_eig = [vec_eig ( abs(eigenparam_matrix{i}(k,k:end)) - abs(eigenparam_matrix1{i}(k,k:end)) )./( abs(eigenparam_matrix{i}(k,k:end)) + abs(eigenparam_matrix1{i}(k,k:end)) )];
            vec_param = [vec_param ( parameter_matrix{i}(k,k:end) - parameter_matrix1{i}(k,k:end) )./( abs(parameter_matrix{i}(k,k:end)) + abs(parameter_matrix1{i}(k,k:end)) )];
            vec_spike = [vec_spike ( spike_matrix{i}(k,k:end) - spike_matrix1{i}(k,k:end) )./( spike_matrix{i}(k,k:end) + spike_matrix1{i}(k,k:end) )];
            % creating shuffled data results:
            if ( (i_file==NOF) && (do_shuffle==1) )
                for i_shuf=1:N_shuff
                    vec_rate2 = [vec_rate2 rates_matrix(k,k:end)];
                    vec_sens2 = [vec_sens2 abs(eigenparam_matrix{i}(k,k:end))];
                    vec_eig2 = [vec_eig2  ( abs(eigenparam_matrix{i}(k,k:end)) - abs(eigenparam_matrix1{shuffle_prop(i_shuf,i)}(k,k:end)) )./( abs(eigenparam_matrix{i}(k,k:end)) + abs(eigenparam_matrix1{shuffle_prop(i_shuf,i)}(k,k:end)) )];
                    vec_spike2 = [vec_spike2 ( spike_matrix{i}(k,k:end) - spike_matrix1{shuffle_prop(i_shuf,i)}(k,k:end) )./( spike_matrix{i}(k,k:end) + spike_matrix1{shuffle_prop(i_shuf,i)}(k,k:end) )];
                    vec_param2 = [vec_param2 ( parameter_matrix{i}(k,k:end) - parameter_matrix1{shuffle_prop(i_shuf,i)}(k,k:end) )./(  abs(parameter_matrix{i}(k,k:end)) + abs(parameter_matrix1{shuffle_prop(i_shuf,i)}(k,k:end)) )];
                end
            end
            % and adding 
            if ( (length(filenames)==1) && (i_file==NOF) && (exist('filename_rr')) )
                for i_res=1:N_res
                    vec_rate3 = [vec_rate3 rates_matrix(k,k:end)];
                    vec_sens3 = [vec_sens3 abs(eigenparam_matrix{i}(k,k:end))];
                    vec_eig3 = [vec_eig3  ( abs(eigenparam_matrix{i}(k,k:end)) - abs(eigenparam_matrix2{i,i_res}(k,k:end)) )./( abs(eigenparam_matrix{i}(k,k:end)) + abs(eigenparam_matrix2{i,i_res}(k,k:end)) )];
                    vec_spike3 = [vec_spike3 ( spike_matrix{i}(k,k:end) - spike_matrix2{i,i_res}(k,k:end) )./( spike_matrix{i}(k,k:end) + spike_matrix2{i,i_res}(k,k:end) )];
                    vec_param3 = [vec_param3 ( parameter_matrix{i}(k,k:end) - parameter_matrix2{i,i_res}(k,k:end) )./(  abs(parameter_matrix{i}(k,k:end)) + abs(parameter_matrix2{i,i_res}(k,k:end)) )];
                end  
            elseif ( (length(filenames)==length(filenames2)) && (i_file==1) && (exist('filename_rr')) )
                for i_res=1:N_res
                    vec_rate3 = [vec_rate3 rates_matrix(k,k:end)];
                    vec_sens3 = [vec_sens3 abs(eigenparam_matrix{i}(k,k:end))];
                    vec_eig3 = [vec_eig3  ( abs(eigenparam_matrix{i}(k,k:end)) - abs(eigenparam_matrix2{i,i_res}(k,k:end)) )./( abs(eigenparam_matrix{i}(k,k:end)) + abs(eigenparam_matrix2{i,i_res}(k,k:end)) )];
                    vec_spike3 = [vec_spike3 ( spike_matrix{i}(k,k:end) - spike_matrix2{i,i_res}(k,k:end) )./( spike_matrix{i}(k,k:end) + spike_matrix2{i,i_res}(k,k:end) )];
                    vec_param3 = [vec_param3 ( parameter_matrix{i}(k,k:end) - parameter_matrix2{i,i_res}(k,k:end) )./(  abs(parameter_matrix{i}(k,k:end)) + abs(parameter_matrix2{i,i_res}(k,k:end)) )];
                end  
            end
        end              
    end
    
    % and save it in cell:
    cell_eigenparam{i_file} = vec_sens;
    cell_rates{i_file} = vec_rate;
    cell_param_diff{i_file} = vec_param;
    cell_spike_diff{i_file} = vec_spike;
    cell_eig_diff{i_file} = vec_eig;
    if ( (i_file==NOF) && (do_shuffle==1) )
        cell_eigenparam{NOF+1} = vec_sens2;
        cell_rates{NOF+1} = vec_rate2;
        cell_param_diff{NOF+1} = vec_param2;
        cell_spike_diff{NOF+1} = vec_spike2;
        cell_eig_diff{NOF+1} = vec_eig2;
    end
    if ( (i_file==NOF) && exist('filename_rr'))
        cell_eigenparam{NOF+2} = vec_sens3;
        cell_rates{NOF+2} = vec_rate3;
        cell_param_diff{NOF+2} = vec_param3;
        cell_spike_diff{NOF+2} = vec_spike3;
        cell_eig_diff{NOF+2} = vec_eig3;
    end
    
end
% clear the cluttering variables no longer needed:
clear vec_rate vec_eig vec_eig1 vec_eig2 vec_param vec_param1 vec_param2 vec_spike vec_spike1 vec_spike2 shuffle_prop prop_r rates_matrix x y 
clear eigenparam_matrix eigenparam_matrix1 spike_matrix spike_matrix1

% if edges for binning are not provided, then resort to percentiles in
% the data: (this might though not be the best method for highly skewed
% distributions of sensitivity values)
if isempty(bin_edges)
    bin_edges = prctile(cell_eigenparam{i_file}(:),[2.5 25 50 75 97.5]);
    bin_edges = [0 bin_edges max(cell_eigenparam{i_file}(:))+0.1];  
end
% set colors:
%load ./functions/colormaps.mat
eval(['col = colorcell_',num2str(NOF+1),';'])
col{end+1} = [0 0 0];  % adding extra color for shuffled data, instead of color of first rec;
figure
set(gcf,'position',[150 150 1200 650])
label_for_legend = [];
string_for_legend = [];
% and now plot the mean and variance change in provided bins:
for i_file=1:NOF    
    % panels with means:
    subplot(231)
    hold on
    [bin,M_param,V_param] = GetBinnedStats(cell_param_diff{i_file}(:),cell_eigenparam{i_file}(:),bin_edges);
    plot(bin,M_param,'lineWidth',2,'color',col{i_file+1});
    subplot(232)
    hold on
    [bin,M_spike,V_spike] = GetBinnedStats(cell_spike_diff{i_file}(:),cell_eigenparam{i_file}(:),bin_edges);
    plot(bin,M_spike,'lineWidth',2,'color',col{i_file+1})
    subplot(233)
    hold on
    [bin,M_eigen,V_eigen] = GetBinnedStats(cell_eig_diff{i_file}(:),cell_eigenparam{i_file}(:),bin_edges);
    h = plot(bin,M_eigen,'lineWidth',2,'color',col{i_file+1});
    label_for_legend(i_file) = h;
    string_for_legend = [string_for_legend,',''rec ',num2str(i_file),' '''];
    % panels with variances:
    subplot(234)
    hold on
    plot(bin,V_param,'lineWidth',2,'color',col{i_file+1})
    subplot(235)
    hold on
    plot(bin,V_spike,'lineWidth',2,'color',col{i_file+1})
    subplot(236)
    hold on
    plot(bin,V_eigen,'lineWidth',2,'color',col{i_file+1})   
        
%     [cc,ccpp] = corr(vec_eig',vec_rate');
%     disp(' ')
%     disp(['Correlation between sensitivity and firing rates: ',num2str(cc),' and its significance: ',num2str(ccpp),' .'])
%     disp(' ')
end
% if shuffled, then add shuffled (in black)
if (do_shuffle==1)
    subplot(231)
    hold on
    [bin,M_param,V_param] = GetBinnedStats(cell_param_diff{NOF+1}(:),cell_eigenparam{NOF+1}(:),bin_edges);
    plot(bin,M_param,'lineWidth',2,'color',col{end})
    subplot(232)
    hold on
    [bin,M_spike,V_spike] = GetBinnedStats(cell_spike_diff{NOF+1}(:),cell_eigenparam{NOF+1}(:),bin_edges);
    plot(bin,M_spike,'lineWidth',2,'color',col{end})
    subplot(233)
    hold on
    [bin,M_eigen,V_eigen] = GetBinnedStats(cell_eig_diff{NOF+1}(:),cell_eigenparam{NOF+1}(:),bin_edges);
    h = plot(bin,M_eigen,'lineWidth',2,'color',col{end});
    label_for_legend = [label_for_legend  h];
    string_for_legend = [string_for_legend,',''shuffled'''];
    % panels with variances:
    subplot(234)
    hold on
    plot(bin,V_param,'lineWidth',2,'color',col{end})
    subplot(235)
    hold on
    plot(bin,V_spike,'lineWidth',2,'color',col{end})
    subplot(236)
    hold on
    plot(bin,V_eigen,'lineWidth',2,'color',col{end})   
end
% if file with resample-refit was provided, also plot average over
% resamples along with the data:
if  (exist('filename_rr'))
    subplot(231)
    hold on
    [bin,M_param,V_param] = GetBinnedStats(cell_param_diff{NOF+2}(:),cell_eigenparam{NOF+2}(:),bin_edges);
    plot(bin,M_param,'lineWidth',2,'color',col{1})
    subplot(232)
    hold on
    [bin,M_spike,V_spike] = GetBinnedStats(cell_spike_diff{NOF+2}(:),cell_eigenparam{NOF+2}(:),bin_edges);
    plot(bin,M_spike,'lineWidth',2,'color',col{1})
    subplot(233)
    hold on
    [bin,M_eigen,V_eigen] = GetBinnedStats(cell_eig_diff{NOF+2}(:),cell_eigenparam{NOF+2}(:),bin_edges);
    h = plot(bin,M_eigen,'lineWidth',2,'color',col{1});
    label_for_legend = [label_for_legend  h];
    string_for_legend = [string_for_legend,',''refit'''];
    % panels with variances:
    subplot(234)
    hold on
    plot(bin,V_param,'lineWidth',2,'color',col{1})
    subplot(235)
    hold on
    plot(bin,V_spike,'lineWidth',2,'color',col{1})
    subplot(236)
    hold on
    plot(bin,V_eigen,'lineWidth',2,'color',col{1})  
end
% and fix plot appearance:
subplot(231)
SetPanelAppearance('Parameters','sensitivity',[],[],'mean relative change',[],[])
subplot(232)
SetPanelAppearance('Spike counts','sensitivity',[],[],'mean relative change',[],[])
subplot(233)
SetPanelAppearance('Eigenparameters','sensitivity',[],[],'mean relative change',[],[])
eval(['legend(label_for_legend',string_for_legend,')'])
subplot(234)
SetPanelAppearance('Parameters','sensitivity',[],[],'var relative change',[],[])
subplot(235)
SetPanelAppearance('Spike counts','sensitivity',[],[],'var relative change',[],[])
subplot(236)
SetPanelAppearance('Eigenparameters','sensitivity',[],[],'var relative change',[],[])
eval(['legend(label_for_legend',string_for_legend,')'])

if scatterplot==1
    edges2 = logspace(-1,1,64);
    col2 = color_map;
    for i_file=1:NOF
        figure
        set(gcf,'position',[100 300 1300 360])
        subplot(131)
        PlotColorScatter(cell_eigenparam{i_file}(:),cell_param_diff{i_file}(:),cell_rates{i_file}(:),edges2,col2)
        SetPanelAppearance('Parameters','sensitivity',[],[],'relative change',[-1:0.5:1],[])
        subplot(132)
        PlotColorScatter(cell_eigenparam{i_file}(:),cell_spike_diff{i_file}(:),cell_rates{i_file}(:),edges2,col2)
        SetPanelAppearance('Spike counts','sensitivity',[],[],'relative change',[-1:0.5:1],[])
        subplot(133)
        PlotColorScatter(cell_eigenparam{i_file}(:),cell_eig_diff{i_file}(:),cell_rates{i_file}(:),edges2,col2)
        SetPanelAppearance('Eigenparameters','sensitivity',[],[],'relative change',[-1:0.5:1],[])
    end
    if ( do_shuffle==1 )
        figure
        set(gcf,'position',[100 300 1300 360])
        subplot(131)
        PlotColorScatter(cell_eigenparam{NOF+1}(:),cell_param_diff{NOF+1}(:),cell_rates{NOF+1}(:),edges2,col2)
        SetPanelAppearance('Parameters','sensitivity',[],[],'relative change',[-1:0.5:1],[])
        subplot(132)
        PlotColorScatter(cell_eigenparam{NOF+1}(:),cell_spike_diff{NOF+1}(:),cell_rates{NOF+1}(:),edges2,col2)
        SetPanelAppearance('Spike counts','sensitivity',[],[],'relative change',[-1:0.5:1],[])
        subplot(133)
        PlotColorScatter(cell_eigenparam{NOF+1}(:),cell_eig_diff{NOF+1}(:),cell_rates{NOF+1}(:),edges2,col2)
        SetPanelAppearance('Eigenparameters','sensitivity',[],[],'relative change',[-1:0.5:1],[])        
    end
end

%-------------END---------------%
