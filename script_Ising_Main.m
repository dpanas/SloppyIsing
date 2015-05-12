% %-------------------------------%
% script:   Ising_Main
%           Main script of the Ising model fitting process - requires both 
%           spike files and a 'key' file with coordinates of groups of
%           neurons that are to be analysed. Automatically saves results in
%           provided output files and plots basic convergence measures.
%           The learning rate and maximum number of iterations can be
%           adjusted. In cases of groups not converging as desired,
%           the group number is saved and script_Ising_Additional can be
%           used with different learning parameters to re-fit such groups.
%           The ising fit is exact (i.e. an iterative fitting method is
%           used), therefore it is not recommended for large groups of
%           neurons.
%
% !!! This is a script, not a function - parameters need to be changed
% manually by the user upon each call !!!
%
% dependency: ReadChannels, FitModelIndep, FitModelIsing, EvaluateModel,
%             InformationMeasures, CalculateModelStats, format of the spike
%             files, format of the group key file
%
% DAP Apr 2015
% %-------------------------------%

% first clear the workspace:
clear

% variables are here:
%-------------------%
% the learning rate for ising fit algorithm:
learn_rate_ising = 0.9;
% maximum iterations for the fitting:
max_iter = 60000;
% names of files with recorded, filtered, and binned spikes:
% (same as for script_SampleRandomGroups)
filenames{1} = './results_spikes/chip136_0_spikes_new_prob10_bin5_filt.mat';
filenames{2} = './results_spikes/chip136_1_spikes_new_prob10_bin5_filt.mat';
filenames{3} = './results_spikes/chip136_2_spikes_new_prob10_bin5_filt.mat';
filenames{4} = './results_spikes/chip136_4_spikes_new_prob10_bin5_filt.mat';
filenames{5} = './results_spikes/chip136_5_spikes_new_prob10_bin5_filt.mat';
filenames{6} = './results_spikes/chip136_6_spikes_new_prob10_bin5_filt.mat';
% the file with groups of neurons, provided by script_SampleRandomGroups:
keyname = './keys/chip136_all';
% and the names for files to be written out (check not to overwrite!)
fileouts{1} =  './results_ising/chip136_0_ising_8n_new_prob10_bin5_filt.mat';
fileouts{2} =  './results_ising/chip136_1_ising_8n_new_prob10_bin5_filt.mat';
fileouts{3} =  './results_ising/chip136_2_ising_8n_new_prob10_bin5_filt.mat';
fileouts{4} =  './results_ising/chip136_4_ising_8n_new_prob10_bin5_filt.mat';
fileouts{5} =  './results_ising/chip136_5_ising_8n_new_prob10_bin5_filt.mat';
fileouts{6} =  './results_ising/chip136_6_ising_8n_new_prob10_bin5_filt.mat';
%-------------------%

% (below no more variables, just code)
%-------------------%
close all
addpath('./functions/')
eval(['load ',keyname])
NOF = length(filenames);
nos = size(key,1);
non = size(key{1},1);
disp(' ')
disp(['Computing Ising model fits for ',num2str(NOF),' files, ',num2str(nos),' groups of ',num2str(non),' neurons each.'])
% going through data files on the list:
for j = 1:NOF;
    % loading the spikes:
    filename_spikes = filenames{j};
    eval(['load ',filename_spikes])
    disp(' ')
    disp(['Loaded spikes: ',' '' ',filename_spikes,' '' '])
    N = rtime./binsize;
    % preparing space for the results:
    fields = cell(nos,1);            % cell for the values of the 'magnetic fields'
    interactions = cell(nos,1);              % cell for the values of the 'interactions' or 'functional connections'
    corrs = cell(nos,1);             % cell for the pairwise correlation (second order stats)
    means = cell(nos,1);             % cell for the means in the group (first order stats)
    rates_data=zeros(nos,2^non);         % array with the rates of pattern occurrence in each group
    rates_model=zeros(nos,2^non);        % array with the rates of pattern occurrence approximated by the ising model for each group
    rates_indep=zeros(nos,2^non);        % array with the rates of pattern occurrence approximated by the independent model
    time_total = [];                   % time total of computations in sec;
    E_ind = zeros(nos,1);            % entropy of the independent model pattern distribution
    E_mod = zeros(nos,1);            % entropy of the ising model pattern distribution
    E_dat = zeros(nos,1);            % entropy of the pattern distribution obtained from data (rates of pattern occurrence)
    SJ_ind = zeros(nos,1);              % Shannon-Jensen divergence between independent distribution and data distribution
    SJ_mod = zeros(nos,1);              % Shannon-Jensen divergence between ising distribution and data distribution
    I = [];                 % a variable to record which groups reched max number of iterations (and did not converge as desired)
    accur_mod = zeros(1,nos);           % variable to approximately measure the accuracy of model fitting the first order marginals
    accur_mod2 = zeros(1,nos);          % variable to approximately measure the accuracy of model fitting the second order marginals
    tstart = tic;
    % loop over number of sets and fit an Independent and an Ising model for each pre-selected group of neurons:
    for i = 1:nos    
        % read in data from chosen channels:
        [D]=ReadChannels(spikes_binned,ceil(rtime./binsize),key{i});
        % compute data stats and save for reference:
        s_data  = mean(D,1);  % assuming that columns hold channels
        ss_data = D'*D./ceil(rtime/binsize);
        corrs{i} = ss_data;
        means{i} = s_data;
        % fit independent model and compute basic information measures:
        [rat,ent] = FitModelIndep(s_data);
        rates_indep(i,:) = rat;
        E_ind(i) = ent;
        % fit the Ising model parameters:
        [H,J,ent2,c] = FitModelIsing(s_data,ss_data,learn_rate_ising,max_iter);
        E_mod(i) = ent2;
        % check if Ising model converged:
        if c<max_iter
            disp(['Set ',num2str(i),' converged.'])
        else
            disp(['Ising model did not converge in set ',num2str(i)])
            I = [I i];
        end   
        % and record the models parameters: 
        fields{i} = H;
        interactions{i} = J;
        % eveluate rates of data and the model:
        [dat_rate,mod_rate,~]=EvaluateModel(H,J,H,D);
        rates_data(i,:) = dat_rate;
        rates_model(i,:) = mod_rate;
        % and basic information measures:
        [~,~,~,sj] = InformationMeasures(rates_data(i,:),rates_indep(i,:));
        SJ_ind(i) = sj;
        [ent,~,~,sj] = InformationMeasures(rates_data(i,:),rates_model(i,:));
        E_dat(i) = ent;
        SJ_mod(i) = sj;
        % and the approximate accuracy of spike fit:
        [s_ising,ss_ising,~] = CalculateModelStats(fields{i},interactions{i});
        accur_mod(i) = max(abs(means{i}-s_ising)*N/2);
        accur_mod2(i) =  max(max( (corrs{i}*N +2*meshgrid((means{i}+1)*N/2) +2*meshgrid((means{i}+1)*N/2)' - N)./4 - (ss_ising*N +2*meshgrid((s_ising+1)*N/2) +2*meshgrid((s_ising+1)*N/2)' - N)./4 ));
    end
    % double-check the convergence (NaNs, very large differences in spike counts, unrealistic multiinformation ratio):
    prop = find(isnan(E_mod));
    disp([' Sets with NaNs: ',num2str(prop)])
    prop = [prop find(((E_ind-E_mod)./(E_ind-E_dat))>1)];
    prop = [prop find(((E_ind-E_mod)./(E_ind-E_dat))<0)];
    disp([' Sets with unrealistic multiinformation ratio: ',num2str(find(((E_ind-E_mod)./(E_ind-E_dat))>1)),' and ',num2str(find(((E_ind-E_mod)./(E_ind-E_dat))<0))])
    % add any groups that appeared to converge within the iterations but have very suspicious behaviour otherwise:
    for i = 1:length(prop)
        prop2 = find(I==prop(i));
        if isempty(prop2)
            I = [I prop(i)];
        end
    end
    I = sort(I);
    % plot the summary on convergence of fits:
    figure
    set(gcf,'position',[100 400 1500 400])
    subplot(131)
    plot((E_ind-E_dat),((E_ind-E_mod)./(E_ind-E_dat)),'.')
    title('Multiinformation ratio','fontsize',13)
    subplot(132)
    plot(SJ_ind,'k.')
    hold on
    plot(SJ_mod,'r.')
    title('ShannonJensen divergence to data','fontsize',13)
    legend('indep','ising')
    subplot(133)
    plot(accur_mod,'r.')
    hold on
    plot(accur_mod2,'c.')
    title('Spike marginals accuracy','fontsize',13)
    legend('spikes','co-spikes')
    % time of computations:
    time_total = toc(tstart);
    disp(' ')
    disp(['Script took ',num2str(time_total./60/60),' hours to compute.'])
    % file info and date of creation:
    Info_Ising = ['This file was created by the script_Ising_Main on ',date];
    % and finally, saving the results:
    filename_ising = fileouts{j};
    save(filename_ising,'Info_Ising','corrs','fields','filename_spikes','filename_ising','interactions','key','learn_rate_ising','max_iter','means','rates_data','rates_model','rates_indep','time_total','I')
    disp(' ')
    disp('Script finished, results written succesfully:')
    disp(filename_ising)
    disp(' ')
end
%---------END----------%