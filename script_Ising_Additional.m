% %-------------------------------%
% script:   Ising_Additional
%           Additional script for Ising model fitting, to take care of
%           groups of neurons that did not converge satisfactorily in the
%           first attempt. Requires the ising results files for correction
%           and new learning rate and / or maximum number of iterations.
%           Overwrites the previous version of the file, changing only the
%           problematic groups. Plots basic convergence measures.
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
% the new learning rate for ising fit algorithm:
learn_rate_ising2 = 0.95;
% maximum iterations for the fitting:
max_iter2 = 60000;
% names of files with unconverged / unsatisfactorily converged groups:
filenames{1} = './results_ising/chip136_0_ising_8n_new_prob10_bin5_filt.mat';
filenames{2} = './results_ising/chip136_1_ising_8n_new_prob10_bin5_filt.mat';
filenames{3} = './results_ising/chip136_2_ising_8n_new_prob10_bin5_filt.mat';
filenames{4} = './results_ising/chip136_4_ising_8n_new_prob10_bin5_filt.mat';
filenames{5} = './results_ising/chip136_5_ising_8n_new_prob10_bin5_filt.mat';
filenames{6} = './results_ising/chip136_6_ising_8n_new_prob10_bin5_filt.mat';
%-------------------%

% (below no more variables, just code)
%-------------------%
close all
addpath('./functions/')
NOF = length(filenames);
disp(' ')
disp(['Computing additional Ising model fits for ',num2str(NOF),' files.'])
% going through the ising files on the list:
for j = 1:NOF;
    % loading the ising result files:
    filename_ising = filenames{j};
    eval(['load ',filename_ising])
    disp(' ')
    disp(['Loaded ising results for additional fits: ',' '' ',filename_ising,' '' '])
    disp(['Original learning rate: ',num2str(learn_rate_ising),', and new learning rate: ',num2str(learn_rate_ising2),'.'])
    disp(['Original max number of iterations: ',num2str(max_iter),', and new max number of iterations: ',num2str(max_iter2),'.'])
    % and loading the corresponding spike file:
    eval(['load ',filename_spikes])
    N = rtime./binsize;
    nos = size(key,1);
    non = size(key{1},1);
    % changing the learning rate and max iterations for groups that did not converge:
    learn_rate_ising = [learn_rate_ising learn_rate_ising2];
    max_iter = [max_iter max_iter2];
    Icurr = I;
    I = [];
    tstart = tic;
    % loop over the unconverged sets and fit again an Ising model:
    for i = Icurr   
        % fit the Ising model parameters:
        [H,J,ent2,c] = FitModelIsing(means{i},corrs{i},learn_rate_ising(end),max_iter(end));
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
        % eveluate rates the model:
        [~,mod_rate,~]=EvaluateModel(H,J,H,zeros(rtime./binsize,non));
        rates_model(i,:) = mod_rate;
    end
    % after the refits, again re-evaluate convergence:
    accur_mod = zeros(1,nos);           % variable to approximately measure the accuracy of model fitting the first order marginals
    accur_mod2 = zeros(1,nos);          % variable to approximately measure the accuracy of model fitting the second order marginals
    for i=1:nos
        % and basic information measures:
        [~,ent,~,sj]= InformationMeasures(rates_data(i,:),rates_indep(i,:));
        E_ind(i) = ent;
        SJ_ind(i) = sj;
        [ent,ent2,~,sj] = InformationMeasures(rates_data(i,:),rates_model(i,:));
        E_dat(i) = ent;
        E_mod(i) = ent2;
        SJ_mod(i) = sj;
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
    time_total_corrections = toc(tstart);
    disp(['Script took ',num2str(time_total_corrections./60/60),' hours to compute.'])
    % add a note to the file info :
    Info_Ising = [Info_Ising,', and modified on ',date];
    % saving the results (overwriting previous version):
    save(filename_ising,'Info_Ising','corrs','fields','filename_spikes','filename_ising','interactions','key','learn_rate_ising','max_iter','means','rates_data','rates_model','rates_indep','time_total','I')
    disp(' ')
    disp('Script finished, results written succesfully:')
    disp(filename_ising)
    disp(' ')
end
%---------END----------%