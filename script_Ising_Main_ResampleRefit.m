% %-------------------------------%
% script:   Ising_Main_ResampleRefit
%           Additional script for the Ising analysis: re-sampling spike
%           trains from the fit distributions and re-fitting them again, in
%           order to verify that the results are not due to overfitting.
%           The number of times each group is resampled 'no_samples' is
%           adjustable. For recordings that are generally similar,
%           re-sampling only one of them should suffice as such a safety
%           check.
%
% !!! This is a script, not a function - parameters need to be changed
% manually by the user upon each call !!!
%
% dependency: SampleFromIsing, FitModelIndep, FitModelIsing, EvaluateModel,
%             InformationMeasures, CalculateModelStats, format of the ising
%             file with results
%
% DAP Apr 2015
% %-------------------------------%

% first clear the workspace:
clear

% variables are here:
%-------------------%
% name of the file (e.g. baseline one) with Ising results to be re-sampled
% and re-fit again:
filename = './results_ising/chip136_0_ising_8n_new_prob10_bin5_filt.mat';
% name of the ouptut file with results of the re-fitting:
fileouts = './results_ising/chip136_0_resample-refit_8n_new_prob10_bin5_filt.mat';
% number of re-samples to be used:
no_samples = 20;
%-------------------%

% (below no more variables, just code)
%-------------------%
close all
addpath('./functions/')
eval(['load ',filename])
eval(['load ',filename_spikes])
nos = size(key,1);
non = size(key{1},1);
disp(' ')
disp(['Computing Ising model fits to re-sampled data, ',num2str(no_samples),' re-samples for ',num2str(nos),' groups of ',num2str(non),' neurons each.'])
disp(' ')
disp(['Loaded Ising results: ',filename_ising])
disp(['   and the spike file ',filename_spikes])
disp(' ')
N = rtime./binsize;
% preparing space for the results:
fields_resample = cell(nos,no_samples);            % cell for the values of the 'magnetic fields'
interactions_resample = cell(nos,no_samples);              % cell for the values of the 'interactions' or 'functional connections'
corrs_resample = cell(nos,no_samples);             % cell for the pairwise correlation (second order stats)
means_resample = cell(nos,no_samples);             % cell for the means in the group (first order stats)
rates_data_resample=cell(nos,no_samples);         % array with the rates of pattern occurrence in each group
rates_model_resample=cell(nos,no_samples);        % array with the rates of pattern occurrence approximated by the ising model for each group
rates_indep_resample=cell(nos,no_samples);        % array with the rates of pattern occurrence approximated by the independent model
time_total_resample = [];                   % time total of computations in sec;
E_ind_resample = zeros(nos,no_samples);            % entropy of the independent model pattern distribution
E_mod_resample = zeros(nos,no_samples);            % entropy of the ising model pattern distribution
E_dat_resample = zeros(nos,no_samples);            % entropy of the pattern distribution obtained from data (rates of pattern occurrence)
SJ_ind_resample = zeros(nos,no_samples);              % Shannon-Jensen divergence between independent distribution and data distribution
SJ_mod_resample = zeros(nos,no_samples);              % Shannon-Jensen divergence between ising distribution and data distribution
I_resample = [];                 % a variable to record which groups reched max number of iterations (and did not converge as desired)
accur_mod_resample = zeros(nos,no_samples);           % variable to approximately measure the accuracy of model fitting the first order marginals
accur_mod2_resample = zeros(nos,no_samples);          % variable to approximately measure the accuracy of model fitting the second order marginals
tstart = tic;  
learn_rate_ising_resample = learn_rate_ising(1);
max_iter_resample = max_iter(1);
% loop over number of sets, re-sample spike trains from model disribution, and fit again Independent and Ising, no_sample times:
for i = 1:nos    
    for k=1:no_samples
        % instead of reading in data, get fake data sampled from the Ising
        % model:
        [~,D] = SampleFromIsing(rates_model(i,:),ceil(N),1);
        % compute data stats and save for reference:
        s_data  = mean(D,1);  % assuming that columns hold channels
        ss_data = D'*D./ceil(rtime/binsize);
        corrs_resample{i,k} = ss_data;
        means_resample{i,k} = s_data;
        % fit independent model and compute basic information measures:
        [rat,ent] = FitModelIndep(s_data);
        rates_indep_resample{i,k} = rat;
        E_ind_resample(i,k) = ent;
        % fit the Ising model parameters:
        [H,J,ent2,c] = FitModelIsing(s_data,ss_data,learn_rate_ising_resample(1),max_iter_resample(1));
        E_mod_resample(i,k) = ent2;
        % check if Ising model converged:
        if c<max_iter_resample
            disp(['Set ',num2str(i),' re-sample ',num2str(k),' converged.'])
        else
            disp(['Ising model did not converge in set ',num2str(i),' re-sample ',num2str(k)])
            I_resample = [I_resample; i k];
        end   
        % and record the models parameters: 
        fields_resample{i,k} = H;
        interactions_resample{i,k} = J;
        % eveluate rates of data and the model:
        [dat_rate,mod_rate,~]=EvaluateModel(H,J,H,D);
        rates_data_resample{i,k} = dat_rate;
        rates_model_resample{i,k} = mod_rate;
        % and basic information measures:
        [~,~,~,sj] = InformationMeasures(rates_data_resample{i,k},rates_indep_resample{i,k});
        SJ_ind_resample(i,k) = sj;
        [ent,~,~,sj] = InformationMeasures(rates_data_resample{i,k},rates_model_resample{i,k});
        E_dat_resample(i,k) = ent;
        SJ_mod_resample(i,k) = sj;
        % and the approximate accuracy of spike fit:
        [s_ising,ss_ising,~] = CalculateModelStats(fields_resample{i,k},interactions_resample{i,k});
        accur_mod_resample(i,k) = max(abs(means_resample{i,k}-s_ising)*N/2);
        accur_mod2_resample(i,k) =  max(max( (corrs_resample{i,k}*N +2*meshgrid((means_resample{i,k}+1)*N/2) +2*meshgrid((means_resample{i,k}+1)*N/2)' - N)./4 - (ss_ising*N +2*meshgrid((s_ising+1)*N/2) +2*meshgrid((s_ising+1)*N/2)' - N)./4 ));
    end
end
% double-check the convergence (NaNs, very large differences in spike counts, unrealistic multiinformation ratio):
for i=1:nos
    prop = find(isnan(E_mod_resample(i,:)));
    if ~isempty(prop)
        disp([' Sets with NaNs: ',num2str(i),' and re-samples ',num2str(prop)])
    end 
    prop = [prop find(((E_ind_resample(i,:)-E_mod_resample(i,:))./(E_ind_resample(i,:)-E_dat_resample(i,:)))>1) find(((E_ind_resample(i,:)-E_mod_resample(i,:))./(E_ind_resample(i,:)-E_dat_resample(i,:)))<0)];
    if ~isempty( [find(((E_ind_resample(i,:)-E_mod_resample(i,:))./(E_ind_resample(i,:)-E_dat_resample(i,:)))>1) find(((E_ind_resample(i,:)-E_mod_resample(i,:))./(E_ind_resample(i,:)-E_dat_resample(i,:)))<0)] )
        disp([' Sets with unrealistic multiinformation ratio: ',num2str(find(((E_ind_resample(i,:)-E_mod_resample(i,:))./(E_ind_resample(i,:)-E_dat_resample(i,:)))>1)),' and ',num2str(find(((E_ind_resample(i,:)-E_mod_resample(i,:))./(E_ind_resample(i,:)-E_dat_resample(i,:)))<0))])
    end
    if ~isempty(prop)
        prop2 = find(I_resample(:,1)==i);
        I2 = I_resample(prop2,2);
        % add any groups that appeared to converge within the iterations but have very suspicious behaviour otherwise:
        for l = 1:length(prop)
            prop3 = find(I2==prop(l));
            if isempty(prop3)
                I_resample = [I_resample; i prop(l)];
            end
        end
    end
end
% plot the summary on convergence of fits:
figure
set(gcf,'position',[100 400 1500 400])
[X,Y]=meshgrid([1:1:no_samples],[1:1:nos]);
subplot(131)
plot3(X,(E_ind_resample-E_dat_resample),((E_ind_resample-E_mod_resample)./(E_ind_resample-E_dat_resample)),'b.')
title('Multiinformation ratio','fontsize',13)
subplot(132)
plot3(X,Y,SJ_ind_resample,'k.')
hold on
plot3(X,Y,SJ_mod_resample,'r.')
title('ShannonJensen divergence to data','fontsize',13)
subplot(133)
plot3(X,Y,accur_mod_resample,'r.')
hold on
plot3(X,Y,accur_mod2_resample,'c.')
title('Spike marginals accuracy','fontsize',13)
% time of computations:
time_total_resample = toc(tstart);
disp(' ')
disp(['Script took ',num2str(time_total_resample./60/60),' hours to compute.'])
% file info and date of creation:
Info_Ising_RefitResample = ['This file was created by the script_Ising_Main_ResampleRefit on ',date,' for the file ',filename_ising];
% and finally, saving the results:
filename_ising_resample = fileouts;
save(filename_ising_resample,'Info_Ising_RefitResample','no_samples','corrs_resample','fields_resample','filename_ising','filename_spikes','filename_ising_resample','interactions_resample','key','learn_rate_ising','max_iter','means_resample','rates_data_resample','rates_model_resample','rates_indep_resample','time_total_resample','I_resample')
disp(' ')
disp('Script finished, results written succesfully:')
disp(filename_ising_resample)
disp(' ')

%---------END----------%