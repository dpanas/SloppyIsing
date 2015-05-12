% %-------------------------------%
% script:   FIM_Decomposition_ResampleRefit
%           Additional script for further analysis of Ising re-sampling
%           results, so that also FIM decompositions are obtained for the
%           re-fit samples and can be compared. The only adjustable
%           parameters are the file names.
%
% !!! This is a script, not a function - parameters need to be changed
% manually by the user upon each call !!!
%
% dependancy: CalculateModelStats, FisherInformationMatrix,
%             GiniCoefficient, format of the ising result files
%
% DAP Apr 2015
% %-------------------------------%

% first clear the workspace:
clear

% variables are here:
%-------------------%
% filename of the re-sampled ising results:
filename = './results_ising/chip136_0_resample-refit_8n_new_prob10_bin5_filt.mat';
% filename for the output 
fileouts = './results_fisher/chip136_0_fisher_resample-refit_8n_new_prob10_bin5_filt.mat';
%-------------------%

% (below no more variables, just code)
%-------------------%
close all
addpath('./functions/')
eval(['load ',filename])
eval(['load ',filename_spikes])
nos = size(key,1);
non = length(key{1});
disp(' ')
disp(['Computing Fisher Information Matrix Decomposition for re-sampled data, ',num2str(no_samples),' re-samples for ',num2str(nos),' sets.'])
disp(' ')
disp(['Loaded Ising re-sampling results: ',filename_ising_resample])
disp(['   and the spike file ',filename_spikes])
disp(' ')
N = ceil(rtime/binsize);
% preparing space for the results:
fisher_matrix = cell(nos,no_samples);             % cell holding the Fisher Information Matrix for each model / set of neurons
parameter_matrix = cell(nos,no_samples);          % cell of matrices with model parameters (fields on diagonal, interactions off diagonal)
eigenparam_matrix = cell(nos,no_samples);         % cell of matrices with eigenparameters, i.e. entries of the first eigenvector, corresponding to each individual model parameter from parameter matrix, scaled by eigenvalue
spike_matrix = cell(nos,no_samples);             % cell of matrices with spike and co-spike counts, again in a layout corresponding to the parameter matrix
gini_FIM = -1*ones(no_samples,nos);             % Gini coefficient of the FIM matrix (i.e. how sparse is the matrix?)
% loop over modelled re-samples:
for k=1:no_samples
    % loop over groups of neurons:
    for i=1:nos
        [s_ising,ss_ising,~] = CalculateModelStats(fields_resample{i,k},interactions_resample{i,k});
        fisher_matrix{i,k} = FisherInformationMatrix(rates_model_resample{i,k},s_ising,ss_ising);
        % perform eigenvalue decomposition of the Fisher Information Matrix
        [eve,eva] = eig(fisher_matrix{i,k});
        % re-format results of the matlab eig:
        eva = eva([1:(length(eva)+1):end]);   
        % sort eigenvalues and eigenvectors in descending order:
        [eva,idx] = sort(eva,'descend');    
        eve = eve(:,idx);
        eigenvalues{i,k} = eva; 
        eigenvectors{i,k} = eve;
        gini_FIM(k,i) = GiniCoefficient(fisher_matrix{i,k}(:));
        % reshape the eigenvector to the format of matrix, with single
        % parameters on the diagonal:
        first_eig(1:(non+1):(non*non)) = eigenvectors{i,k}(1:non,1);
        first_eig = reshape(first_eig,[non non]);
        counter = non;
        for l=1:non
            first_eig(l,l+1:end) = eigenvectors{i,k}(counter+1:counter+non-l,1) ;
            first_eig(l+1:end,l) = eigenvectors{i,k}(counter+1:counter+non-l,1) ;
            counter = counter + non - l;
        end   
        % and renormalize the first eigenvector by first eigenvalue:
        eigenparam_matrix{i,k} = eigenvalues{i,k}(1).*first_eig;  
        % get the spike matrix:
        prop = (corrs_resample{i,k}*N +2*meshgrid((means_resample{i,k}+1)*N/2) +2*meshgrid((means_resample{i,k}+1)*N/2)' - N)./4 ;
        prop = prop.*(ones(non,non)-eye(non))+(eye(non).*(((means_resample{i,k}'+1)*N/2)*ones(1,non)));
        spike_matrix{i,k} = prop;
        % and the parameter matrix:
        prop = interactions_resample{i,k};
        prop = prop.*(ones(non,non)-eye(non))+eye(non).*(fields_resample{i,k}'*ones(1,non));
        parameter_matrix{i,k} = prop;
    end
end
% file info and date of creation:
Info_Fisher = ['This file was created by the script FIM_Decomposition_ResampleRefit on ',date];   
% and save the results:
filename_fisher_resample = fileouts;
save(filename_fisher_resample,'Info_Fisher','filename_fisher_resample','filename_spikes','filename_ising_resample','fisher_matrix','eigenvalues','eigenvectors','spike_matrix','parameter_matrix','gini_FIM','eigenparam_matrix')
disp(' ')
disp(['Script finished, results written succesfully:  ',filename_fisher_resample])
disp(' ')
%---------END----------%
