% %-------------------------------%
% script:   FIM_Decomposition
%           Script for further analysis of Ising results - computing the
%           Fisher Information Matrix for each of the fitted models and
%           decomposing the matrix into eigenvectors, to see which
%           parameters / combinations of parameters are key to behaviour of
%           the group. Only adjustable parameters are file names.
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
% names of files with ising results:
filenames{1} = './results_ising/chip136_0_ising_8n_new_prob10_bin5_filt.mat';
filenames{2} = './results_ising/chip136_1_ising_8n_new_prob10_bin5_filt.mat';
filenames{3} = './results_ising/chip136_2_ising_8n_new_prob10_bin5_filt.mat';
filenames{4} = './results_ising/chip136_4_ising_8n_new_prob10_bin5_filt.mat';
filenames{5} = './results_ising/chip136_5_ising_8n_new_prob10_bin5_filt.mat';
filenames{6} = './results_ising/chip136_6_ising_8n_new_prob10_bin5_filt.mat';
% and the names for files to be written out (check not to overwrite!)
fileouts{1} =  './results_fisher/chip136_0_fisher_8n_new_prob10_bin5_filt.mat';
fileouts{2} =  './results_fisher/chip136_1_fisher_8n_new_prob10_bin5_filt.mat';
fileouts{3} =  './results_fisher/chip136_2_fisher_8n_new_prob10_bin5_filt.mat';
fileouts{4} =  './results_fisher/chip136_4_fisher_8n_new_prob10_bin5_filt.mat';
fileouts{5} =  './results_fisher/chip136_5_fisher_8n_new_prob10_bin5_filt.mat';
fileouts{6} =  './results_fisher/chip136_6_fisher_8n_new_prob10_bin5_filt.mat';
%-------------------%

% (below no more variables, just code)
%-------------------%
close all
addpath('./functions/')
NOF = length(filenames);
disp(' ')
disp(['Fisher Information Matrix Decomposition script for ',num2str(NOF),' files.'])
% going through the files on the list:
for j=1:NOF
    % loading the result ising files:
    filename_ising = filenames{j};
    eval(['load ',filename_ising])
    disp(' ')
    disp(['Loaded ising results: ',' '' ',filename_ising,' '' '])
    eval(['load ',filename_spikes])
    disp(['Loaded spikes: ',' '' ',filename_spikes,' '' '])
    nos = size(key,1);
    non = length(key{1});
    N = ceil(rtime/binsize);
    % preparing space for the results:
    fisher_matrix = cell(nos,1);             % cell holding the Fisher Information Matrix for each model / set of neurons
    parameter_matrix = cell(nos,1);          % cell of matrices with model parameters (fields on diagonal, interactions off diagonal)
    eigenparam_matrix = cell(nos,1);         % cell of matrices with eigenparameters, i.e. entries of the first eigenvector, corresponding to each individual model parameter from parameter matrix, scaled by eigenvalue
    spike_matrix = cell(nos,1);             % cell of matrices with spike and co-spike counts, again in a layout corresponding to the parameter matrix
    gini_FIM = -1*ones(1,nos);             % Gini coefficient of the FIM matrix (i.e. how sparse is the matrix?)
    % loop over modelled groups of neurons:
    for i=1:nos
        [s_ising,ss_ising,~] = CalculateModelStats(fields{i},interactions{i});
        fisher_matrix{i} = FisherInformationMatrix(rates_model(i,:),s_ising,ss_ising);
        % perform eigenvalue decomposition of the Fisher Information Matrix
        [eve,eva] = eig(fisher_matrix{i});
        % re-format results of the matlab eig:
        eva = eva([1:(length(eva)+1):end]);   
        % sort eigenvalues and eigenvectors in descending order:
        [eva,idx] = sort(eva,'descend');    
        eve = eve(:,idx);
        eigenvalues{i} = eva; 
        eigenvectors{i} = eve;
        gini_FIM(i) = GiniCoefficient(fisher_matrix{i}(:));
        % reshape the eigenvector to the format of matrix, with single
        % parameters on the diagonal:
        first_eig(1:(non+1):(non*non)) = eigenvectors{i}(1:non,1);
        first_eig = reshape(first_eig,[non non]);
        counter = non;
        for k=1:non
            first_eig(k,k+1:end) = eigenvectors{i}(counter+1:counter+non-k,1) ;
            first_eig(k+1:end,k) = eigenvectors{i}(counter+1:counter+non-k,1) ;
            counter = counter + non - k;
        end   
        % and renormalize the first eigenvector by first eigenvalue:
        eigenparam_matrix{i} = eigenvalues{i}(1).*first_eig;  
        % get the spike matrix:
        prop = (corrs{i}*N +2*meshgrid((means{i}+1)*N/2) +2*meshgrid((means{i}+1)*N/2)' - N)./4 ;
        prop = prop.*(ones(non,non)-eye(non))+(eye(non).*(((means{i}'+1)*N/2)*ones(1,non)));
        spike_matrix{i} = prop;
        % and the parameter matrix:
        prop = interactions{i};
        prop = prop.*(ones(non,non)-eye(non))+eye(non).*(fields{i}'*ones(1,non));
        parameter_matrix{i} = prop;
    end
    % file info and date of creation:
    Info_Fisher = ['This file was created by the script FIM_Decomposition on ',date];   
    % and save the results:
    filename_fisher = fileouts{j};
    save(filename_fisher,'Info_Fisher','filename_fisher','filename_spikes','filename_ising','fisher_matrix','eigenvalues','eigenvectors','spike_matrix','parameter_matrix','gini_FIM','eigenparam_matrix')
    disp(' ')
    disp(['Script finished, results written succesfully:  ',filename_fisher])
    disp(' ')
end
%---------END----------%
