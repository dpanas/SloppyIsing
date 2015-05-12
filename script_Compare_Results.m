% %-------------------------------%
% script:   Compare_Results
%           Script for some basic plots about the results, showing some
%           info about fit quality, sloppiness, sparsity, as well as
%           comparing similarity between groups in various measures.
%
% !!! This is a script, not a function - parameters need to be changed
% manually by the user upon each call !!!
%
% dependancy: SetUpFigure, PlotWithShade, SetPanelAppearance, PlotBoxplot, 
%             InformationMeasures, , R2Coefficient, PlotOverlappingHistograms 
%             ./functions/colormaps.mat, format of the result files 
%
% DAP Apr 2015
% %-------------------------------%

% first clear the workspace:
clear

% variables are here:
%-------------------%
filenames{1} = './results_fisher/chip136_0_fisher_8n_new_prob10_bin5_filt.mat';
filenames{2} = './results_fisher/chip136_1_fisher_8n_new_prob10_bin5_filt.mat';
filenames{3} = './results_fisher/chip136_2_fisher_8n_new_prob10_bin5_filt.mat';
filenames{4} = './results_fisher/chip136_4_fisher_8n_new_prob10_bin5_filt.mat';
filenames{5} = './results_fisher/chip136_5_fisher_8n_new_prob10_bin5_filt.mat';
filenames{6} = './results_fisher/chip136_6_fisher_8n_new_prob10_bin5_filt.mat';
%-------------------%

% (below no more variables, just code)
%-------------------%
close all
addpath('./functions/')
NOF = length(filenames);
eval(['load ./functions/colormaps'])
eval(['col = colorcell_',num2str(NOF),';'])
disp(' ')
disp(['Script comparing the ising and fisher results between groups of neurons from different recordings.'])
disp(' ')
% first setting up figures:
[fhandle_ising,subplots_ising] = SetUpFigure(4);     % figure for comparing ising results
[fhandle_sloppy,subplots_sloppy] = SetUpFigure(2);    % figure for comparing sloppiness and sparsity 
string_legend = ['legend('];
string_legend2 = ['legend('];
% loop for reading in, renaming / reformatting some data for plots, getting stats etc:
for i_file=1:NOF
    % loading the result files:
    filename_fisher = filenames{i_file};
    eval(['load ',filename_fisher])
    disp(' ')
    disp(['Loaded FIM decomposition results: ',' '' ',filename_fisher,' '' '])
    eval(['load ',filename_ising])
    disp([' and corresponding ising results: ',' '' ',filename_ising,' '' '])
    eval(['load ',filename_spikes])
    disp(['   and corresponding spikes file: ',' '' ',filename_spikes,' '' '])
    nos = length(key);
    non = length(key{1});
    % rename variables with overlapping names between files:
    eval(['spike_matrix_',num2str(i_file),' = spike_matrix;';])
    eval(['parameter_matrix_',num2str(i_file),' = parameter_matrix;';])
    eval(['eigparam_matrix_',num2str(i_file),' = eigenparam_matrix;';])
    eval(['eigenvalues_',num2str(i_file),' = eigenvalues;'])
    eval(['gini_',num2str(i_file),' = gini_FIM;'])
    eval(['fisher_matrix_',num2str(i_file),' = fisher_matrix;'])
    rates_m{i_file} = rates_model; 
    gini_coefficients{i_file} = gini_FIM;
    SJ{i_file} = zeros(1,nos);     % a cell for Shannon-Jensen divergences, for easy boxplotting
    eigval_mat = [];
    E_dat = zeros(1,nos);
    E_mod = zeros(1,nos);
    E_ind = zeros(1,nos);
    % a loop over the groups of neurons:
    for j=1:nos
        stats = corrs{j};
        stats = stats.*(ones(non,non)-eye(non)) + eye(non).*(means{j}'*ones(1,non));
        eval(['marginals_',num2str(i_file),'{',num2str(j),'} = stats;'])
        eigval_mat = [eigval_mat; eigenvalues{j}];
        % getting entropies for calculating multiinformation and SJ for
        % baseline recording between model and data:
        [e1,e2,~,~] = InformationMeasures(rates_data(j,:),rates_indep(j,:));
        E_dat = [E_dat e1];
        E_ind = [E_ind e2];
        [e1,e2,~,sj] = InformationMeasures(rates_data(j,:),rates_model(j,:));
        E_mod = [E_mod e2];
        if i_file==1
            SJ{i_file}(j) = sj;
            string_labels{i_file} = 'mod-dat';
        end
        % getting SJ divergence between models of further recordings and
        % baseline:
        if i_file>1
            [~,~,~,sj] = InformationMeasures(rates_m{i_file}(j,:),rates_m{1}(j,:));
            SJ{i_file}(j) = sj;
            string_labels{i_file} = ['mod',num2str(i_file),'vs1'];
        end
    end
    % if the first file was already read, comparison measures to it can be
    % computed:    
    if i_file>1
        % set up space:
        r2_param{i_file-1} = zeros(1,nos);
        r2_marg{i_file-1} = zeros(1,nos);
        r2_eig{i_file-1} = zeros(1,nos);
        r2_fim{i_file-1} = zeros(1,nos);
        ref = 1;
        % and again loop over groups of neurons, to compute R2 coefficient
        % of determination between values in a group from two different
        % time points:
        for j=1:nos
            eval(['r2_param{',num2str(i_file-1),'}(',num2str(j),') = R2Coefficient(parameter_matrix_',num2str(ref),'{',num2str(j),'}(:),parameter_matrix_',num2str(i_file),'{',num2str(j),'}(:));'])
            eval(['r2_marg{',num2str(i_file-1),'}(',num2str(j),') = R2Coefficient(abs(marginals_',num2str(ref),'{',num2str(j),'}(:)),abs(marginals_',num2str(i_file),'{',num2str(j),'}(:)));'])
            eval(['r2_eig{',num2str(i_file-1),'}(',num2str(j),') = R2Coefficient(eigparam_matrix_',num2str(ref),'{',num2str(j),'}(:),eigparam_matrix_',num2str(i_file),'{',num2str(j),'}(:));'])
            eval(['r2_fim{',num2str(i_file-1),'}(',num2str(j),') = R2Coefficient(fisher_matrix_',num2str(ref),'{',num2str(j),'}(:),fisher_matrix_',num2str(i_file),'{',num2str(j),'}(:));'])
        end
        edges = 0:0.1:1;
        hist_r2_param = histc(r2_param{i_file-1},edges);
        % plot a histogram of parameter similarity between recordings:
        figure(fhandle_ising)
        eval(subplots_ising{3})
        PlotOverlappingHistogram(hist_r2_param,edges,col{i_file});
        % and scatterplot of shannon-jensen divergence versus parameter
        % similarity
        eval(subplots_ising{4})
        plot(r2_param{i_file-1},SJ{i_file},'.','MarkerSize',6,'color',col{i_file})
        hold on
        string_legend2 = [string_legend2,'''Rec ',num2str(i_file),'vs1'','];
    end
    % plot multiinformation ratio:
    figure(fhandle_ising)
    eval(subplots_ising{1})
    plot((E_ind-E_dat),(E_ind-E_mod)./(E_ind-E_dat),'.','MarkerSize',6,'color',col{i_file})
    hold on
    figure(fhandle_sloppy)
    % plot the distributions of eigenvalues:
    eval(subplots_sloppy{1})
    PlotWithShade([],mean(eigval_mat,1),mean(eigval_mat,1)-std(eigval_mat,1),mean(eigval_mat,1)+std(eigval_mat,1),col{i_file},0.2);
    string_legend = [string_legend,'''Rec ',num2str(i_file),''','];
    clear prop stats params eigmatrix j k fisher_matrix rates_model means fields interactions corrs
end
% figure with ising results, fix up the appearance and legends etc:
string_legend = string_legend(1:end-1);
string_legend = [string_legend,',''location'',''Southeast'')'];
string_legend2 = string_legend2(1:end-1);
string_legend2 = [string_legend2,',''location'',''Northeast'')'];
figure(fhandle_ising)
eval(subplots_ising{1})
eval(string_legend)
SetPanelAppearance('Multiinformation ratio','Multiinformation',[],[],'Multiinf ratio',[0:0.2:1],[])
eval(subplots_ising{3})
SetPanelAppearance('Parameter similarity distributions','r2 coefficient',[0:0.2:1],[],'Histogram counts',[],[])
eval(subplots_ising{4})
SetPanelAppearance('Distribution vs parameter similarity','r2 coefficient',[0:0.2:1],[],'Shannon-Jensen divergence',[],[])
eval(string_legend2)
% plot a boxplot of shannon-jensen divergences:
eval(subplots_ising{2})
PlotBoxplot(SJ,string_labels,[]);
SetPanelAppearance('SJ div boxplots','    ',[],[],'Shannon-Jensen divergence',[],[])
% figure with fisher results, fix up the appearance and legends etc:
figure(fhandle_sloppy)
eval(subplots_sloppy{1})
SetPanelAppearance('Eigenvalue distributions','Rank',[0:10:50],[],'Eigenvalue',[],[])
eval(subplots_sloppy{2})
% plot a boxlpot of Gini coefficients:
PlotBoxplot(gini_coefficients,[],[]);
SetPanelAppearance('FIM gini coefficients','Recording',[0:1:NOF+1],[],'Gini coefficient',[0:.2:1],[])
clear rates_data rates_indep eigenvectors eigenvalues i time_total max_iter I Info_Fisher Info_Ising filename_fisher filename_ising filename_spikes
% figure with boxplots of r2 coefficients of marginals, parameters and
% FIMs:
figure
set(gcf,'position',[200 400 1200 400])
for i=2:NOF
    string_labels2{i-1} = string_labels{i};
    col2{i-1} = col{i};
end
subplot('position',[.09 .18 .23 .68])
PlotBoxplot(r2_marg,string_labels2,col2);
SetPanelAppearance('Marginals','',[],[],'r2 coefficient',[0:0.2:1],[])
subplot('position',[.41 .18 .23 .68])
PlotBoxplot(r2_param,string_labels2,col2);
SetPanelAppearance('Parameters','',[],[],'r2 coefficient',[0:0.2:1],[])
subplot('position',[.73 .18 .23 .68])
PlotBoxplot(r2_fim,string_labels2,col2);
SetPanelAppearance('FIM entries','',[],[],'r2 coefficient',[0:0.2:1],[])
%---------END----------%